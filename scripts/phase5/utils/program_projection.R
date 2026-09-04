# =============================================================================
# Project cells onto already-derived cNMF program spectra.
#
# Used to place cells that were deliberately EXCLUDED from program discovery
# (non-malignant and Ambiguous historical fibroblasts) into the program space
# without letting them influence the programs. Projection is descriptive: it
# never reclassifies a cell.
#
# Method: non-negative least squares by multiplicative update with W FIXED -
# exactly the H-step of Frobenius NMF. Deterministic given the initialisation,
# needs no new package, and cannot move the spectra.
#
#   X  genes x cells   TPM of the cells being projected
#   W  genes x K       consensus program spectra in TPM units
#   H  K x cells       solved for, >= 0
# =============================================================================
project_programs <- function(X, W, n_iter = 300L, eps = 1e-10) {
  stopifnot(nrow(X) == nrow(W))
  X <- as.matrix(X); W <- as.matrix(W)
  WtW <- crossprod(W)                       # K x K
  WtX <- crossprod(W, X)                    # K x cells
  H <- matrix(mean(X) / (mean(W) + eps), nrow = ncol(W), ncol = ncol(X))
  rownames(H) <- colnames(W); colnames(H) <- colnames(X)
  for (i in seq_len(n_iter)) H <- H * (WtX / (WtW %*% H + eps))
  H[!is.finite(H)] <- 0
  t(H)                                      # cells x K
}

# counts (genes x cells) -> TPM-like column-normalised matrix on a gene subset
tpm_on_genes <- function(counts, genes, scale_factor = 1e6) {
  g <- intersect(genes, rownames(counts))
  ls <- Matrix::colSums(counts)
  ls[ls == 0] <- 1
  m <- counts[g, , drop = FALSE]
  m <- as.matrix(m) %*% Matrix::Diagonal(x = scale_factor / ls)
  m <- as.matrix(m); dimnames(m) <- list(g, colnames(counts))
  # genes present in the spectra but absent here contribute zero
  out <- matrix(0, nrow = length(genes), ncol = ncol(counts),
                dimnames = list(genes, colnames(counts)))
  out[g, ] <- m
  out
}

relative_usage <- function(H) {
  s <- rowSums(H); s[s == 0] <- 1
  H / s
}
