# Phase 6 regulatory / pathway networks

Cached so Phase 6 is reproducible without a live network call, and so the exact
network used is inspectable rather than implied.

| file | rows | source | tracked in git |
| --- | ---: | --- | :--: |
| `collectri_omnipath_raw.tsv` | 64,515 | OmniPath REST, `datasets=collectri&genesymbols=yes&organisms=9606` | **no** (17 MB, re-downloadable) |
| `collectri_human.tsv` | 41,674 | processed from the raw dump: source/target gene symbols + `mor` | yes |
| `progeny_human_top500.tsv` | 7,000 | `decoupleR::get_progeny(organism = "human", top = 500)` | yes |

`mor` follows the CollecTRI convention: **−1 only for purely repressive edges**
(`is_inhibition` true and `is_stimulation` false); +1 otherwise, including the
4,776 edges annotated as both.

## Why the REST endpoint rather than the R wrapper

`decoupleR::get_collectri()` and `OmnipathR::collectri()` both **fail** in this
project's frozen stack — OmnipathR 3.14.0's `unnest_evidences()` errors on the
CollecTRI static table (`if (.keep) . else select(., -!!evs_col)`).

**OmnipathR was not upgraded.** It is a frozen Phase 3/4 dependency and moving it
to fix a convenience wrapper would risk the CCC results that depend on it. The
identical data was taken from OmniPath's documented REST endpoint instead, which
is what the wrapper calls anyway.

Checksums: `NETWORK_CHECKSUMS.md5`.
