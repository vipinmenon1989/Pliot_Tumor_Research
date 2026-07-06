# Conda Environment Comparison Report
**Before vs After Snakemake Installation in R_env**

## 1. Summary of Changes
- **Packages Added**: 81
- **Packages Removed**: 0
- **Packages Upgraded**: 2
- **Packages Downgraded**: 0

## 2. Critical Package Integrity Verification
The following table lists the status of critical R and single-cell packages:

| Package | Status | Version | Change details |
| --- | --- | --- | --- |
| R | Unchanged | 4.4.3 | None |
| Seurat | Unchanged | 5.3.0 | None |
| SeuratObject | Unchanged | 5.1.0 | None |
| Matrix | Unchanged | 1.7_3 | None |
| sctransform | Unchanged | 0.4.2 | None |
| SingleCellExperiment | Unchanged | 1.28.0 | None |
| SummarizedExperiment | Unchanged | 1.36.0 | None |
| BiocGenerics | Unchanged | 0.52.0 | None |
| future | Unchanged | 1.67.0 | None |
| future.apply | Unchanged | 1.20.0 | None |
| ggplot2 | Unchanged | 3.5.2 | None |
| patchwork | Unchanged | 1.3.1 | None |
| data.table | Unchanged | 1.17.8 | None |
| harmony | Unchanged | 1.2.3 | None |
| DoubletFinder | Unchanged | NOT_INSTALLED | None |
| scDblFinder | Unchanged | 1.20.2 | None |

## 3. Detailed Package Transactions

### Packages Added
| Package | Version | Build | Channel |
| --- | --- | --- | --- |
| amply | 0.1.7 | pyhd8ed1ab_0 | conda-forge |
| annotated-doc | 0.0.4 | pyhcf101f3_0 | conda-forge |
| annotated-types | 0.7.0 | pyhd8ed1ab_1 | conda-forge |
| argparse-dataclass | 2.0.0 | pyhd8ed1ab_1 | conda-forge |
| attrs | 26.1.0 | pyhcf101f3_0 | conda-forge |
| backports.zstd | 1.6.0 | py311h6b1f9c4_0 | conda-forge |
| brotli-python | 1.2.0 | py311h66f275b_1 | conda-forge |
| certifi | 2026.6.17 | pyhd8ed1ab_0 | conda-forge |
| charset-normalizer | 3.4.7 | pyhd8ed1ab_0 | conda-forge |
| coin-or-cbc | 2.10.13 | h4d16d09_0 | conda-forge |
| coin-or-cgl | 0.60.10 | hc46dffc_0 | conda-forge |
| coin-or-clp | 1.17.11 | hc03379b_0 | conda-forge |
| coin-or-osi | 0.108.12 | hf4fecb4_0 | conda-forge |
| coin-or-utils | 2.11.13 | hc93afbd_0 | conda-forge |
| colorama | 0.4.6 | pyhd8ed1ab_1 | conda-forge |
| coloredlogs | 15.0.1 | pyhd8ed1ab_4 | conda-forge |
| conda-inject | 1.3.2 | pyhd8ed1ab_0 | conda-forge |
| configargparse | 1.7.5 | pyhcf101f3_0 | conda-forge |
| connection_pool | 0.0.3 | pyhd3deb0d_0 | conda-forge |
| docutils | 0.22.4 | pyhd8ed1ab_0 | conda-forge |
| dpath | 2.2.0 | pyha770c72_1 | conda-forge |
| eido | 0.2.5 | pyhd8ed1ab_0 | conda-forge |
| gitdb | 4.0.12 | pyhd8ed1ab_0 | conda-forge |
| gitpython | 3.1.50 | pyhd8ed1ab_0 | conda-forge |
| greenlet | 3.5.3 | py311hc665b79_0 | conda-forge |
| h2 | 4.3.0 | pyhcf101f3_0 | conda-forge |
| hpack | 4.2.0 | pyhd8ed1ab_0 | conda-forge |
| humanfriendly | 10.0 | pyh707e725_8 | conda-forge |
| hyperframe | 6.1.0 | pyhd8ed1ab_0 | conda-forge |
| idna | 3.18 | pyhcf101f3_0 | conda-forge |
| immutables | 0.21 | py311h49ec1c0_2 | conda-forge |
| jinja2 | 3.1.6 | pyhcf101f3_1 | conda-forge |
| jsonschema | 4.26.0 | pyhcf101f3_0 | conda-forge |
| jsonschema-specifications | 2025.9.1 | pyhcf101f3_0 | conda-forge |
| jupyter_core | 5.9.1 | pyhc90fa1f_0 | conda-forge |
| liblapacke | 3.9.0 | 32_he2f377e_openblas | conda-forge |
| logmuse | 0.3.0 | pyhcf101f3_0 | conda-forge |
| markdown-it-py | 4.2.0 | pyhd8ed1ab_0 | conda-forge |
| markupsafe | 3.0.3 | py311h3778330_1 | conda-forge |
| mdurl | 0.1.2 | pyhd8ed1ab_1 | conda-forge |
| nbformat | 5.10.4 | pyhd8ed1ab_1 | conda-forge |
| packaging | 25.0 | pyh29332c3_1 | conda-forge |
| pephubclient | 0.4.4 | pyhd8ed1ab_1 | conda-forge |
| peppy | 0.40.8 | pyhd8ed1ab_0 | conda-forge |
| platformdirs | 4.10.0 | pyhcf101f3_0 | conda-forge |
| psutil | 7.2.2 | py311haee01d2_0 | conda-forge |
| pulp | 2.8.0 | py311h77a8cca_3 | conda-forge |
| pydantic | 2.13.4 | pyhcf101f3_0 | conda-forge |
| pydantic-core | 2.46.4 | py311h902ca64_0 | conda-forge |
| pyparsing | 3.3.2 | pyhcf101f3_0 | conda-forge |
| pysocks | 1.7.1 | pyha55dd90_7 | conda-forge |
| python-fastjsonschema | 2.21.2 | pyhe01879c_0 | conda-forge |
| referencing | 0.37.0 | pyhcf101f3_0 | conda-forge |
| requests | 2.34.2 | pyhcf101f3_0 | conda-forge |
| rich | 15.0.0 | pyhcf101f3_0 | conda-forge |
| rpds-py | 2026.6.3 | py311h1baac5b_0 | conda-forge |
| shellingham | 1.5.4 | pyhd8ed1ab_2 | conda-forge |
| slack-sdk | 3.43.0 | pyhcf101f3_0 | conda-forge |
| slack_sdk | 3.43.0 | pyh9dfb50f_0 | conda-forge |
| smart_open | 7.7.1 | pyhcf101f3_0 | conda-forge |
| smmap | 5.0.3 | pyhcf101f3_1 | conda-forge |
| snakemake | 9.23.1 | hdfd78af_1 | bioconda |
| snakemake-interface-common | 1.23.0 | pyhdfd78af_1 | bioconda |
| snakemake-interface-executor-plugins | 9.4.0 | pyh84498cf_0 | bioconda |
| snakemake-interface-logger-plugins | 2.1.0 | pyhdfd78af_0 | bioconda |
| snakemake-interface-report-plugins | 1.3.0 | pyhd4c3c12_0 | bioconda |
| snakemake-interface-scheduler-plugins | 2.0.2 | pyhd4c3c12_0 | bioconda |
| snakemake-interface-storage-plugins | 4.4.1 | pyh84498cf_0 | bioconda |
| snakemake-minimal | 9.23.1 | pyhdfd78af_1 | bioconda |
| sqlalchemy | 2.0.51 | py311haee01d2_0 | conda-forge |
| sqlmodel | 0.0.37 | pyhcf101f3_0 | conda-forge |
| tabulate | 0.10.0 | pyhcf101f3_0 | conda-forge |
| tenacity | 9.1.4 | pyhcf101f3_0 | conda-forge |
| throttler | 1.2.2 | pyhd8ed1ab_0 | conda-forge |
| typer | 0.26.8 | pyhcf101f3_0 | conda-forge |
| typing-extensions | 4.15.0 | h396c80c_0 | conda-forge |
| typing-inspection | 0.4.2 | pyhcf101f3_2 | conda-forge |
| ubiquerg | 0.9.3 | pyhd8ed1ab_0 | conda-forge |
| urllib3 | 2.7.0 | pyhd8ed1ab_0 | conda-forge |
| wrapt | 1.17.3 | py311h49ec1c0_1 | conda-forge |
| yte | 1.9.4 | pyhd8ed1ab_0 | conda-forge |

### Packages Updated
| Package | Before Version | After Version | Channel |
| --- | --- | --- | --- |
| ca-certificates | 2026.1.4 | 2026.6.17 | conda-forge |
| openssl | 3.6.0 | 3.6.3 | conda-forge |
