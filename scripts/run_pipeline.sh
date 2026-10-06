#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data results/robustness figures/supplemental

Rscript scripts/export_author_metadata.R
python scripts/preprocess_qc.py
python scripts/cluster_nuclei.py
python scripts/annotate_clusters.py
python scripts/pseudobulk.py

Rscript scripts/differential_abundance.R
Rscript scripts/sex_de_edger.R
SEXDE_MODE=authorlike Rscript scripts/sex_de_edger.R
python scripts/sex_de_pydeseq2.py

ls data/pseudobulk_meta_sub_*.csv | sed 's#data/pseudobulk_meta##; s#\.csv##' | \
  xargs -P 3 -I{} sh -c 'SEXDE_TAGS={} SEXDE_OUT={} Rscript scripts/sex_de_edger.R'

python scripts/make_figures.py
