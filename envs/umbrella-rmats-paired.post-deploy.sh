#!/usr/bin/env bash
# Snakemake activates this conda environment before invoking the hook.
# Preserve upstream authorship and install the package recommended by rMATS.
set -euo pipefail
: "${CONDA_PREFIX:?Snakemake must activate the paired rMATS environment}"
pairadise_tmp=$(mktemp -d)
trap 'rm -rf "$pairadise_tmp"' EXIT
git clone --quiet https://github.com/Xinglab/PAIRADISE.git "$pairadise_tmp/source"
git -C "$pairadise_tmp/source" checkout --quiet --detach eff57cb3cc3fa23675d07d6c80b63dc6d864278e
"$CONDA_PREFIX/bin/R" CMD INSTALL --library="$CONDA_PREFIX/lib/R/library" \
    "$pairadise_tmp/source/pairadise/src/pairadise_model"
"$CONDA_PREFIX/bin/Rscript" -e 'stopifnot("my.data" %in% names(formals(PAIRADISE::pairadise)))'
