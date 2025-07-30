#!/bin/bash

set -euo pipefail

MAX_RETRIES=10
RETRY=0

while [ $RETRY -lt $MAX_RETRIES ]; do
    echo ">>> Starte Snakemake Rerun #$((RETRY+1))"
    snakemake --jobs 100 --rerun-incomplete --rerun-triggers mtime -k "$@" && break
    RETRY=$((RETRY+1))
    echo ">>> Snakemake failed, new try ($RETRY/$MAX_RETRIES)..."
    sleep 30
done

if [ $RETRY -eq $MAX_RETRIES ]; then
    echo ">>> Snakemake failed after $MAX_RETRIES tries."
    exit 1
fi
