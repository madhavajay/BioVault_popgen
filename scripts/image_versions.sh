#!/usr/bin/env bash
# Shared container image defaults for local shell runners.

: "${BIOSYNTH_IMAGE:=ghcr.io/openmined/biosynth:0.1.27}"
: "${BIOVAULT_IMAGE:=ghcr.io/madhavajay/biovault-popgen:0.1.7}"
: "${BIOVAULT_FAST_IMAGE:=ghcr.io/madhavajay/biovault-popgen:0.1.7-fast}"
