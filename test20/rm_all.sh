#!/usr/bin/env bash
set -euo pipefail

IMAGE="bootstrapcollaboration/sdpb:master"
WORKDIR="$(pwd)"
DATA_DIR="/data"
N_PMP_JSON="${DATA_DIR}/n_pmp.json"
OUT_PREFIX="${DATA_DIR}/out"

# 0) Run PMP
# wolframscript -file test17-1.m
# 1) Remove old checkpoint file
docker run --rm \
  -v "${WORKDIR}:${DATA_DIR}" \
  "${IMAGE}" \
  sh -lc "rm -rf ${OUT_PREFIX}*"
