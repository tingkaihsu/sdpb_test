#!/usr/bin/env bash
set -euo pipefail

IMAGE="bootstrapcollaboration/sdpb:master"
WORKDIR="$(pwd)"
DATA_DIR="/data"
N_PMP_JSON="${DATA_DIR}/n_pmp.json"
OUT_PREFIX="${DATA_DIR}/out"

# 1) Remove old checkpoint file
docker run --rm \
  -v "${WORKDIR}:${DATA_DIR}" \
  "${IMAGE}" \
  sh -lc "rm -rf ${OUT_PREFIX}.ck"

# 2) Convert PMP to SDP
docker run --shm-size=16384m \
  -v "${WORKDIR}:${DATA_DIR}" \
  "${IMAGE}" \
  pmp2sdp 2048 -i "${N_PMP_JSON}" -o "${OUT_PREFIX}"

# 3) Run SDPB
# NOTE: --procsPerNode was removed — it is obsolete in SDPB 3.x and the MPI
# process layout is determined automatically from the MPI environment.
docker run --shm-size=16384m \
  -v "${WORKDIR}:${DATA_DIR}" \
  "${IMAGE}" \
  mpirun --allow-run-as-root -n 24 sdpb \
    --writeSolution="x,y,z,X,Y" \
    --maxComplementarity=1e1000 \
    --dualityGapThreshold=1e-30 \
    --stepLengthReduction=0.7 \
    --primalErrorThreshold=1e-30 \
    --dualErrorThreshold=1e-30 \
    --precision=2048 \
    --maxIterations=50000 \
    --noFinalCheckpoint \
    -s "${OUT_PREFIX}"
