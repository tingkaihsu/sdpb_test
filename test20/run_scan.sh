#!/usr/bin/env bash
# ============================================================================
# run_scan.sh -- loop over mAsq values:
#     Mathematica (PMP)  ->  pmp2sdp  ->  sdpb  ->  save out_out/z.txt
#
# usage (background, survives logout):
#     nohup bash ./run_scan.sh mAsq_list.txt > scan.log 2>&1 < /dev/null &
#     echo $!        # PID of the loop
#     disown
#
# input file: one mAsq per line, e.g.
#     1/5
#     0.15
#     # comments and blank lines are ignored
#
# per mass, in this directory:
#     z_mAsq_<tag>.txt     copy of out_out/z.txt   (only if SDPB found the optimum)
#     out_mAsq_<tag>.txt   copy of out_out/out.txt (terminateReason, objectives)
#     logs/mAsq_<tag>.log  full Mathematica / pmp2sdp / sdpb output
#     scan_summary.tsv     one line per mass: status, primalObjective, bound
# <tag> is the value as written in the input file, with "/" -> "_over_".
# Re-running the script skips every mass whose z_mAsq_<tag>.txt already exists.
#
# overridable from the environment:
#     MFILE=AAAA_test20.m  IMAGE=bootstrapcollaboration/sdpb:master  NPROC=24
# ============================================================================
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

INPUT="${1:-mAsq_list.txt}"
MFILE="${MFILE:-BBBB_test20.m}"
IMAGE="${IMAGE:-bootstrapcollaboration/sdpb:master}"
NPROC="${NPROC:-24}"

WORKDIR="$SCRIPT_DIR"
DATA_DIR="/data"
PMP_HOST="n_pmp.json"                  # written by the .m next to itself
N_PMP_JSON="${DATA_DIR}/n_pmp.json"
OUT_PREFIX="${DATA_DIR}/out"           # sdp dir: ./out ; results: ./out_out
LOGDIR="logs"
FAILDIR="failed"
SUMMARY="scan_summary.tsv"
CNAME="sdpb_scan_$$"                   # container names, so they can be killed

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

# ---------------------------------------------------------------- checks ----
[[ -f "$INPUT" ]] || { log "input file '$INPUT' not found"; exit 1; }
[[ -f "$MFILE" ]] || { log "Mathematica file '$MFILE' not found"; exit 1; }
for c in wolframscript docker flock setsid; do
  command -v "$c" >/dev/null 2>&1 || { log "'$c' not found in PATH"; exit 1; }
done

# one scan per directory: they would overwrite each other's n_pmp.json / out_out
exec 9> .run_scan.lock
flock -n 9 || { log "another run_scan.sh is already running in $SCRIPT_DIR"; exit 1; }

mkdir -p "$LOGDIR" "$FAILDIR"
[[ -f "$SUMMARY" ]] || printf 'mAsq\texact\tstatus\tprimalObjective\tbound(-primalObjective)\tfinished\n' > "$SUMMARY"

# ------------------------------------------------- clean stop on kill -------
CHILD=""
exec 8>&1                              # main log, also reachable inside redirected calls
cleanup() {
  [[ -n "$CHILD" ]] && kill -TERM -- "-$CHILD" 2>/dev/null   # whole process group
  docker rm -f "${CNAME}_pmp2sdp" "${CNAME}_sdpb" >/dev/null 2>&1
}
trap 'log "received signal, stopping" >&8; cleanup; exit 130' INT TERM

# run a command in the background and wait for it, so that a kill of this
# script is handled immediately (bash defers traps during foreground jobs)
run_bg() {
  setsid "$@" < /dev/null 8>&- 9>&- &   # own process group; no lock/log fds
  CHILD=$!
  wait "$CHILD"
  local rc=$?
  CHILD=""
  return $rc
}

# ------------------------------------------------------- one mass point -----
REASON=""
run_one() {   # $1 exact value passed to Mathematica, $2 per-mass log file
  local exact="$1" mlog="$2"

  # 0) remove everything a previous mass left behind (out*, owned by root,
  #    are removed from inside the container, as in run_sdpb.sh)
  rm -f "$PMP_HOST"
  run_bg docker run --rm -v "${WORKDIR}:${DATA_DIR}" "$IMAGE" \
      sh -lc "rm -rf ${OUT_PREFIX}.ck ${OUT_PREFIX} ${OUT_PREFIX}_out" >> "$mlog" 2>&1 \
    || { REASON="cleanup container failed"; return 1; }
  if [[ -e "$PMP_HOST" || -e out_out ]]; then REASON="old output still present"; return 1; fi

  # 1) PMP from Mathematica
  run_bg wolframscript -file "$MFILE" "$exact" >> "$mlog" 2>&1
  local rc=$?
  if (( rc != 0 )); then REASON="wolframscript exit $rc"; return 1; fi
  [[ -s "$PMP_HOST" ]] || { REASON="no $PMP_HOST produced"; return 1; }

  # 2) PMP -> SDP
  run_bg docker run --rm --name "${CNAME}_pmp2sdp" --shm-size=16384m \
      -v "${WORKDIR}:${DATA_DIR}" "$IMAGE" \
      pmp2sdp 2048 -i "$N_PMP_JSON" -o "$OUT_PREFIX" >> "$mlog" 2>&1 \
    || { REASON="pmp2sdp failed"; return 1; }

  # 3) SDPB (same options as run_sdpb.sh)
  run_bg docker run --rm --name "${CNAME}_sdpb" --shm-size=16384m \
      -v "${WORKDIR}:${DATA_DIR}" "$IMAGE" \
      mpirun --allow-run-as-root -n "$NPROC" sdpb \
        --writeSolution="x,y,z,X,Y" \
        --maxComplementarity=1e1000 \
        --dualityGapThreshold=1e-30 \
        --stepLengthReduction=0.7 \
        --primalErrorThreshold=1e-30 \
        --dualErrorThreshold=1e-30 \
        --precision=2048 \
        --maxIterations=50000 \
        --noFinalCheckpoint \
        -s "$OUT_PREFIX" >> "$mlog" 2>&1 \
    || { REASON="sdpb failed"; return 1; }

  [[ -s out_out/out.txt ]] || { REASON="no out_out/out.txt"; return 1; }
  [[ -s out_out/z.txt   ]] || { REASON="no out_out/z.txt";   return 1; }
  return 0
}

# ------------------------------------------------------------- main loop ----
mapfile -t LINES < "$INPUT"            # read the whole list once
log "scan started: $INPUT (${#LINES[@]} lines), $MFILE, PID $$"

for raw in ${LINES[@]+"${LINES[@]}"}; do
  v="${raw%%#*}"                       # drop comments
  v="${v//[[:space:]]/}"               # drop spaces, tabs and \r
  [[ -z "$v" ]] && continue

  # exact form for Mathematica: decimals become fractions (0.15 -> 15/100)
  if   [[ "$v" =~ ^[0-9]+/[0-9]+$ ]]; then exact="$v"
  elif [[ "$v" =~ ^[0-9]+$ ]];        then exact="$v"
  elif [[ "$v" =~ ^([0-9]*)\.([0-9]+)$ ]]; then
    num="${BASH_REMATCH[1]}${BASH_REMATCH[2]}"
    num="${num#"${num%%[!0]*}"}"; [[ -z "$num" ]] && num=0      # strip leading zeros
    exact="${num}/1$(printf '%0*d' "${#BASH_REMATCH[2]}" 0)"
  else
    log "skip '$raw': not an integer, decimal or fraction"
    printf '%s\t-\tinvalid input\t-\t-\t%s\n' "$v" "$(date '+%F %T')" >> "$SUMMARY"
    continue
  fi

  tag="${v//\//_over_}"
  result="z_mAsq_${tag}.txt"
  mlog="$LOGDIR/mAsq_${tag}.log"
  if [[ -s "$result" ]]; then log "skip mAsq=$v: $result exists"; continue; fi

  log "mAsq=$v (exact $exact): start"
  : > "$mlog"
  t0=$SECONDS

  if run_one "$exact" "$mlog"; then
    reason=$(sed -n 's/^ *terminateReason *= *"\(.*\)";.*/\1/p' out_out/out.txt)
    pobj=$(sed -n 's/^ *primalObjective *= *\([^;]*\);.*/\1/p' out_out/out.txt)
    if [[ "$pobj" == -* ]]; then bound="${pobj#-}"; else bound="-$pobj"; fi
    if [[ "$reason" == "found primal-dual optimal solution" ]]; then
      cp out_out/out.txt "out_mAsq_${tag}.txt"
      cp out_out/z.txt "${result}.part" && mv "${result}.part" "$result"   # atomic
      status="ok"
    else
      mkdir -p "$FAILDIR/$tag"
      cp out_out/out.txt out_out/z.txt "$FAILDIR/$tag/" 2>/dev/null
      status="not optimal: $reason"
    fi
  else
    pobj="-"; bound="-"; status="error: $REASON"
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$v" "$exact" "$status" "$pobj" "$bound" \
    "$(date '+%F %T')" >> "$SUMMARY"
  log "mAsq=$v: $status (bound $bound) in $((SECONDS - t0)) s"
done

log "scan finished"