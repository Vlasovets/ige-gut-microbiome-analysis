#!/bin/bash
#SBATCH --job-name=rev_monitor
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --output=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/monitor_%j.out
#SBATCH --error=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/slurm/logs/monitor_%j.err

# Watchdog for the revision compute batch.
#
# Runs as a SLURM job rather than a background shell ON PURPOSE: a background
# process dies with the login session, and the point of this is to keep working
# after the user logs off.
#
# WHAT IT DOES
#   Every POLL seconds it reads sacct for the tracked jobs and classifies every
#   array task. Tasks that ended in FAILED / TIMEOUT / OUT_OF_MEMORY /
#   NODE_FAIL / PREEMPTED are resubmitted for just those indices, with
#   escalated resources chosen by failure mode:
#     OUT_OF_MEMORY -> double the memory
#     TIMEOUT       -> double the wall time
#     other         -> same resources (transient node/filesystem faults)
#   Retries are capped at MAX_RETRY rounds per job so a deterministic bug
#   cannot cause an infinite resubmission loop.
#
# WHAT IT DELIBERATELY DOES NOT DO
#   If the ASV SMOKE TEST fails, the dependent array is left cancelled and NOT
#   resubmitted. A smoke failure means the configuration is wrong, and
#   resubmitting 45 tasks against a broken config wastes shared cluster time.
#   That case is reported loudly instead.
#
# Usage:
#   sbatch revision/slurm/99_monitor_and_retry.sh <jobid> [<jobid> ...]

set -uo pipefail

REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
STATE="${REPO}/revision/slurm/logs/monitor_state.tsv"
REPORT="${REPO}/revision/slurm/logs/monitor_report.md"
POLL=300
MAX_RETRY=2
cd "${REPO}"

TRACKED=("$@")
if [ ${#TRACKED[@]} -eq 0 ]; then echo "FATAL: no job ids given" >&2; exit 1; fi

declare -A RETRIES
: > "${STATE}"

log() { echo "[$(date '+%F %T')] $*"; }

# Resolve each job's role from its OWN JobName, not from its position in the
# argument list.
#
# ⚠ This previously assigned roles positionally - arg1 = smoke, arg2 = asv,
# arg3 = ext. Two invocations passed a different number of ids, so both arrays
# were labelled "smoke". Consequences, had anything failed: a job labelled
# "smoke" is deliberately never retried, so NEITHER array had retry coverage;
# and script_for() would have resubmitted extended-DA indices using the ASV
# script, writing into the wrong results directory. Nothing failed, so
# "No task ever needed a retry" reflected luck, not demonstrated capability.
job_name() { sacct -j "$1" --noheader -o JobName%32 2>/dev/null | head -1 | tr -d ' '; }

role_of() {
  case "$(job_name "$1")" in
    rev_da_asv_smoke)  echo "smoke" ;;
    rev_ml_smoke)      echo "smoke" ;;
    rev_nested_smoke)  echo "smoke" ;;
    rev_da_asv)        echo "asv" ;;
    rev_da_ext)        echo "ext" ;;
    rev_alpha)         echo "alpha" ;;
    rev_ml_primary)    echo "ml_primary" ;;
    rev_ml_xgb)        echo "ml_xgb" ;;
    rev_nested)        echo "nested" ;;
    rev_net_warm)      echo "net_warm" ;;
    *)                 echo "unknown" ;;
  esac
}

script_for() {
  case "$(role_of "$1")" in
    smoke)      echo "" ;;   # smoke failures are never auto-retried
    asv)        echo "revision/slurm/02c_consensus_da_asv.sh" ;;
    ext)        echo "revision/slurm/02d_extended_da.sh" ;;
    alpha)      echo "revision/slurm/05c_alpha_diversity.sh" ;;
    ml_primary) echo "revision/slurm/04b_predictive_ml_primary.sh" ;;
    ml_xgb)     echo "revision/slurm/04c_predictive_ml_xgb_backfill.sh" ;;
    nested)     echo "revision/slurm/07b_nested_models.sh" ;;
    net_warm)   echo "revision/slurm/06b_warm_caches.sh" ;;
    *)          echo "" ;;
  esac
}

SMOKE_ID=""; ASV_ID=""; EXT_ID=""
for j in "${TRACKED[@]}"; do
  case "$(role_of "$j")" in
    smoke) SMOKE_ID="$j" ;;
    asv)   ASV_ID="$j" ;;
    ext)   EXT_ID="$j" ;;
    alpha|ml_primary|ml_xgb|nested|net_warm) : ;;  # known roles; retried via script_for
    *) echo "WARNING: job $j has unrecognised name '$(job_name "$j")'; it will be" \
            "watched but never auto-retried (no script mapping)." ;;
  esac
done

log "=== monitor start ==="
log "tracked: ${TRACKED[*]}"
log "  smoke=${SMOKE_ID}  asv_array=${ASV_ID}  ext_array=${EXT_ID}"

terminal_states="COMPLETED|FAILED|TIMEOUT|OUT_OF_MEMORY|CANCELLED|NODE_FAIL|PREEMPTED|DEADLINE"

all_done() {
  for j in "${TRACKED[@]}"; do
    [ -z "$j" ] && continue
    # A job is still live if squeue knows about it.
    if squeue -j "$j" -h -o "%i" 2>/dev/null | grep -q .; then return 1; fi
  done
  return 0
}

# Return failed array indices for a job, one per line.
failed_indices() {
  sacct -j "$1" --noheader --parsable2 --format=JobID,State 2>/dev/null \
    | grep -Ev '\.(batch|extern|[0-9]+)\|' \
    | grep -E "\|(FAILED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|PREEMPTED)" \
    | sed 's/|.*//' \
    | sed 's/.*_//' \
    | grep -E '^[0-9]+$' | sort -un
}

dominant_failure() {
  sacct -j "$1" --noheader --parsable2 --format=State 2>/dev/null \
    | grep -Eo 'FAILED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|PREEMPTED' \
    | sort | uniq -c | sort -rn | head -1 | awk '{print $2}'
}

while true; do
  if all_done; then log "all tracked jobs have left the queue"; break; fi

  for j in "${TRACKED[@]}"; do
    [ -z "$j" ] && continue
    squeue -j "$j" -h -o "%i" 2>/dev/null | grep -q . && continue   # still running

    mapfile -t bad < <(failed_indices "$j")
    [ ${#bad[@]} -eq 0 ] && continue

    n="${RETRIES[$j]:-0}"
    if [ -n "${SMOKE_ID}" ] && [ "$j" = "${SMOKE_ID}" ]; then
      log "SMOKE TEST ${j} FAILED. Not auto-retrying, and NOT resubmitting the"
      log "  dependent array ${ASV_ID} - a smoke failure means the config is wrong."
      log "  Inspect: revision/slurm/logs/da_asv_smoke_${j}.err"
      printf '%s\t%s\t%s\n' "$j" "smoke_failed" "manual_intervention_required" >> "${STATE}"
      TRACKED=("${TRACKED[@]/$j}")
      continue
    fi
    if [ "$n" -ge "$MAX_RETRY" ]; then
      log "job ${j}: ${#bad[@]} task(s) still failing after ${n} retries - giving up"
      printf '%s\t%s\t%s\n' "$j" "exhausted" "${bad[*]}" >> "${STATE}"
      TRACKED=("${TRACKED[@]/$j}")
      continue
    fi

    scr="$(script_for "$j")"
    if [ -z "$scr" ] || [ ! -f "$scr" ]; then
      log "job ${j}: cannot map to a script; skipping retry"
      TRACKED=("${TRACKED[@]/$j}")
      continue
    fi

    mode="$(dominant_failure "$j")"
    extra=()
    case "$mode" in
      OUT_OF_MEMORY) cur=$(sacct -j "$j" --noheader -o ReqMem | head -1 | tr -d ' Gcn')
                     extra=(--mem="$(( ${cur:-32} * 2 ))G")
                     log "job ${j}: dominant failure OOM - doubling memory to $(( ${cur:-32} * 2 ))G" ;;
      TIMEOUT)       extra=(--time=20:00:00)
                     log "job ${j}: dominant failure TIMEOUT - raising wall time to 20h" ;;
      *)             log "job ${j}: dominant failure ${mode} - retrying with same resources" ;;
    esac

    idx=$(IFS=,; echo "${bad[*]}")
    log "job ${j}: resubmitting ${#bad[@]} failed task(s) [${idx}] from ${scr} (retry $((n+1)))"
    new=$(sbatch --parsable --array="${idx}" "${extra[@]}" "$scr" 2>&1)
    if [[ "$new" =~ ^[0-9]+$ ]]; then
      log "  -> new job ${new}"
      printf '%s\t%s\t%s\n' "$j" "retried_as_${new}" "${idx}" >> "${STATE}"
      RETRIES[$new]=$((n+1))
      TRACKED=("${TRACKED[@]/$j}")
      TRACKED+=("$new")
      [ "$j" = "${ASV_ID}" ] && ASV_ID="$new"
      [ "$j" = "${EXT_ID}" ] && EXT_ID="$new"
    else
      log "  -> RESUBMIT FAILED: ${new}"
      printf '%s\t%s\t%s\n' "$j" "resubmit_error" "${new}" >> "${STATE}"
      TRACKED=("${TRACKED[@]/$j}")
    fi
  done
  sleep "${POLL}"
done

# ── Final report ─────────────────────────────────────────────────────────
{
  echo "# Compute batch report"
  echo
  echo "Generated $(date '+%F %T') by monitor job ${SLURM_JOB_ID}."
  echo
  echo '## Output files produced'
  echo
  printf '%s\n' "- consensus DA, ASV rank:   $(ls -1 revision/results/da/da_long_asv_*.csv 2>/dev/null | wc -l) of 45 files"
  printf '%s\n' "- extended DA (module 02d):  $(ls -1 revision/results/da_extended/*.csv 2>/dev/null | wc -l) files"
  echo
  echo '### ASV-rank consensus DA — expected 45'
  ls -1 revision/results/da/da_long_asv_*.csv 2>/dev/null | wc -l
  echo
  echo '### Extended DA — expected 180'
  ls -1 revision/results/da_extended/daext_long_*.csv 2>/dev/null | wc -l
  echo
  echo '## Retry ledger'
  echo
  if [ -s "${STATE}" ]; then
    echo '```'; cat "${STATE}"; echo '```'
  else
    echo 'No task ever needed a retry.'
  fi
} > "${REPORT}"

log "=== monitor done - report at ${REPORT} ==="
cat "${REPORT}"
