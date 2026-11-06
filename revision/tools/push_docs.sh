#!/bin/bash
# Publish the Markdown briefing to a PRIVATE GitHub repo, so it can be fetched
# and read on a local machine.
#
# WHY A SEPARATE PRIVATE REPO
#   These documents state, in plain language, that the submitted manuscript's
#   headline DA result was unadjusted, that four of five TaxSEA supplementary
#   rows are mathematically impossible, and that the differential networks were
#   never tested. All of that will be disclosed to the editor in the cover
#   letter - but it should reach the editor BEFORE it reaches the open web.
#   A private repo keeps the fetch-from-anywhere convenience without publishing
#   the manuscript's problems ahead of the response letter.
#
#   Make it public later if you want - `gh repo edit --visibility public`.
#
# WHAT IS PUBLISHED
#   revision/../docs/*.md only. No code, no CSVs, no individual-level data.
#   The docs carry aggregate, per-taxon and per-analysis statistics.
#
# PREREQUISITE (one time, interactive - this script cannot do it for you):
#   gh auth login
#
# Usage:
#   bash revision/tools/push_docs.sh            # create or update the repo
#   REPO_NAME=my-other-name bash revision/tools/push_docs.sh

set -euo pipefail

REPO=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis
REPO_NAME="${REPO_NAME:-ige-revision-briefing}"
STAGE="/ictstr01/scratch/users/${USER}/${REPO_NAME}"

cd "${REPO}"

if ! gh auth status >/dev/null 2>&1; then
  cat >&2 <<'MSG'
ERROR: gh is not authenticated.

Run this once, interactively, then re-run this script:

    gh auth login

Choose: GitHub.com -> HTTPS -> authenticate with a browser (or paste a token).
This script deliberately does not read stored credentials or ask you for a
token.
MSG
  exit 1
fi

echo "=== regenerating docs from the current results ==="
set +u; source ~/miniconda3/etc/profile.d/conda.sh; conda activate slr-env; set -u
python revision/analysis/10_export_markdown.py

if [ ! -d docs ] || [ -z "$(ls -A docs 2>/dev/null)" ]; then
  echo "FATAL: docs/ is empty - nothing to publish" >&2; exit 1
fi

echo "=== staging at ${STAGE} ==="
mkdir -p "${STAGE}"
# --exclude '.git' is load-bearing: without it --delete removes the staging
# repo's own .git (it does not exist under docs/), so every run re-inits a
# fresh history and the push is rejected as unrelated. Cost me one failed push.
rsync -a --delete --exclude '.git' docs/ "${STAGE}/"

# Guard: refuse to publish anything that is not Markdown. Cheap insurance
# against a stray CSV or figure being swept in by a future edit.
# Git's own internals are not publication content, so they are exempt - but
# nothing else is. The exemption is a path prefix, not a general loosening.
stray() { find "${STAGE}" -type f -not -path "${STAGE}/.git/*" ! -name '*.md'; }
if stray | grep -q .; then
  echo "FATAL: non-Markdown files staged for publication:" >&2
  stray >&2
  echo "Refusing to push. Remove them or update this guard deliberately." >&2
  exit 1
fi

cd "${STAGE}"
[ -d .git ] || git init -q -b main

if gh repo view "${REPO_NAME}" >/dev/null 2>&1; then
  echo "=== repo exists - pushing update ==="
  URL="https://github.com/$(gh repo view "${REPO_NAME}" --json nameWithOwner -q .nameWithOwner).git"
  if git remote get-url origin >/dev/null 2>&1; then
    git remote set-url origin "${URL}"
  else
    git remote add origin "${URL}"
  fi
  # Rebuild local history on top of whatever the remote already has, so an
  # earlier divergence resolves without force-pushing over someone else's work.
  if git fetch -q origin main 2>/dev/null; then
    git reset -q --soft FETCH_HEAD 2>/dev/null || true
  fi
  git add -A
  if git diff --cached --quiet 2>/dev/null; then
    echo "No changes to publish."
  else
    git commit -q -m "Update revision briefing ($(date '+%Y-%m-%d %H:%M'))"
  fi
  git push -q origin HEAD:main
else
  echo "=== creating PRIVATE repo ${REPO_NAME} ==="
  git add -A
  git diff --cached --quiet 2>/dev/null || \
    git commit -q -m "Revision briefing ($(date '+%Y-%m-%d %H:%M'))"
  gh repo create "${REPO_NAME}" --private --source=. --remote=origin --push
fi

echo
echo "=== published ==="
gh repo view "${REPO_NAME}" --json nameWithOwner,visibility,url \
  -q '"  " + .nameWithOwner + "  [" + .visibility + "]\n  " + .url'
echo
echo "Fetch it locally with:"
echo "  gh repo clone $(gh repo view "${REPO_NAME}" --json nameWithOwner -q .nameWithOwner)"
