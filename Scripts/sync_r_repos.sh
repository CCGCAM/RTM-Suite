#!/usr/bin/env bash
# Publish ToolsRTM and SCOPEinR from RTM-Suite to their standalone GitLab repos.
#
# Each publish creates a fresh one-commit repository and force-pushes it as
# `main`. This intentionally replaces the destination history. Make package
# changes in RTM-Suite, not in the standalone GitLab repositories.
#
# Usage:
#   Scripts/sync_r_repos.sh              # publish both packages
#   Scripts/sync_r_repos.sh toolsrtm     # publish only ToolsRTM
#   Scripts/sync_r_repos.sh scopeinr     # publish only SCOPEinR
#
# Requires SSH push access to git@gitlab.com:caminoccg/...

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRATCH_DIR="${TMPDIR:-/tmp}/rtm_suite_r_sync"

declare -A SOURCES=(
  [toolsrtm]="ToolsRTM"
  [scopeinr]="SCOPEinR"
)
declare -A REMOTES=(
  [toolsrtm]="git@gitlab.com:caminoccg/toolsrtm.git"
  [scopeinr]="git@gitlab.com:caminoccg/scopeinr.git"
)
declare -A COMMIT_MSGS=(
  [toolsrtm]="ToolsRTM: publish current RTM-Suite package"
  [scopeinr]="SCOPEinR: publish current RTM-Suite package"
)

sync_one() {
  local pkg="$1"
  local src="$REPO_ROOT/${SOURCES[$pkg]}"
  local dst="$SCRATCH_DIR/$pkg"
  local remote="${REMOTES[$pkg]}"

  if [ ! -f "$src/DESCRIPTION" ]; then
    echo "ERROR: R package DESCRIPTION not found at $src" >&2
    return 1
  fi

  echo "== ${SOURCES[$pkg]} -> $remote =="
  rm -rf "$dst"
  mkdir -p "$dst"
  cp -r "$src/." "$dst/"

  # Remove repository metadata and local R build/check artifacts.
  find "$dst" -type d \( \
    -name ".git" -o \
    -name ".Rcheck" -o \
    -name "revdep" -o \
    -name "check" \
  \) -exec rm -rf {} + 2>/dev/null || true
  find "$dst" -type f \( \
    -name ".Rhistory" -o \
    -name ".RData" -o \
    -name "*.Rproj.user" -o \
    -name "*.tar.gz" \
  \) -delete

  (
    cd "$dst"
    git init -q
    git add -A
    git commit -q -m "${COMMIT_MSGS[$pkg]}"
    git branch -M main
    git remote add origin "$remote"
    git push --force -u origin main
  )

  echo "== ${SOURCES[$pkg]} published =="
  echo
}

if [ "$#" -eq 0 ]; then
  targets=(toolsrtm scopeinr)
else
  targets=("$@")
fi

for pkg in "${targets[@]}"; do
  if [ -z "${REMOTES[$pkg]+x}" ]; then
    echo "ERROR: unknown package \"$pkg\" -- expected one of: ${!REMOTES[*]}" >&2
    exit 1
  fi
  sync_one "$pkg"
done
