#!/usr/bin/env bash
# Publish the current state of python/toolsrtm and python/scopeinpython to
# their own standalone GitHub repos (ToolsRTMinPython, scopeinpython).
#
# Run this from anywhere, any time you want the split repos to reflect
# RTM-Suite's current state -- there is NO automatic sync, this script IS
# the sync step.
#
# What it does, per package: copies the package folder to a scratch
# directory (skipping build artifacts), creates a FRESH single commit
# there, and force-pushes it as `main` to the package's own repo. Each
# publish therefore REPLACES that repo's entire history with one new
# commit -- intentional, so nothing needs merging/reconciling, but it
# means: don't make edits directly in ToolsRTMinPython/scopeinpython on
# GitHub (a UI edit, a PR merge, ...) expecting them to survive -- the
# next sync overwrites them. Always edit inside RTM-Suite (python/toolsrtm,
# python/scopeinpython) and re-run this script to publish.
#
# Usage:
#   scripts/sync_python_repos.sh                  # sync both packages
#   scripts/sync_python_repos.sh toolsrtm          # just one
#   scripts/sync_python_repos.sh scopeinpython
#
# Requires: an SSH remote you can push to (git@github.com:CCGCAM/...).

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRATCH_DIR="${TMPDIR:-/tmp}/rtm_suite_python_sync"

declare -A REMOTES=(
  [toolsrtm]="git@github.com:CCGCAM/ToolsRTMinPython.git"
  [scopeinpython]="git@github.com:CCGCAM/scopeinpython.git"
)
declare -A COMMIT_MSGS=(
  [toolsrtm]="toolsrtm: Python port of ToolsRTM (leaf/canopy RTMs, sensor convolution, spectral indices)"
  [scopeinpython]="scopeinpython: Python port of SCOPEinR (soil, canopy BRDF, biochemistry, energy balance, fluorescence)"
)

sync_one() {
  local pkg="$1"
  local src="$REPO_ROOT/python/$pkg"
  local dst="$SCRATCH_DIR/$pkg"
  local remote="${REMOTES[$pkg]}"

  if [ ! -d "$src" ]; then
    echo "ERROR: $src not found -- skipping $pkg" >&2
    return 1
  fi

  echo "== $pkg -> $remote =="
  rm -rf "$dst"
  mkdir -p "$dst"
  cp -r "$src/." "$dst/"

  # Strip local build artifacts -- never want these in the published repo.
  find "$dst" -type d \( -name "__pycache__" -o -name "*.egg-info" -o -name ".pytest_cache" \) -exec rm -rf {} + 2>/dev/null || true

  printf '__pycache__/\n*.egg-info/\n.pytest_cache/\n*.pyc\nbuild/\ndist/\n' > "$dst/.gitignore"

  ( cd "$dst"
    git init -q
    git add -A
    git commit -q -m "${COMMIT_MSGS[$pkg]}"
    git branch -M main
    git remote add origin "$remote"
    git push --force -u origin main
  )
  echo "== $pkg published =="
  echo
}

if [ "$#" -eq 0 ]; then
  targets=(toolsrtm scopeinpython)
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
