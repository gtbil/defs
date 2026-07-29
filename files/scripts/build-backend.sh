#!/usr/bin/env bash
# Build the self-contained backend. Run once on a machine WITH internet.
# Requires only `micromamba` (or conda/mamba) available on PATH for the BUILD.
# The end user needs NOTHING installed — the output is fully relocatable.
#
#   1. create the conda env from environment.yml
#   2. verify every Perl module the CGI needs; cpanm-install any missing
#   3. conda-pack the env into backend/bundle/
#
# `web/` (the patched cgi-bin + htdocs) is committed directly to the repo and
# is NOT regenerated here — see the top-level README's "Layout" section. It
# used to be extracted from a vendored Springer supplementary archive on
# every build, but that regenerated web/ from scratch each time, silently
# discarding any fixes made directly to the CGI scripts. web/ is now the
# source of truth; edit it in place.
set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
ENV_NAME="batchprimer3"

MAMBA="${MAMBA:-micromamba}"
command -v "$MAMBA" >/dev/null || { echo "Need micromamba/conda on PATH (set \$MAMBA)"; exit 1; }

# 1. env  (create if missing; otherwise install into it so a re-run just adds
#          new packages and keeps any modules cpanm already installed)
echo ">> creating/updating conda env"
if "$MAMBA" env list 2>/dev/null | grep -qE "/envs/${ENV_NAME}\b|^${ENV_NAME}\b"; then
  "$MAMBA" install -y -n "$ENV_NAME" -f "$HERE/backend/environment.yml"
else
  "$MAMBA" create -y -n "$ENV_NAME" -f "$HERE/backend/environment.yml"
fi
ENV_PREFIX="$("$MAMBA" run -n "$ENV_NAME" printenv CONDA_PREFIX 2>/dev/null | tail -n1 || true)"
[[ -n "$ENV_PREFIX" && -d "$ENV_PREFIX" ]] || { echo "Could not resolve env prefix"; exit 1; }

run() { "$MAMBA" run -n "$ENV_NAME" "$@"; }

# 2. verify / backfill Perl modules (handles wrong bioconda package names)
echo ">> verifying Perl modules"
REQUIRED_MODULES=(
  CGI GD GD::Graph::bars GD::Graph::colour GD::Text
  Archive::Zip Email::Valid IPC::Open3 POSIX Socket
  Plack Plack::App::CGIBin CGI::Emulate::PSGI CGI::Compile
)
MISSING=()
for m in "${REQUIRED_MODULES[@]}"; do
  if ! run perl -M"$m" -e1 >/dev/null 2>&1; then
    echo "   missing: $m"
    MISSING+=("$m")
  fi
done
if (( ${#MISSING[@]} )); then
  echo ">> installing missing modules with cpanm"
  run cpanm --notest "${MISSING[@]}"
fi

# 3. pack
# conda-pack is a Python tool (not a Perl module). Install it into the env,
# then call the binary directly by prefix so we DON'T activate the env
# (conda-pack refuses to pack the currently-active environment).
echo ">> conda-pack"
"$MAMBA" install -y -n "$ENV_NAME" -c conda-forge conda-pack
rm -rf "$HERE/backend/bundle"
mkdir -p "$HERE/backend/bundle"
"$ENV_PREFIX/bin/conda-pack" -p "$ENV_PREFIX" \
  --output "$HERE/backend/bundle.tar.gz" --force
tar -xzf "$HERE/backend/bundle.tar.gz" -C "$HERE/backend/bundle"
rm -f "$HERE/backend/bundle.tar.gz"

echo
echo "Done. Backend packed to backend/bundle/"
echo "Next:  npm install  &&  npm start     (dev run)"
echo "  or:  npm run dist                    (build AppImage / dir)"
