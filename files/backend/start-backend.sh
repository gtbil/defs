#!/usr/bin/env bash
# Started by the Electron main process with these env vars set:
#   BP3_BUNDLE  path to the conda-packed environment (read-only, in resources)
#   BP3_PSGI    path to app.psgi
#   BP3_WEB     writable copy of the web tree (cgi-bin + htdocs)
#   BP3_PORT    localhost port to bind
set -euo pipefail

: "${BP3_BUNDLE:?}"; : "${BP3_PSGI:?}"; : "${BP3_WEB:?}"; : "${BP3_PORT:?}"

# --- One-time relocation of the packed conda env --------------------------
# conda-pack requires running conda-unpack once after the env is extracted at
# a new prefix. We copy the bundle into a writable, stable per-user location
# and unpack there (the bundle in resources is read-only).
ENV_DIR="${BP3_WEB%/web}/env"
if [[ ! -f "$ENV_DIR/.unpacked" ]]; then
  # Start clean: a previously failed unpack can leave a half-rewritten tree
  # whose read-only files would break a plain re-copy.
  rm -rf "$ENV_DIR"
  mkdir -p "$ENV_DIR"
  cp -a "$BP3_BUNDLE/." "$ENV_DIR/"
  # cpanm installs some .pm files read-only (mode 0444); conda-unpack rewrites
  # the build-prefix placeholder in place (opens files rb+), so the whole tree
  # must be user-writable first.
  chmod -R u+w "$ENV_DIR"
  # conda-unpack rewrites absolute prefixes baked into the packed env.
  if [[ -x "$ENV_DIR/bin/conda-unpack" ]]; then
    "$ENV_DIR/bin/conda-unpack"
  fi
  touch "$ENV_DIR/.unpacked"
fi

export PATH="$ENV_DIR/bin:$PATH"
# Help the loader find bundled shared libs (libgd, libpng, freetype, ...).
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
# GD font lookups
export GDFONTPATH="$ENV_DIR/fonts:${GDFONTPATH:-}"
# Let the CGI scripts find the bundled *.pm modules that sit beside them.
export PERL5LIB="$BP3_WEB/cgi-bin:${PERL5LIB:-}"

# Point the (now absolute) primer3 reference at the bundled binary.
P3="$ENV_DIR/bin/primer3_core"; [[ -x "$P3" ]] || P3="$ENV_DIR/bin/primer3"
if [[ -x "$P3" ]]; then
  ln -sf "$P3" "$BP3_WEB/cgi-bin/primer3_core"
else
  echo "WARNING: primer3 binary not found in bundle ($ENV_DIR/bin)" >&2
fi

# --- Bake runtime paths into the CGI copies -------------------------------
# web/cgi-bin/*.cgi carry placeholder tokens (__BP3_CGI__ etc.); we substitute the
# real values into the WRITABLE copy each launch (idempotent because we re-copy
# from the pristine tree only on first run).
substitute_token () { # $1 token  $2 value
  if grep -rlq "$1" "$BP3_WEB" 2>/dev/null; then
    grep -rl "$1" "$BP3_WEB" | while read -r f; do
      sed -i "s#$1#$2#g" "$f"
    done
  fi
}
substitute_token "__BP3_PORT__"   "$BP3_PORT"
substitute_token "__BP3_CGI__"    "$BP3_WEB/cgi-bin"
substitute_token "__BP3_HTDOCS__" "$BP3_WEB/htdocs"

# --- Serve ---------------------------------------------------------------
# Single-process is fine for a desktop app; bump --workers if you use Starman.
exec plackup \
  --host 127.0.0.1 \
  --port "$BP3_PORT" \
  "$BP3_PSGI"
