# BatchPrimer3 desktop app

A self-contained desktop wrapper for the legacy BatchPrimer3 CGI application
(You et al. 2008). One **Linux x86_64** build runs natively on Linux and on
Windows through **WSLg** — there is no separate Windows build.

## How it works

The old setup ran the Perl CGI scripts under Apache2 + mod_perl inside a
writable Singularity sandbox started with `sudo`. This replaces all of that:

- **No container, no runtime to install.** Dependencies (Perl, `primer3`,
  `GD`/`libgd`, Plack, …) are built once with bioconda and packed into a
  relocatable directory (`conda-pack`) that carries its own shared libraries.
  The end user installs nothing.
- **No Apache, no root.** The same `.cgi` scripts are served by `plackup`
  (`Plack::App::CGIBin`) on a localhost port, as an ordinary user process.
- **A real window.** An Electron shell spawns the server, waits for it, opens
  the UI, and shuts the server down on quit. Results are written to a normal
  writable per-user directory instead of `chmod 777 /var/www/html`.

## Build (once, on a machine with internet + micromamba)

```bash
npm install
npm run build:backend      # create env, verify Perl modules, conda-pack
```

`build:backend` runs `scripts/build-backend.sh`, which:
1. creates the conda env from `backend/environment.yml`,
2. verifies each required Perl module and `cpanm`-installs any the channel
   named differently,
3. packs everything into `backend/bundle/`.

`web/` (the patched cgi-bin + htdocs) is committed straight to the repo and
isn't touched by this build — see "Layout" below. Building requires no
network access beyond fetching conda packages/Perl modules; there's no
journal archive to download or extract anymore.

## Run / package

```bash
npm start          # dev run
npm run dist       # AppImage + unpacked dir (electron-builder)
```

On Windows, run the Linux build inside WSL (WSLg provides the display). WSL2
forwards `localhost` to Windows, but the window renders through WSLg directly
so that isn't even needed.

## Source-specific fixes already applied

These are baked directly into the committed `web/` tree and `backend/app.psgi`
— there's no separate patch step to re-run:

- **primer3 binary.** The scripts call `'./primer3_core'` (relative), which
  won't resolve under Plack. It's rewritten to an absolute path and
  `start-backend.sh` symlinks the bundled bioconda `primer3_core` into place.
- **writable dirs.** Results go to `htdocs/batch_primers` (created by the
  script) and per-user parameter files go to `cgi-bin/parameters` (must
  pre-exist); both ship pre-created in `web/`.
- **`use Thread;`.** Vestigial and unused, and the old `Thread.pm` was removed
  from modern Perl — commented out.
- **`Thread->new(...)` in `batchprimer3_results.cgi`.** The actual primer-
  picking work was launched via the same long-removed `Thread` module and
  then immediately `->join()`ed, so it bought no concurrency — replaced with
  a direct subroutine call.
- **In-process CGI compilation.** `Plack::App::CGIBin` defaults to compiling
  each script into a subroutine and reusing it across requests. These scripts
  rely on old-style "one process per request" CGI semantics (named subs
  closing over top-level lexicals, e.g. `input_screen()` in
  `batchprimer3.cgi`), which breaks under reuse — primer defaults and help
  URLs would render blank after the first request. `app.psgi` now passes
  `exec_cb => sub { 1 }` to force real fork+exec per request instead.
- **Shebang lines.** `#!/usr/bin/env perl -w` isn't parseable by modern `env`
  (no multi-word args without `-S`), and `#!/usr/bin/perl` pointed at the bare
  system Perl (no `CGI.pm`). Both matter now that scripts are actually
  exec'd rather than compiled in-process. All three `.cgi` entry points use
  `#!/usr/bin/env perl`, which resolves to the bundled Perl via `$PATH`.

## Known limitation

The "email me the results" path uses `|/usr/sbin/sendmail`. There's no MTA in
the bundle, so that one feature won't work; all results still render in the
browser and are saved under `batch_primers/` regardless.

## Layout

```
main.js                  Electron shell (spawn backend, window, cleanup)
loading.html             boot splash
backend/
  environment.yml        conda dependency set
  app.psgi               Plack routing (cgi-bin + htdocs)
  start-backend.sh       first-run unpack, token substitution, plackup
  bundle/                <- generated: relocatable conda env
scripts/
  build-backend.sh       creates the conda env, verifies Perl modules, packs
web/                     committed: patched cgi-bin + htdocs (source of truth)
```

`web/` was originally generated from the BatchPrimer3 authors' Springer
supplementary archive (You et al. 2008) by a one-time extract-and-patch
script, but that step regenerated the tree from scratch on every build,
silently discarding any fixes applied afterward. It's now committed directly
and edited in place; there's nothing left to fetch or re-extract.
