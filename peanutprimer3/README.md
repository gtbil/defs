# PeanutPrimer3

A native, cross-platform desktop tool for designing KASP/PACE allele-specific
markers -- a from-scratch reimplementation of BatchPrimer3's "Allele-specific
primers and allele-flanking primers" mode, which doesn't work in the legacy
tool. See `../files/` for the legacy BatchPrimer3 wrapper this replaces.

For each SNP or indel marker, PeanutPrimer3 designs and automatically pairs:

- **Two allele-specific primers** (one per allele), each anchored so its 3'
  end sits on the polymorphic site.
- **One common/flanking primer**, in the orientation opposite the
  allele-specific pair.

Unlike the legacy tool, the common primer and its matching allele-specific
primers are always grouped and correctly oriented in the output -- there's
never a need to manually check orientation and re-pair rows by hand.

## Input format

Sequences are FASTA-like; a `>` header is optional (records are split by
blank lines / header-like lines otherwise). Mark the variant one of three
ways:

- An embedded IUPAC ambiguity code (`R`, `Y`, `S`, `W`, `K`, `M` for
  biallelic; `V`, `H`, `D`, `B` for triallelic), e.g. `...ACGT**R**ACGT...`.
- Bracket notation for a SNP: `...ACGT[A/G]ACGT...`.
- Bracket notation for an indel (insertion or deletion, `-` = empty allele):
  `...ACGT[A/ATT]ACGT...` or `...ACGT[AA/-]ACGT...`.

At least 100bp of flanking sequence on each side of the variant is
recommended (the app warns, but doesn't hard-block, when a marker has less).

## Workspace layout

- `core/` (`peanutprimer3-core`): all domain logic (parsing, variant model,
  primer candidate generation/scoring, ARMS second-mismatch, common-primer
  design and orientation pairing, KASP tails, batch orchestration, CSV/TSV/
  Excel export). No UI dependency -- see its module docs for the algorithm.
- `app/` (`peanutprimer3-app`): the `egui`/`eframe` desktop UI.
- `vendor/primer3-sys-patched/`: a locally patched copy of the `primer3-sys`
  crate (see below).

Primer thermodynamics and the common/flanking primer pair are delegated to
the real [primer3](https://primer3.org/) C library via the
[`primer3`](https://github.com/fg-labs/primer3-rs) Rust crate (vendored, so
no system primer3 install is required); the allele-specific-primer design,
ARMS logic, and orientation pairing are peanutprimer3's own domain logic.

## Building

Requires a C compiler for the host target (primer3's C source is compiled
from source) and, on Linux, `libclang` for `bindgen`.

```sh
cargo build --release          # native (Linux)
cargo test --workspace         # run all tests
cargo run -p peanutprimer3-app # launch the desktop app
```

## Cross-compiling to Windows (from Linux/RHEL9)

No Docker/Podman/`cross` needed -- just mingw-w64 and the Rust target:

```sh
sudo dnf install mingw64-gcc mingw64-gcc-c++      # if not already installed
rustup target add x86_64-pc-windows-gnu

CC_x86_64_pc_windows_gnu=x86_64-w64-mingw32-gcc \
CXX_x86_64_pc_windows_gnu=x86_64-w64-mingw32-g++ \
AR_x86_64_pc_windows_gnu=x86_64-w64-mingw32-ar \
cargo build --release --target x86_64-pc-windows-gnu -p peanutprimer3-app
```

The resulting `.exe` is at
`target/x86_64-pc-windows-gnu/release/peanutprimer3-app.exe`.

### Why `vendor/primer3-sys-patched/`

Upstream `primer3-sys` vendors and compiles primer3's C source, including
`masker.c`, which `#include`s `<sys/mman.h>` and uses `getline()` --
neither exists under mingw. Both are only used by primer3's list-file-based
repeat-masking feature (which peanutprimer3 never invokes -- we don't pass a
repeat library), but the whole file still has to compile since
`choose_primers()` unconditionally calls into it. `vendor/primer3-sys-patched/`
is a copy of the crate with `build.rs` patched to swap the mmap/`getline()`
calls for portable `fopen`/`fread`-based equivalents -- a straight
functional replacement (not a Windows-only stub), wired in via
`[patch.crates-io]` in the workspace `Cargo.toml`. If upstream fixes this,
the patch (and this vendored copy) can be dropped.

The same patched `build.rs` also statically links `libstdc++` on Windows
(`cargo:rustc-link-lib=static=stdc++`, paired with `-static-libgcc` in
`.cargo/config.toml`): `libprimer3.cc` is compiled as C++, and a plain
dynamic link pulls in `libstdc++-6.dll` at runtime, which isn't present on
a stock Windows install or in Wine without the mingw runtime -- the app
would otherwise fail to launch with `Library libstdc++-6.dll ... not
found`. Verified with `objdump -p` that the resulting `.exe` depends only
on standard Windows system DLLs, and confirmed it starts under Wine.

## macOS

Not yet targeted. Nothing in this architecture (a pure-Rust GUI, no
webview) blocks adding it later -- the recommended path is a GitHub Actions
`macos-14` runner (sidesteps Apple SDK licensing for local cross-compiling
entirely).

## License

GPL-2.0-or-later, inherited from the primer3 C library this links against.
