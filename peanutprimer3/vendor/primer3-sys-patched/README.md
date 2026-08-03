[![Build](https://github.com/fg-labs/primer3-rs/actions/workflows/check.yml/badge.svg)](https://github.com/fg-labs/primer3-rs/actions/workflows/check.yml)
[![Crates.io](https://img.shields.io/crates/v/primer3)](https://crates.io/crates/primer3)
[![Documentation](https://img.shields.io/docsrs/primer3)](https://docs.rs/primer3)
[![License](http://img.shields.io/badge/license-GPL--2.0--or--later-blue.svg)](https://github.com/fg-labs/primer3-rs/blob/main/LICENSE)

# primer3-rs

Safe, idiomatic Rust bindings to the [primer3](https://primer3.org/) C library for PCR primer design and oligonucleotide thermodynamic calculations.

<p>
<a href="https://fulcrumgenomics.com">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset=".github/logos/fulcrumgenomics-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset=".github/logos/fulcrumgenomics-light.svg">
  <img alt="Fulcrum Genomics" src=".github/logos/fulcrumgenomics-light.svg" height="100">
</picture>
</a>
</p>

<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img src="https://img.shields.io/badge/Email_us-brightgreen.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>
<a href="https://www.fulcrumgenomics.com"><img src="https://img.shields.io/badge/Visit_Us-blue.svg?&style=for-the-badge&logo=wordpress&logoColor=white"/></a>

> [!WARNING]
> **primer3-rs is licensed under GPL-2.0-or-later.**
>
> This is unusual for our crates, and it is non-negotiable: the upstream
> [primer3](https://github.com/primer3-org/primer3) C library is GPL-2.0-or-later,
> and primer3-rs statically compiles and links it. Any application, library, or
> pipeline that depends on `primer3`, `primer3-sys`, or `primer3-tool` — whether
> via static or dynamic linking — forms a **combined work** subject to the terms
> of the GNU General Public License. In practice this means:
>
> - Distributing a binary that links primer3-rs requires releasing the complete
>   corresponding source of that binary under a GPL-compatible license.
> - You **cannot** combine primer3-rs with proprietary, closed-source, or
>   permissively-licensed-but-non-GPL-compatible code in a distributed product.
> - Internal/in-house use is unrestricted; the GPL's obligations trigger on
>   *distribution*, not on use.
>
> If these constraints are a problem for your use case, consider invoking the
> upstream `primer3_core` binary as a subprocess instead — the GPL does not
> reach across process boundaries in the same way. See the [License](#license)
> section below for details.
>
> **Not legal advice.** The above is our good-faith reading of the GPL and how
> it applies to primer3-rs. We are not lawyers, this is an opinion rather than
> a legal determination, and the line between "mere aggregation" and a
> "combined work" can be fact-specific. If the licensing implications matter
> to your project, please seek qualified legal counsel.

## Overview

primer3-rs provides three crates:

- **`primer3-sys`** -- Raw FFI bindings generated via `bindgen` against the primer3 C headers
- **`primer3`** -- Safe, idiomatic Rust API with builders, strong types, and proper error handling
- **`primer3-tool`** -- CLI binary exposing `tm`, `hairpin`, `homodimer`, `heterodimer`, and `design` subcommands

### Features

- Melting temperature (`Tm`) calculation with multiple methods (Breslauer, SantaLucia, SantaLucia 2004)
- Thermodynamic analysis: hairpin, homodimer, heterodimer, and 3' end stability
- Full primer design via `design_primers()` with typed settings and results
- Builder patterns for all complex configuration
- Thread-safe: thermodynamic functions are safe to call concurrently; `design_primers()` and `align()` are mutex-guarded
- Vendored C source compiled automatically -- no system dependencies required

## Installation

Supported platforms: Linux and macOS. Windows is not supported.

Add to your `Cargo.toml`:

```toml
[dependencies]
primer3 = "0.1"
```

The primer3 C library source is vendored and compiled automatically during build.
You will need a C/C++ compiler and `libclang` (for `bindgen`).

### Command-line tool

To install the `primer3-tool` CLI (which exposes `tm`, `hairpin`, `homodimer`, `heterodimer`, and `design` subcommands):

```console
cargo install primer3-tool
```

### System library (advanced)

To link against a system-installed `libprimer3` instead of compiling from source:

```toml
[dependencies]
primer3 = { version = "0.1", default-features = false, features = ["system"] }
```

## Quick Start

### Melting Temperature

```rust,no_run
use primer3::calc_tm;

let tm = calc_tm("GTAAAACGACGGCCAGT").unwrap();
println!("Tm = {tm:.1} C");
```

### Thermodynamic Analysis

```rust,no_run
use primer3::{calc_hairpin, calc_heterodimer};

let hairpin = calc_hairpin("CCCCCATCCGATCAGGGGG").unwrap();
println!("Hairpin Tm = {:.1} C, dG = {:.0} cal/mol", hairpin.tm(), hairpin.dg());

let dimer = calc_heterodimer("AAAAAAAAAA", "TTTTTTTTTT").unwrap();
println!("Heterodimer Tm = {:.1} C", dimer.tm());
```

### Primer Design

```rust,no_run
use primer3::{design_primers, SequenceArgs, PrimerSettings};

let seq_args = SequenceArgs::builder()
    .sequence("ATCGATCGATCGATCGATCG...")
    .target(100, 200)
    .build()
    .unwrap();

let settings = PrimerSettings::builder()
    .primer_opt_tm(60.0)
    .primer_min_tm(57.0)
    .primer_max_tm(63.0)
    .product_size_range(75, 150)
    .build()
    .unwrap();

let result = design_primers(&seq_args, &settings, None, None).unwrap();
for pair in result.pairs() {
    let left = pair.left();
    let right = pair.right();
    println!(
        "Left: {} (Tm={:.1}), Right: {} (Tm={:.1}), Penalty={:.2}",
        left.sequence(), left.tm(),
        right.sequence(), right.tm(),
        pair.penalty(),
    );
}
```

## Development

### Prerequisites

- Rust (see `rust-toolchain.toml` for the pinned version)
- C/C++ compiler (for compiling the vendored primer3 C library)
- `libclang` (for `bindgen` to generate FFI bindings)

### Getting Started

```console
git clone --recurse-submodules https://github.com/fg-labs/primer3-rs.git
cd primer3-rs
cargo build
cargo test
```

### Git Hooks

Install pre-commit hooks to run formatting and linting checks before each commit:

```console
./scripts/install-hooks.sh
```

### CI Commands

```console
cargo ci-fmt    # Check formatting
cargo ci-lint   # Run clippy with strict settings
cargo ci-test   # Run all tests
```

## License

primer3-rs is licensed under [GPL-2.0-or-later](LICENSE). All three crates in
this workspace (`primer3-sys`, `primer3`, `primer3-tool`) carry the same
license. See the warning at the top of this README for the practical
implications.

### Why GPL?

The upstream [primer3](https://github.com/primer3-org/primer3) C library is
licensed under GPL-2.0-or-later and its copyright holders have declined
requests to relicense it under a more permissive license. Because primer3-rs
statically compiles and links primer3's C source into the resulting artifacts,
the entire workspace is a derivative work and must be distributed under a
GPL-compatible license. We chose GPL-2.0-**or-later** (rather than GPL-2.0
only) to give downstream users the option of combining primer3-rs with
GPL-3.0-compatible code, including Apache-2.0–licensed code, at their
discretion.

### Escape hatches

If the GPL's copyleft obligations are incompatible with your project, the only
clean path is to not link against primer3-rs at all. Instead, invoke the
upstream `primer3_core` binary as a subprocess and communicate via the Boulder
I/O format on stdin/stdout. Running a GPL program as a separate process does
not impose GPL terms on the calling program, provided there is no intimate
linking (shared address space, shared data structures, etc.).

## Acknowledgments

- [primer3](https://primer3.org/) -- the underlying C library for primer design
- [primer3-py](https://github.com/libnano/primer3-py) -- Python bindings that inspired this project's API surface
