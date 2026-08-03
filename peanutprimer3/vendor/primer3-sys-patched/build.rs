use std::env;
use std::fs;
use std::path::{Path, PathBuf};

/// Rewrites a vendored C source file into `OUT_DIR` after applying a set of
/// exact-string replacements, and returns the path to the patched copy.
///
/// Used to carry targeted thread-safety patches against upstream primer3
/// without mutating the vendored submodule tree. Each replacement is applied
/// exactly once and asserts the target string is present, so a stale patch
/// fails the build loudly as soon as upstream fixes the underlying issue.
#[cfg(feature = "vendored")]
fn patch_source(src: &Path, out_dir: &Path, replacements: &[(&str, &str)], label: &str) -> PathBuf {
    let filename = src.file_name().expect("source path has no file name");
    let dst = out_dir.join(filename);
    let mut contents =
        fs::read_to_string(src).unwrap_or_else(|e| panic!("failed to read {}: {e}", src.display()));
    for (needle, replacement) in replacements {
        assert!(
            contents.contains(needle),
            "{label}: patch target not found in {} — upstream may have fixed this; \
             remove the corresponding patch from build.rs.",
            src.display()
        );
        contents = contents.replacen(needle, replacement, 1);
    }
    fs::write(&dst, contents).unwrap_or_else(|e| panic!("failed to write {}: {e}", dst.display()));
    println!("cargo:rerun-if-changed={}", src.display());
    dst
}

#[allow(clippy::too_many_lines)]
fn main() {
    let src_dir = PathBuf::from("vendor/primer3/src");
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());

    // Thread-safety patch for oligotm.c: the Owczarzy salt correction branch
    // declares its polynomial coefficients as function-local `static`
    // variables, racing under concurrent oligotm() calls. Dropping `static`
    // is a pure win (values are overwritten on every call anyway). Pending
    // upstream at primer3-org/primer3#96.
    #[cfg(feature = "vendored")]
    let oligotm_src = patch_source(
        &src_dir.join("oligotm.c"),
        &out_dir,
        &[(
            "    static double a = 0,b = 0,c = 0,d = 0,e = 0,f = 0,g = 0;",
            "    double a = 0,b = 0,c = 0,d = 0,e = 0,f = 0,g = 0;",
        )],
        "oligotm.c thread-safety",
    );

    // NOTE: dpal.c's `_dpal_generic()` declares its ~40 MB scoring and
    // traceback matrices as function-local `static` arrays, which is a
    // concurrency bug. We tried promoting them to `static _Thread_local`
    // here, but on Linux with glibc the resulting 40 MB per-thread static
    // TLS block overflows Rust's default 2 MB test thread stack, so the fix
    // has to live at the Rust boundary instead — see the ALIGN_MUTEX in
    // primer3/src/alignment.rs. A proper upstream fix would heap-allocate
    // these matrices lazily via thread-local pointers.

    // Windows/mingw portability patch (peanutprimer3): masker.c
    // unconditionally #includes <sys/mman.h> and uses mmap()/munmap() to
    // load repeat-masking "list files", which doesn't exist under mingw.
    // `choose_primers()` in libprimer3.cc unconditionally calls into
    // masker.c (create_output_sequence/read_and_mask_sequence/
    // delete_input_sequence/delete_output_sequence), so masker.c can't
    // simply be dropped from the build -- it has to actually compile.
    // Since the mmap usage here is just "read a whole file into memory"
    // (no shared-memory or partial-mapping semantics are relied upon),
    // swapping it for a portable fopen/fread(+free instead of munmap)
    // implementation is a straight functional equivalent, not a
    // Windows-only stub -- the list-file-based masking feature (which
    // peanutprimer3 doesn't use; we never pass a repeat library) keeps
    // working identically on every platform. Upstream:
    // https://github.com/fg-labs/primer3-rs (report if not already fixed).
    #[cfg(feature = "vendored")]
    let masker_src = patch_source(
        &src_dir.join("masker.c"),
        &out_dir,
        &[
            ("#include <sys/mman.h>\n", ""),
            (
                "const char *\nmmap_by_filename (const char *filename, size_t *size)\n{\n  struct stat st;\n  int status, handle;\n  const char *data;\n\n  status = stat (filename, &st);\n  if (status < 0) {\n    return NULL;\n  }\n\n  handle = open (filename, O_RDONLY);\n  if (handle < 0) {\n    return NULL;\n  }\n\n  data = (const char *) mmap (NULL, st.st_size, PROT_READ, MAP_PRIVATE, handle, 0);\n  if (data == (const char *) -1) {\n    return NULL;\n  } else {\n    *size = (size_t)st.st_size;\n  }\n\n  close (handle);\n  return data;\n}",
                "const char *\nmmap_by_filename (const char *filename, size_t *size)\n{\n  FILE *f;\n  long file_size;\n  char *data;\n\n  f = fopen (filename, \"rb\");\n  if (!f) {\n    return NULL;\n  }\n  if (fseek (f, 0, SEEK_END) != 0) {\n    fclose (f);\n    return NULL;\n  }\n  file_size = ftell (f);\n  if (file_size < 0) {\n    fclose (f);\n    return NULL;\n  }\n  rewind (f);\n\n  data = (char *) malloc ((size_t) file_size);\n  if (!data) {\n    fclose (f);\n    return NULL;\n  }\n  if (file_size > 0 && fread (data, 1, (size_t) file_size, f) != (size_t) file_size) {\n    free (data);\n    fclose (f);\n    return NULL;\n  }\n  fclose (f);\n\n  *size = (size_t) file_size;\n  return data;\n}",
            ),
            (
                "if (fp[i]->pointer) munmap ((void *) fp[i]->pointer, fp[i]->size);",
                "if (fp[i]->pointer) free ((void *) fp[i]->pointer);",
            ),
            // getline() is a POSIX/GNU extension the mingw-w64 runtime
            // doesn't provide at all (not just an unmapped declaration --
            // the symbol itself isn't in its C library), used elsewhere in
            // this file to read the list-file format line by line. Same
            // portable-equivalent approach as above: a small drop-in
            // implementation with the same growing-buffer semantics,
            // scoped to Windows only so every other platform keeps using
            // its native libc getline() unchanged.
            (
                "#include \"libprimer3.h\"\n\nunsigned int glistmaker_code_match = 'G' << 24 | 'T' << 16 | '4' << 8 | 'C';",
                "#include \"libprimer3.h\"\n\n#ifdef _WIN32\nstatic long\npeanutprimer3_getline (char **lineptr, size_t *n, FILE *stream)\n{\n  size_t pos = 0;\n  int c;\n\n  if (lineptr == NULL || n == NULL || stream == NULL) {\n    return -1;\n  }\n  if (*lineptr == NULL || *n == 0) {\n    *n = 128;\n    *lineptr = (char *) malloc (*n);\n    if (*lineptr == NULL) {\n      return -1;\n    }\n  }\n  while ((c = fgetc (stream)) != EOF) {\n    if (pos + 1 >= *n) {\n      size_t new_size = *n * 2;\n      char *new_ptr = (char *) realloc (*lineptr, new_size);\n      if (new_ptr == NULL) {\n        return -1;\n      }\n      *n = new_size;\n      *lineptr = new_ptr;\n    }\n    (*lineptr)[pos++] = (char) c;\n    if (c == '\\n') {\n      break;\n    }\n  }\n  if (pos == 0 && c == EOF) {\n    return -1;\n  }\n  (*lineptr)[pos] = '\\0';\n  return (long) pos;\n}\n#define getline peanutprimer3_getline\n#endif\n\nunsigned int glistmaker_code_match = 'G' << 24 | 'T' << 16 | '4' << 8 | 'C';",
            ),
        ],
        "masker.c mingw portability (no sys/mman.h)",
    );

    #[cfg(feature = "vendored")]
    {
        // C sources (compiled as C)
        // MAX_PRIMER_LENGTH defaults to 36 in the C source but can be up to
        // THAL_MAX_ALIGN (60). Set it to 60 so oligotm/seqtm can handle
        // sequences up to 60 bp with the nearest-neighbor model.
        cc::Build::new()
            .files([
                src_dir.join("thal.c"),
                oligotm_src.clone(),
                src_dir.join("dpal.c"),
                src_dir.join("p3_seq_lib.c"),
                masker_src,
                src_dir.join("read_boulder.c"),
                src_dir.join("print_boulder.c"),
                src_dir.join("format_output.c"),
            ])
            .include(&src_dir)
            .define("MAX_PRIMER_LENGTH", "60")
            .opt_level(2)
            .warnings(false)
            // thal.c needs -ffloat-store for reproducible floating point
            .flag_if_supported("-ffloat-store")
            .compile("primer3_c");

        // C++ source (libprimer3.cc requires C++11)
        cc::Build::new()
            .file(src_dir.join("libprimer3.cc"))
            .include(&src_dir)
            .define("MAX_PRIMER_LENGTH", "60")
            .cpp(true)
            .std("c++11")
            .opt_level(2)
            .warnings(false)
            .compile("primer3_cc");

        // Link the C++ standard library
        let target = env::var("TARGET").unwrap();
        if target.contains("apple") || target.contains("freebsd") {
            println!("cargo:rustc-link-lib=c++");
        } else if target.contains("windows") {
            // Static link (peanutprimer3 patch): a plain "-lstdc++" resolves
            // to the shared libstdc++-6.dll under mingw, which isn't present
            // on a stock Windows install or in Wine without the mingw
            // runtime -- the app fails to launch with a missing-DLL error.
            // Pair with `-static-libgcc` in .cargo/config.toml for the
            // implicit libgcc dependency too. rustc's own `-L` search paths
            // don't include the mingw compiler's internal lib directory, so
            // ask the compiler itself where libstdc++.a lives rather than
            // hardcoding a distro-specific path.
            let cxx = env::var("CXX_x86_64_pc_windows_gnu")
                .or_else(|_| env::var("CXX"))
                .unwrap_or_else(|_| "x86_64-w64-mingw32-g++".to_string());
            if let Ok(output) = std::process::Command::new(&cxx).arg("-print-file-name=libstdc++.a").output() {
                let path = String::from_utf8_lossy(&output.stdout).trim().to_string();
                if let Some(dir) = Path::new(&path).parent() {
                    println!("cargo:rustc-link-search=native={}", dir.display());
                }
            }
            println!("cargo:rustc-link-lib=static=stdc++");
        } else {
            println!("cargo:rustc-link-lib=stdc++");
        }
    }

    #[cfg(feature = "system")]
    {
        println!("cargo:rustc-link-lib=primer3");
    }

    // Generate bindings
    let bindings = bindgen::Builder::default()
        .header("wrapper.h")
        .clang_arg(format!("-I{}", src_dir.display()))
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .allowlist_function("thal")
        .allowlist_function("set_thal_default_args")
        .allowlist_function("set_thal_oligo_default_args")
        .allowlist_function("thal_set_null_parameters")
        .allowlist_function("thal_load_parameters")
        .allowlist_function("thal_free_parameters")
        .allowlist_function("get_thermodynamic_values")
        .allowlist_function("destroy_thal_structures")
        .allowlist_function("oligotm")
        .allowlist_function("seqtm")
        .allowlist_function("long_seq_tm")
        .allowlist_function("oligodg")
        .allowlist_function("end_oligodg")
        .allowlist_function("symmetry")
        .allowlist_function("divalent_to_monovalent")
        .allowlist_function("p3_create_global_settings")
        .allowlist_function("p3_destroy_global_settings")
        .allowlist_function("create_seq_arg")
        .allowlist_function("destroy_seq_args")
        .allowlist_function("choose_primers")
        .allowlist_function("destroy_p3retval")
        .allowlist_function("p3_set_gs_.*")
        .allowlist_function("p3_set_sa_.*")
        .allowlist_function("p3_sa_add_to_.*")
        .allowlist_function("read_and_create_seq_lib")
        .allowlist_function("destroy_seq_lib")
        .allowlist_function("seq_lib_num_seq")
        .allowlist_function("create_empty_seq_lib")
        .allowlist_function("add_seq_and_rev_comp_to_seq_lib")
        .allowlist_type("thal_args")
        .allowlist_type("thal_results")
        .allowlist_type("thal_parameters")
        .allowlist_type("thal_alignment_type")
        .allowlist_type("thal_mode")
        .allowlist_type("tm_ret")
        .allowlist_type("tm_method_type")
        .allowlist_type("salt_correction_type")
        .allowlist_type("p3_global_settings")
        .allowlist_type("seq_args")
        .allowlist_type("p3retval")
        .allowlist_type("primer_rec")
        .allowlist_type("primer_pair")
        .allowlist_type("oligo_array")
        .allowlist_type("oligo_stats")
        .allowlist_type("pair_stats")
        .allowlist_type("seq_lib")
        .allowlist_type("args_for_one_oligo_or_primer")
        .allowlist_type("oligo_weights")
        .allowlist_type("pair_weights")
        .allowlist_type("task")
        // masker types
        .allowlist_type("masker_parameters")
        .allowlist_type("formula_parameters")
        .allowlist_type("input_sequence")
        .allowlist_type("output_sequence")
        .allowlist_type("masking_direction")
        .allowlist_var("THAL_MAX_ALIGN")
        .allowlist_var("THAL_MAX_SEQ")
        .allowlist_var("_INFINITY")
        .allowlist_var("ABSOLUTE_ZERO")
        .allowlist_var("MAX_LOOP")
        .allowlist_var("MIN_LOOP")
        .generate()
        .expect("Unable to generate bindings");

    bindings.write_to_file(out_dir.join("bindings.rs")).expect("Couldn't write bindings!");
}
