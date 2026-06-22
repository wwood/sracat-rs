use std::env;

fn main() {
    // The ncbi-vdb headers and libraries come from the pixi environment.
    // CONDA_PREFIX is set automatically inside `pixi run`.
    let prefix = env::var("CONDA_PREFIX")
        .expect("CONDA_PREFIX not set; build inside `pixi run` (e.g. `pixi run build`)");

    // Compile the C shim against the ncbi-vdb headers.
    cc::Build::new()
        .file("shim.c")
        .include(format!("{prefix}/include"))
        .compile("sracatshim");

    println!("cargo:rustc-link-search=native={prefix}/lib");

    // How to link ncbi-vdb:
    //   dylib  (default) - link the shared lib and bake an rpath into the pixi
    //                      env, so the binary runs without LD_LIBRARY_PATH but
    //                      only on a machine where that env still exists. Fast
    //                      to link; right for local dev / `pixi run build`.
    //   static           - link libncbi-vdb.a so the binary carries no
    //                      libncbi-vdb.so dependency and no conda-path rpath:
    //                      self-contained and relocatable, for releases
    //                      (cargo-dist etc.). ncbi-vdb is C++ internally and
    //                      uses zlib, so its transitive deps are pulled in.
    let link = env::var("SRACAT_VDB_LINK").unwrap_or_else(|_| "dylib".into());
    match link.as_str() {
        "static" => {
            println!("cargo:rustc-link-lib=static=ncbi-vdb");
            // ncbi-vdb is C++ internally: link the platform C++ runtime
            // (libc++ on macOS, libstdc++ elsewhere) plus zlib.
            let cxx = match env::var("CARGO_CFG_TARGET_OS").as_deref() {
                Ok("macos") => "c++",
                _ => "stdc++",
            };
            println!("cargo:rustc-link-lib=dylib={cxx}");
            println!("cargo:rustc-link-lib=dylib=z");
            // libm / libpthread / libdl / libc are linked by default.
        }
        "dylib" => {
            println!("cargo:rustc-link-lib=dylib=ncbi-vdb");
            println!("cargo:rustc-link-arg=-Wl,-rpath,{prefix}/lib");
        }
        other => panic!("SRACAT_VDB_LINK must be 'dylib' or 'static', got {other:?}"),
    }

    println!("cargo:rerun-if-changed=shim.c");
    println!("cargo:rerun-if-changed=shim.h");
    println!("cargo:rerun-if-env-changed=CONDA_PREFIX");
    println!("cargo:rerun-if-env-changed=SRACAT_VDB_LINK");
}
