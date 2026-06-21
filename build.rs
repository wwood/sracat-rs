use std::env;

fn main() {
    // The ncbi-vdb headers and shared library come from the pixi environment.
    // CONDA_PREFIX is set automatically inside `pixi run`.
    let prefix = env::var("CONDA_PREFIX")
        .expect("CONDA_PREFIX not set; build inside `pixi run` (e.g. `pixi run build`)");

    // Compile the C shim against the ncbi-vdb headers.
    cc::Build::new()
        .file("shim.c")
        .include(format!("{prefix}/include"))
        .compile("sracatshim");

    // Link the ncbi-vdb shared library, and bake an rpath so the resulting
    // binary finds it without LD_LIBRARY_PATH.
    println!("cargo:rustc-link-search=native={prefix}/lib");
    println!("cargo:rustc-link-lib=dylib=ncbi-vdb");
    println!("cargo:rustc-link-arg=-Wl,-rpath,{prefix}/lib");

    println!("cargo:rerun-if-changed=shim.c");
    println!("cargo:rerun-if-changed=shim.h");
    println!("cargo:rerun-if-env-changed=CONDA_PREFIX");
}
