

fn main() {
    cc::Build::new()
        .file("src/Faddeeva.c")
        .compile("Faddeeva");
    println!("cargo:rerun-if-changed=src/Faddeeva.c");
}
