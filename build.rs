

fn main() {
    cc::Build::new()
        .file("src/faddeeva/Faddeeva.c")
        .compile("Faddeeva");
    println!("cargo:rerun-if-changed=src/faddeeva/Faddeeva.c");
}
