use std::io::Write;
use std::process::{Command, Output, Stdio};

const ASL: &str = env!("CARGO_BIN_EXE_asl");

pub(crate) fn interpret(source: &str) -> Output {
    let mut child = Command::new(ASL)
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(source.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}

pub(crate) fn example(name: &str) -> Output {
    Command::new(ASL)
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .arg(format!("examples/{name}.asl"))
        .output()
        .unwrap()
}

pub(crate) fn stdout(output: Output) -> String {
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

pub(crate) fn diagnostic(output: Output) -> String {
    assert!(!output.status.success(), "{output:?}");
    let error = String::from_utf8(output.stderr).unwrap();
    assert!(!error.is_empty());
    assert!(!error.contains("panicked"), "{error}");
    error
}
