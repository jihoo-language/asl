mod adams;
mod bracket;
mod cobar;
mod cohomology;
mod comodule;
mod complex;
mod cyclic;
mod engine;
mod error;
mod linear;
mod map_bracket;
mod map_cochain;
mod map_defining;
mod map_ext;
mod map_may;
mod map_transfer;
mod may;
mod milnor;
mod model;
mod observations;
mod profile;
mod steenrod;
mod syntax;
mod transfer;

use std::env;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Read, Write};
use std::process::ExitCode;

use crate::engine::Engine;
use crate::error::{Error, Result, Span};
use crate::syntax::Parser;

const MAX_SOURCE_BYTES: usize = 512 * 1024;

const HELP: &str = "\
ASL, the Alexandrov Section Lambda interpreter

Usage: asl <PATH>
       asl -
       asl --help

Pass a UTF-8 ASL source file or use - to read standard input.
The show, card, betti, extdim, extandim, xextdim, mapextdim, maydim, mayandim, mapmaydim, and homotopy statements print to standard output.
A syntax error, invalid program, or failed assertion exits unsuccessfully.
";

fn main() -> ExitCode {
    let args: Vec<OsString> = env::args_os().skip(1).collect();
    match args.as_slice() {
        [arg] if arg == "-h" || arg == "--help" => {
            print!("{HELP}");
            ExitCode::SUCCESS
        }
        [arg] => {
            let source_name = if arg == "-" {
                "<stdin>".into()
            } else {
                arg.to_string_lossy()
            };
            match run(arg) {
                Ok(()) => ExitCode::SUCCESS,
                Err(error) => {
                    eprintln!("{source_name}{error}");
                    ExitCode::FAILURE
                }
            }
        }
        _ => {
            eprintln!("asl: expected one source path or '-'. See asl --help.");
            ExitCode::from(2)
        }
    }
}

fn run(path: &OsString) -> Result<()> {
    let source = if path == "-" {
        read_source(io::stdin().lock())?
    } else {
        let file = File::open(path)
            .map_err(|error| Error::plain(format!("cannot open source: {error}")))?;
        read_source(file)?
    };
    let mut parser = Parser::new(&source)?;
    let mut engine = Engine::new();
    let mut stdout = io::stdout().lock();
    while let Some(statement) = parser.statement()? {
        if let Some(line) = engine.execute(statement)? {
            writeln!(stdout, "{line}").map_err(output_error)?;
        }
    }
    stdout.flush().map_err(output_error)
}

fn output_error(error: io::Error) -> Error {
    Error::plain(format!("cannot write output: {error}"))
}

fn read_source(reader: impl Read) -> Result<String> {
    let mut bytes = Vec::new();
    reader
        .take((MAX_SOURCE_BYTES + 1) as u64)
        .read_to_end(&mut bytes)
        .map_err(|error| Error::plain(format!("cannot read source: {error}")))?;
    if bytes.len() > MAX_SOURCE_BYTES {
        return Err(Error::plain("source exceeds the size limit of 512 KiB"));
    }
    String::from_utf8(bytes).map_err(|error| {
        let prefix = &error.as_bytes()[..error.utf8_error().valid_up_to()];
        let line = prefix.iter().filter(|&&byte| byte == b'\n').count() + 1;
        let column = prefix.rsplit(|&byte| byte == b'\n').next().unwrap().len() + 1;
        Error::new(Span { line, column }, "source is not valid UTF-8")
    })
}
