mod adams;
mod binary_field;
mod binary_linear;
mod bp;
mod bp_classes;
mod bp_ext;
mod bp_height;
mod bp_higher;
mod bp_lattice;
mod bp_mod_ext;
mod bp_page;
mod bp_polynomial;
mod bp_rel;
mod bp_smith;
mod bp_triple;
mod bp_vpage;
mod bracket;
mod ce;
mod cobar;
mod cobar_squares;
mod cohomology;
mod comodule;
mod complex;
mod cyclic;
mod dependent_sets;
mod engine;
mod eo;
mod error;
mod integer_smith;
mod jet_hom;
mod jet_page;
mod khovanov;
mod lee;
mod lee_linear;
mod legendre;
mod legendre_closed;
mod legendre_cm;
mod linear;
mod map_bracket;
mod map_cochain;
mod map_defining;
mod map_ext;
mod map_may;
mod map_system;
mod map_transfer;
mod may;
mod may_lift;
mod milnor;
mod mixed;
mod mixed_bars;
mod mixed_cone;
mod mixed_hom;
mod mixed_homology;
mod mixed_minimal;
mod mixed_page;
mod mixed_tensor;
mod model;
mod motivic;
mod nakamura;
mod observations;
mod odd_cobar;
mod odd_ext;
mod odd_massey;
mod odd_may;
mod odd_profile;
mod odd_restriction_transfer;
mod odd_transfer;
mod page_transition;
mod prime_linear;
mod profile;
mod q_koszul;
mod q_module;
mod rank_jets;
mod rank_variety;
mod real_bars;
mod real_ext;
mod real_extensions;
mod real_massey;
mod real_motivic;
mod real_rho;
mod real_transfer;
mod restriction_transfer;
mod rho_bock;
mod rho_frobenius;
mod rho_localized;
mod rho_minimal;
mod rv_arithmetic;
mod rv_as_cartier;
mod rv_as_characters;
mod rv_as_cover;
mod rv_as_derham;
mod rv_as_odd;
mod rv_as_odd_cartier;
mod rv_as_odd_derham;
mod rv_as_odd_tate;
mod rv_as_two;
mod rv_as_two_cartier;
mod rv_bar;
mod rv_bar_dual;
mod rv_bar_page;
mod rv_barcobar;
mod rv_cartier_arithmetic;
mod rv_conformal;
mod rv_crystal;
mod rv_curve;
mod rv_cyclic_homology;
mod rv_cyclotomic;
mod rv_descent;
mod rv_field;
mod rv_flag;
mod rv_fusion;
mod rv_gauss;
mod rv_hh;
mod rv_hh_homology;
mod rv_hh_mixed;
mod rv_hive;
mod rv_hochster;
mod rv_hodge;
mod rv_hom;
mod rv_ideal;
mod rv_iwahori;
mod rv_iwhecke;
mod rv_iwkl;
mod rv_koszul;
mod rv_laurent;
mod rv_lens;
mod rv_negative;
mod rv_newton;
mod rv_ore;
mod rv_pcurv;
mod rv_period;
mod rv_phodge;
mod rv_poly;
mod rv_qconnection;
mod rv_qgr;
mod rv_qp_quad;
mod rv_qregular;
mod rv_rees;
mod rv_res_space;
mod rv_satake;
mod rv_scheme;
mod rv_shintani;
mod rv_sphhecke;
mod rv_t_bracket;
mod rv_t_system;
mod rv_t_transfer;
mod rv_tl;
mod rv_tree;
mod rv_unram;
mod rv_verlinde;
mod rv_witt_cartier;
mod rv_witt_cover;
mod rv_witt_tower;
mod rv_wz;
mod rv_wzw;
mod set_family;
mod set_j;
mod set_kan;
mod set_logic;
mod set_monodromy;
mod set_topology;
mod sheaf;
mod sheaf_base_change;
mod sheaf_bc_comma;
mod sheaf_cellular;
mod sheaf_deform;
mod sheaf_ext;
mod sheaf_higher;
mod sheaf_hyper;
mod sheaf_image;
mod sheaf_kuranishi;
mod sheaf_leray;
mod sheaf_page;
mod sheaf_serre;
mod sheaf_stalk;
mod sheaf_support;
mod sheaf_truncate;
mod stable_marks;
mod steenrod;
mod syntax;
mod transfer;
mod tropical_matroid;
mod v_bars;
mod v_koszul;
mod vbock;
mod wide_integer;

use std::env;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Read, Write};
use std::panic;
use std::process::ExitCode;
use std::thread;

use crate::engine::Engine;
use crate::error::{Error, Result, Span};
use crate::syntax::Parser;

const MAX_SOURCE_BYTES: usize = 512 * 1024;
// In debug builds the expression checker can use about 600 KiB per nested
// level. Reserve enough stack for the 64-level syntax limit on platforms
// whose default main-thread stack is smaller.
const INTERPRETER_STACK_BYTES: usize = 64 * 1024 * 1024;

const HELP: &str = "\
ASL, the Alexandrov Section Lambda interpreter

Usage: asl <PATH>
       asl -
       asl --help

Pass a UTF-8 ASL source file or use - to read standard input.
The show, card, betti, shdim, shhyperdim, shExtDim, shedim, shserreshift, bpext, bprelext, bpHeightPage, bpHeightImage, bporder, extdim, mextdim, rmextdim, rmrhoextdim, rhoorder, rhodepth, rhohidden, rmTripleDim, rmTripleContains, rmbars, rmbarRank, rhobockdim, pextdim, extandim, pextandim, xextdim, mapextdim, qextdim, rvdim, rvInvariantFactors, rvNewton, rvHodge, rvPosition, rvMuAt, rvIw, rvIwLength, rvIwPosition, rvIwAt, rvIwKL, rvIwSchubertKL, rvHeckeCoeff, rvHeckeAt, rvSphCoeff, rvSphAt, rvHiveCount, rvHiveAt, rvHorn, rvQGrCoeff, rvQGrGW, rvCTExact, rvCTWZCheck, rvOreCoeff, rvQGrPeriod, rvQGrRegularized, rvQGrRegularizedRhs, rvQGrCurvatureChar, rvQGrCurvatureCoeff, rvQGrCurvatureHorizontal, rvFusionCoeff, rvCycloCoeff, rvConformalDim, rvConformalSeries, rvConformalNum, rvConformalDen, rvConformalMod, rvConformalExact, rvSatakeCoeff, rvSatakeAt, rvTreeShell, rvTreeBall, rvTreeCount, rvTreeClosureCount, rvspectrum, rvorbit, rvclosed, rvzeta, rvjet, rvminors, rvorders, rvFrobOrders, rvOrdAt, rvASGenus, rvASPrank, rvASPoints, rvASLPolynomial, rvASNewton, rvASCartierRank, rvASaNumber, rvWittBreaks, rvWittGenus, rvWittPoints, rvWittLPolynomial, rvWittNewton, rvWittCharacterL, rvWittGaps, rvWittBasis, rvWittCartierRank, rvWittaNumber, rvWittTowerBreaks, rvWittTowerGenus, rvWittTowerPoints, rvWittTowerCharacterL, rvWittTowerLPolynomial, rvWittTowerNewton, rvbetti, rvmbetti, rvhochster, rvbardim, rvbaredim, rvpathextdim, rvHHdim, rvNCdim, rvHNdim, rvHPdim, rvHCdim, mixdim, mixmindim, mixbars, mixbardim, mixbarrank, mixbarperiodicdim, rvgroebner, rvhilbert, vextdim, vbars, vbardim, vbockdim, maydim, pmaydim, mayandim, pmayandim, mapmaydim, ceandim, pceandim, and homotopy statements print to standard output.
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
    thread::scope(|scope| {
        let worker = thread::Builder::new()
            .name("asl-interpreter".into())
            .stack_size(INTERPRETER_STACK_BYTES)
            .spawn_scoped(scope, || interpret(path))
            .map_err(|error| Error::plain(format!("cannot start interpreter: {error}")))?;
        match worker.join() {
            Ok(result) => result,
            Err(payload) => panic::resume_unwind(payload),
        }
    })
}

fn interpret(path: &OsString) -> Result<()> {
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
