mod common;

use common::{diagnostic, example, interpret, stdout};
use std::fmt::Write;

#[test]
fn height_inclusion_delays_the_bockstein_and_doubles_the_integral_class() {
    assert_eq!(
        stdout(example("bp-height")),
        "BPBock(2, 2, 4, 1, 8, 4){g0}\nBPExt(2, 2, 1, 8){2*g0}\n"
    );
}

#[test]
fn raw_raising_reindexes_every_slot_and_commutes_with_cofaces() {
    let mut source = String::from(
        r"
        let v = bpv[3,2,2,2];
        let t = bpt[3,2,2,1,2];
        let u = bpt[3,2,2,2,1];
        let z = add(bpMul(v,t),bpPow[2](u));
        let widened = add(bpMul(bpv[3,4,2,2],bpt[3,4,2,1,2]),bpPow[2](bpt[3,4,2,2,1]));
        assert bpRaise[4](z) == widened;
        assert bpRaise[4](bpRaise[3](z)) == widened;
        assert bpRaise[2](z) == z;
        assert bpRaise[4](bpD(z)) == bpD(widened);
        assert bpRaise[4](bpNormalize(z)) == bpNormalize(widened);
        assert bpRaise[4](bpDegeneracy[0](z)) == bpDegeneracy[0](widened);
        assert bpRaise[4](bpReduce[2](z)) == bpReduce[2](widened);
        assert bpRaise[4](bpCup(bpt[3,2,1,1,2], z)) == bpCup(bpt[3,4,1,1,2],widened);
        assert bpRaise[2](bpint[3,0,1,-2]) == bpint[3,2,1,-2];
    ",
    );
    for i in 0..=3 {
        writeln!(
            source,
            "assert bpRaise[4](bpCoface[{i}](z)) == bpCoface[{i}](widened);"
        )
        .unwrap();
        writeln!(source, "assert bpRaise[4](bpCoface[{i}](bpReduce[2](z))) == bpCoface[{i}](bpReduce[2](widened));").unwrap();
    }
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn page_raising_commutes_with_projection_differential_and_product() {
    let source = r"
        let t = bpt[3,1,1,1,1];
        let v = bpv[3,1,0,1];
        let x = bpNovClass[0,4,0](t);
        let y = bpNovClass[0,4,1](v);
        assert bpRaise[2](x) == bpNovClass[0,4,0](bpRaise[2](t));
        assert bpRaise[2](y) == bpNovClass[0,4,1](bpRaise[2](v));
        let z = bpCup(x,y);
        assert bpRaise[2](z) == bpCup(bpRaise[2](x),bpRaise[2](y));
        assert bpRaise[2](bpNovD(z)) == bpNovD(bpRaise[2](z));
        assert bpRaise[3](bpRaise[2](z)) == bpRaise[3](z);
        let b = bpBockgen[3,1,1,0,4,1000,0];
        assert bpRaise[2](b) == bpBockgen[3,2,1,0,4,1000,0];
        assert bpBockD(bpRaise[2](b)) == bpRaise[2](bpBockD(b));
        assert bpRaise[2](add(b,b)) == pscale[2](bpRaise[2](b));
    ";
    assert_eq!(stdout(interpret(source)), "");
}

#[test]
fn height_is_checked_on_cochain_and_class_types() {
    for (term, explanation) in [
        ("bpRaise[1](zero[BP(2,2,0)])", "cannot decrease"),
        ("bpRaise[1](zero[BPI(2,2,0,1)])", "cannot decrease"),
        ("bpRaise[1](zero[BPBock(2,2,1,0,0,0)])", "cannot decrease"),
        ("bpRaise[1](zero[BPNov(2,2,1,0,0,0)])", "cannot decrease"),
        ("bpRaise[1](zero[BPExt(2,2,1,8)])", "cannot decrease"),
        ("bpRaise[9](zero[BP(2,2,0)])", "height exceeds"),
        (
            "bpRaise[2](zero[Ext(0,0)])",
            "a BP cochain, page class, or BPExt class",
        ),
    ] {
        let error = diagnostic(interpret(&format!("show {term};")));
        assert!(error.contains(explanation), "{term}\n{error}");
    }
}
