mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn marking() -> &'static str {
    include_str!("../examples/hyperelliptic-kummer-mukai.asl")
        .split_once("// Mukai calculations below.")
        .unwrap()
        .0
}

#[test]
fn tensor_and_ideal_pairings_braid_and_mixed_intersections() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-mukai")),
        concat!(
            "rvMukaiPair = -4\n",
            "rvMukaiPair = 4\n",
            "rvMukaiPair = -1\n",
            "rvMukaiPair = -1\n",
            "rvMukaiHilbFujiki[2, 2] = 12\n",
            "rvMukaiHilbFujiki[3, 3] = -384\n",
        )
    );
}

#[test]
fn the_embedding_round_trips_all_tropes_under_two_actual_markings() {
    let mut script = marking().to_owned();
    script.push_str("rvjactorsionbasis Twin[L,2](tw){t1,e2,t0,f2};\n");
    for basis in ["B", "Twin"] {
        for a in 0..16 {
            writeln!(script, "assert rvK3ToKummerDiv(rvK3FromKummerDiv(rvKummerTropeDiv[{basis},{a}]))==rvKummerTropeDiv[{basis},{a}];").unwrap();
        }
    }
    script.push_str(
        r"
        let h=rvKummerHyperplane[B];
        let e=rvKummerNodeDiv[B,0];
        let t=rvKummerTropeDiv[B,0];
        let gh=rvK3FromKummerDiv(h);
        assert gh==rvK3Lattice[B]{[0,2,0,0,2,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]};
        assert rvK3FromKummerDiv(add(e,t))==add(rvK3FromKummerDiv(e),rvK3FromKummerDiv(t));
        rvK3LatticePair(gh,gh);
        rvK3LatticePair(rvK3FromKummerDiv(e),rvK3FromKummerDiv(t));
    ",
    );
    assert_eq!(
        stdout(interpret(&script)),
        "rvK3LatticePair = 4\nrvK3LatticePair = 1\n"
    );
}

#[test]
fn reflections_pairings_and_fujiki_allow_full_k3_directions() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{}
            let full=rvK3Lattice[B]{{[2,0,0,0,0,2],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}};
            let y=rvMukai[0,0](full);
            assert rvMukaiC(y)==full;
            assert add(y,neg(y))==zero[RVMukai(B)];
            rvMukaiHilbFujiki[2,0](y,y);
            let e0=rvK3FromKummerDiv(rvKummerNodeDiv[B,0]);
            let e1=rvK3FromKummerDiv(rvKummerNodeDiv[B,1]);
            let r0=rvMukai[0,0](e0);
            let r1=rvMukai[0,0](e1);
            let line=rvMukaiLine(rvKummerTropeDiv[B,0]);
            assert rvMukaiTwist(r0,rvMukaiTwist(r1,line))==rvMukaiTwist(r1,rvMukaiTwist(r0,line));
            assert rvMukaiC(rvMukaiTwist(r0,line))==rvK3Reflect(e0,rvMukaiC(line));
            assert rvMukaiTwist(r0,r0)==neg(r0);
            rvMukaiPair(line,r0);
            rvMukaiPair(rvMukaiTwist(r1,line),rvMukaiTwist(r1,r0));
            ",
            marking()
        ))),
        "rvMukaiHilbFujiki[2, 0] = 48\nrvMukaiPair = 1\nrvMukaiPair = 1\n"
    );
}

#[test]
fn wide_cancellations_and_signed_boundary_values_are_exact() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{}
            let h=rvKummerHyperplane[B];
            let big=rvKummerDivScale[42535295865117307932921825928971026432](h);
            let center=rvK3FromKummerDiv(big);
            let isotropic=rvMukai[42535295865117307932921825928971026432,85070591730234615865843651857942052864](center);
            rvMukaiPair(isotropic,isotropic);
            let before=rvMukai[85070591730234615865843651857942052864,0](neg(center));
            let after=rvMukai[85070591730234615865843651857942052864,0](center);
            assert rvMukaiTensor(before,h)==after;
            assert rvMukaiTensor(after,neg(h))==before;
            let minimum=rvMukaiIdeal[B,170141183460469231731687303715884105729];
            assert minimum==rvMukai[1,-170141183460469231731687303715884105728](zero[RVK3Lattice(B)]);
            ",
            marking()
        ))),
        "rvMukaiPair = 0\n"
    );
}

#[test]
fn rejects_wrong_markings_off_subspace_inputs_nonroots_and_bad_degrees() {
    let setup = format!(
        "{} let h=rvKummerHyperplane[B]; let v=rvMukaiLine(h); let i2=rvMukaiIdeal[B,2];",
        marking()
    );
    for (code, expected) in [
        ("show rvK3ToKummerDiv(rvK3LatticeBasis[B,0]);", "outside the marked Kummer divisor sublattice"),
        ("show rvK3ToKummerDiv(v);", "requires an RVK3Lattice vector"),
        ("show rvMukai[1,0](h);", "requires an RVK3Lattice vector"),
        ("show rvMukaiTwist(i2,v);", "requires a vector of square"),
        ("rvMukaiHilbFujiki[2,1](i2,i2);", "orthogonal to the ideal sheaf vector"),
        ("rvMukaiHilbFujiki[1,1](i2,i2);", "requires n at least two"),
        ("rvMukaiHilbFujiki[2,5](i2,i2);", "between zero and 2n"),
        ("show rvMukaiIdeal[B,170141183460469231731687303715884105730];", "Mukai coordinate exceeds the signed i128 range"),
        ("rvjactorsionbasis Twin[L,2](tw){t0,e2,t1,f2}; rvMukaiPair(v,rvMukaiIdeal[Twin,0]);", "expected type RVMukai(B)"),
        ("rvjactorsionbasis Twin[L,2](tw){t0,e2,t1,f2}; show rvMukaiTensor(v,rvKummerHyperplane[Twin]);", "expected type RVKummerDiv(B)"),
    ] {
        let problem = diagnostic(interpret(&format!("{setup}{code}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
