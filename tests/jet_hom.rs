mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn named_arcs_reuse_jet_homology_and_literal_correction_lifts() {
    assert_eq!(
        stdout(example("jet-hom")),
        concat!(
            "JetHom(One, Cusp, Axis, 0){u^2*e(0,0, 4,0)}\n",
            "rvjet[One, Cusp, Axis] = dimension 6, blocks {1: 2, 2: 2}\n",
            "rvjet[One, Cusp, Curve] = dimension 18, blocks {2: 2, 7: 2}\n",
        )
    );
}

#[test]
fn raw_jet_composition_is_bilinear_and_obeys_leibniz() {
    let source = r"
        qmodule T[0] {a:1,b:0} {0(a)=>[b]};
        rvpath Line[[0],7,2,[[0,2]]];
        let f = jetterm[One,0,0,T,1,0,Line,0];
        let g = jetterm[T,0,0,One,0,0,Line,0];
        let h = jetterm[T,1,0,One,0,0,Line,0];
        let e = jetterm[One,0,0,One,0,0,Line,0];
        let u = jetterm[One,0,0,One,0,0,Line,1];
        assert ecomp(h,f) == e;
        assert ecomp(u,u) == zero[JetHom(One,One,Line,0)];
        assert ecomp(rvScale[2](h),rvScale[3](f)) == e;
        assert rvScale[0](f) == zero[JetHom(One,T,Line,1)];
        assert jetD(f) == rvScale[2](jetterm[One,0,0,T,0,0,Line,1]);
        assert jetD(g) == rvScale[2](jetterm[T,1,0,One,0,0,Line,1]);
        assert jetD(jetD(f)) == zero[JetHom(One,T,Line,1)];
        assert ecomp(g,jetD(f)) == rvScale[2](u);
        assert add(ecomp(jetD(g),f),ecomp(g,jetD(f))) == jetD(ecomp(g,f));
        let fId = \x : JetHom(One,T,Line,1) => x;
        assert fId(f) == f;
        card(\x : JetHom(One,One,Line,0) => S.hi);
    ";
    assert_eq!(stdout(interpret(source)), "card = 16 / 16\n");
}

#[test]
fn named_arcs_validate_parity_coefficients_endpoints_and_truncation() {
    let declaration = "rvpath A[[0],7,2,[[1]]]; rvpath B[[0],7,2,[[1]]];";
    for (source, reason) in [
        ("rvpath A[[0],7,2,[[1]]];", "already defined"),
        (
            "show zero[JetHom(One,One,A,2)];",
            "parity must be zero or one",
        ),
        (
            "show jetterm[One,0,0,One,0,0,A,2];",
            "outside the arc truncation",
        ),
        (
            "show zero[JetHom(One,One,Missing,0)];",
            "unknown Margolis arc",
        ),
        (
            "show rvScale[4](zero[JetHom(One,One,A,0)]);",
            "field coefficient 4 is out of range",
        ),
        (
            "show ecomp(zero[JetHom(One,One,A,0)],zero[JetHom(One,One,B,0)]);",
            "same Margolis arc",
        ),
        (
            "show ecomp(zero[JetHom(One,One,A,0)],zero[JetHom(One,S,A,0)]);",
            "composable JetHom endpoints",
        ),
        (
            "qmodule Empty[] {} {}; show zero[JetHom(Empty,One,A,0)];",
            "not declared for this module",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{declaration}\n{source}")));
        assert!(error.contains(reason), "{source}\n{error}");
    }
}
