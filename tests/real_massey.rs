mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn full_affine_coset_and_sparse_tower_have_actual_witnesses() {
    assert_eq!(
        stdout(example("real-motivic-massey")),
        concat!(
            "rmTripleDim = 1\n",
            "rmTripleContains = yes\n",
            "rmTripleContains = yes\n",
            "rmTripleContains = yes\n",
            "rmTripleContains = no\n",
            "rmTripleDim = 0\n",
            "rhoorder = 4\n",
            "rhoorder = 32\n",
        ),
    );
}

#[test]
fn coset_does_not_change_when_a_raw_representative_changes_by_a_boundary() {
    let source = r"
        let h0 = rmcword[0,0]{([],[1])};
        let h1 = rmcword[0,0]{([1],[])};
        let boundary = cobarD(rmcword[0,1]{});
        let from_boundary = rmTriple(boundary,h0,h1);
        let from_class = rmTriple(zero[RMExt(1,0,-1)],h0,h1);
        assert from_boundary == from_class;
        rmTripleContains(from_boundary,zero[RMExt(2,3,0)]);
        let A = rmTriple(rmcword[1,0]{},zero[RMCobar(0,0,0)],h0);
        assert rmTripleBase(A) == zero[RMExt(0,0,-1)];
        rmTripleDim(A);
    ";
    assert_eq!(
        stdout(interpret(source)),
        "rmTripleContains = yes\nrmTripleDim = 0\n"
    );
}

#[test]
fn masseys_diagnose_noncycles_nonzero_products_and_invalid_directions() {
    for (source, reason) in [
        (
            "let h=rmcword[0,0]{([],[1])}; show rmTriple(h,h,h);",
            "A cup B is not a boundary",
        ),
        (
            "let h=rmcword[0,0]{([],[1])}; show rmTriple(rmcword[1,0]{},h,h);",
            "B cup C is not a boundary",
        ),
        (
            "let u=rmcword[0,1]{}; show rmMassey(u,zero[RMCobar(1,1,0)],u);",
            "require closed inputs",
        ),
        (
            "let x=rmcword[0,0]{([],[1])}; show rmTripleDir[0](rmTriple(x,rmcword[1,0]{},x));",
            "outside the indeterminacy basis",
        ),
        (
            "show rmMasseyNull(rmcword[0,0]{},rmcword[0,0]{});",
            "requires positive total cobar length",
        ),
        (
            "let x=rmcword[0,0]{([],[1])}; rmTripleContains(rmTriple(x,rmcword[1,0]{},x),x);",
            "an RMCobar cochain or RMExt class in degree",
        ),
        (
            r"card(\x : RMTriple(1,1,-1) => S.hi);",
            "RMTriple types do not enumerate ambient cosets",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(reason), "{source}\n{error}");
    }
}
