mod common;

use common::{diagnostic, example, interpret, stdout};

const GENUS_ONE: &str = "rvunram F[17,1]{0,1};
    rvhyperelliptic E[F]{rvUnramElem[F]{1},rvUnramElem[F]{0},rvUnramElem[F]{0},rvUnramElem[F]{1}};
    rvunram W[17,1]{14,0,1}; rvjacbase L[E,W];
    rvunram R[3,2]{0,1};
    let p=rvJacPoint[L](rvUnramElem[W]{1},rvUnramElem[W]{6});
    let q=rvJacPoint[L](rvUnramElem[W]{5},rvUnramElem[W]{0,5});
    let r=rvJacPoint[L](rvUnramElem[W]{16},rvUnramElem[W]{0});
    let root=rvUnramElem[W]{15,1};";

#[test]
fn ordinary_genus_two_recovers_the_full_three_primary_group_at_height_two() {
    let output = stdout(example("hyperelliptic-kummer-brauer"));
    assert!(output.contains(
        "rvJacPowerFrob[B, S, 1] = level=9, multiplier=1, matrix=[[7, 0, 3, 0], [0, 7, 0, 3], [3, 0, 4, 0], [0, 3, 0, 4]], inverse=[[4, 0, 6, 0], [0, 4, 0, 6], [6, 0, 7, 0], [0, 6, 0, 7]], characteristic=[1, 5, 6, 5, 1]"
    ), "{output}");
    for expected in [
        "rvKummerK3Brauer[B3, S, 1] = level=3, prime=3, height=1, exceptional_rank=10, abelian_free_rank=4, smith_exponents=[1, 1, 1, 1, 1, 1], brauer_invariant_factors=[3, 3], complete_primary=false",
        "rvKummerK3Brauer[B, S, 1] = level=9, prime=3, height=2, exceptional_rank=10, abelian_free_rank=4, smith_exponents=[1, 1, 2, 2, 2, 2], brauer_invariant_factors=[3, 3], complete_primary=true",
        "rvKummerK3Brauer[B, S, 3] = level=9, prime=3, height=2, exceptional_rank=10, abelian_free_rank=4, smith_exponents=[2, 2, 2, 2, 2, 2], brauer_invariant_factors=[9, 9], complete_primary=false",
        "rvKummerK3BrauerArtinTate[B, S, 3] = level=9, prime=3, height=2, exceptional_rank=10, abelian_free_rank=4, smith_exponents=[2, 2, 2, 2, 2, 2], brauer_invariant_factors=[9, 9], complete_primary=false, artin_tate_product_magnitude=10036224, artin_tate_valuation=4, detected_order_valuation=4, complete_with_artin_tate=true",
        "rvKummerK3BrauerArtinTate[B, S, 9] = level=9, prime=3, height=2, exceptional_rank=10, abelian_free_rank=4, smith_exponents=[2, 2, 2, 2, 2, 2], brauer_invariant_factors=[9, 9], complete_primary=false, artin_tate_product_magnitude=236759030424576, artin_tate_valuation=6, detected_order_valuation=4, complete_with_artin_tate=false",
        "9216",
        "rvKummerK3Brauer[Z, X, 2] = level=3, prime=3, height=1, exceptional_rank=16, abelian_free_rank=6, smith_exponents=[1, 1, 1, 1, 1, 1], brauer_invariant_factors=[], complete_primary=true",
    ] {
        assert!(output.contains(expected), "{expected}\n{output}");
    }
}

#[test]
fn another_genuine_level_nine_marking_gives_the_same_brauer_group() {
    let script = include_str!("../examples/hyperelliptic-kummer-brauer.asl");
    let output = stdout(interpret(&format!(
        "{script}
         rvjacpowerbasis D[SL,3,2](rvUnramElem[L]{{33}}){{e1,e2,add(e1,f1),f2}};
         rvJacPowerFrob[D,S,1];
         rvKummerK3Brauer[D,S,1];"
    )));
    assert!(
        output.contains("matrix=[[4, 0, 3, 0], [0, 7, 0, 3], [3, 0, 7, 0], [0, 3, 0, 4]]"),
        "{output}"
    );
    let for_basis = |basis| {
        output
            .lines()
            .find(|line| line.starts_with(&format!("rvKummerK3Brauer[{basis}, S, 1]")))
            .unwrap()
            .split_once(" = ")
            .unwrap()
            .1
    };
    assert_eq!(for_basis("B"), for_basis("D"));
}

#[test]
fn prime_power_coordinates_and_similitudes_over_zmod_nine() {
    let output = stdout(interpret(&format!(
        "{GENUS_ONE}
         rvjacpowerbasis B[L,3,2](root){{p,q}};
         rvjacpowerbasis C[L,3,2](root){{p,add(p,q)}};
         rvjacpowerbasis D[L,3,2](rvUnramElem[W]{{7,13}}){{p,rvJacScale[2](q)}};
         assert rvJacWeil[9](p,q)==root;
         let t=add(rvJacScale[7](p),rvJacScale[8](q));
         assert rvJacPowerCoord[B,0,R](t)==rvUnramElem[R]{{7}};
         assert rvJacPowerCoord[B,1,R](t)==rvUnramElem[R]{{8}};
         assert rvJacPowerCoord[C,0,R](t)==rvUnramElem[R]{{8}};
         assert rvJacPowerCoord[C,1,R](t)==rvUnramElem[R]{{8}};
         rvJacPowerTransition[B,C]; rvJacPowerTransition[B,D];
         rvJacPowerFrob[B,E,1]; rvJacPowerFrob[C,E,1];
         rvJacPowerFrob[C,E,2]; rvJacPowerFrob[C,E,1000001];"
    )));
    assert_eq!(
        output,
        concat!(
            "rvJacPowerTransition[B, C] = root_exponent=1, matrix=[[1, 8], [0, 1]], inverse=[[1, 1], [0, 1]]\n",
            "rvJacPowerTransition[B, D] = root_exponent=5, matrix=[[1, 0], [0, 5]], inverse=[[1, 0], [0, 2]]\n",
            "rvJacPowerFrob[B, E, 1] = level=9, multiplier=8, matrix=[[1, 0], [0, 8]], inverse=[[1, 0], [0, 8]], characteristic=[1, 0, 8]\n",
            "rvJacPowerFrob[C, E, 1] = level=9, multiplier=8, matrix=[[1, 2], [0, 8]], inverse=[[1, 2], [0, 8]], characteristic=[1, 0, 8]\n",
            "rvJacPowerFrob[C, E, 2] = level=9, multiplier=1, matrix=[[1, 0], [0, 1]], inverse=[[1, 0], [0, 1]], characteristic=[1, 7, 1]\n",
            "rvJacPowerFrob[C, E, 1000001] = level=9, multiplier=8, matrix=[[1, 2], [0, 8]], inverse=[[1, 2], [0, 8]], characteristic=[1, 0, 8]\n",
        )
    );
}

#[test]
fn rejects_nonprimitive_markings_wrong_curves_and_incompatible_levels() {
    for (query, expected) in [
        ("rvjacpowerbasis A[L,3,0](root){p,q};", "height must be positive"),
        ("rvjacpowerbasis A[L,3,100](root){p,q};", "modulus exceeds the supported integer range"),
        ("rvjacpowerbasis A[L,17,2](root){p,q};", "differ from the residue characteristic"),
        ("rvjacpowerbasis A[L,3,2](rvUnramElem[W]{8,15}){p,q};", "primitive root of unity of order 9"),
        ("rvjacpowerbasis A[L,3,2](root){p,r};", "must be killed by 9"),
        ("rvjacpowerbasis A[L,3,2](root){q,p};", "incorrect symplectic pairing"),
        ("show rvJacPowerCoord[B,2,R](p);", "index is out of range"),
        ("show rvJacPowerCoord[B,0,R](r);", "requires a class killed by 9"),
        ("rvunram R2[3,1]{0,1}; show rvJacPowerCoord[B,0,R2](p);", "degree-one residue ring"),
        ("rvJacTorsionFrob[B,E,1];", "unknown Jacobian torsion basis"),
        ("rvhyperelliptic T[F]{rvUnramElem[F]{1},rvUnramElem[F]{1},rvUnramElem[F]{0},rvUnramElem[F]{1}}; rvJacPowerFrob[B,T,1];", "not the declared base extension"),
    ] {
        let error = diagnostic(interpret(&format!(
            "{GENUS_ONE} rvjacpowerbasis B[L,3,2](root){{p,q}}; {query}"
        )));
        assert!(error.contains(expected), "{error}");
    }
}

#[test]
fn excludes_the_two_primary_gluing_and_degree_zero() {
    let script = include_str!("../examples/hyperelliptic-jacobian-level.asl");
    let error = diagnostic(interpret(&format!("{script}\nrvKummerK3Brauer[B2,C,1];")));
    assert!(error.contains("requires an odd torsion prime"), "{error}");
    let script = include_str!("../examples/hyperelliptic-kummer-brauer.asl");
    let error = diagnostic(interpret(&format!("{script}\nrvKummerK3Brauer[B,S,0];")));
    assert!(
        error.contains("extension degree must be positive"),
        "{error}"
    );
}
