mod common;

use common::{diagnostic, example, interpret, stdout};

const CURVE: &str = r"
rvunram W[3,3]{1,-1,0,1};
rvhyperelliptic C[W]{rvUnramElem[W]{0,1,2},rvUnramElem[W]{1,2},rvUnramElem[W]{1,1},
    rvUnramElem[W]{2,2,1},rvUnramElem[W]{2},rvUnramElem[W]{1}};
";

#[test]
fn the_ordered_closed_point_factor_agrees_with_independent_embedded_counts() {
    assert_eq!(
        stdout(example("closed-hasse-witt")),
        concat!(
            "RVHWWMatrix(W, 2){[[[3, 3, 7], [17, 8, 1]], [[15, 2, 4], [26, 24, 20]]]}\n",
            "rvHWWPoints[C, 1] = 26, trace=2\n",
            "rvHWWPoints[C, 2] = 794, trace=-64\n",
            "rvHWWWeil[C] = [1, -2, 34, -54, 729]\n",
            "rvHWWVerify[C] mod 27 = unit=[1, 25, 7], quotient=[1, 0, 0]\n",
            "rvHWWVerify[Lift] mod 27 = unit=[1, 25, 7], quotient=[1, 0, 0]\n",
        )
    );
    let wrong = format!(
        "{CURVE} {}",
        r"
        let M = rvHWWRatio[C];
        let M1 = rvHWWFrob[1](M);
        let M2 = rvHWWFrob[2](M);
        let reverse = rvHWWCompose(rvHWWCompose(M2,M1),M);
        rvHWWVerify[C](reverse);
    "
    );
    let error = diagnostic(interpret(&wrong));
    assert!(
        error.contains("does not divide the point-count Weil polynomial"),
        "{error}"
    );
    let nonscalar = diagnostic(interpret(&format!(
        "{CURVE} show rvHWWDescend(rvHWWUnit[C]);"
    )));
    assert!(
        nonscalar.contains("expected descent to the unramified base ring"),
        "{nonscalar}"
    );
}

#[test]
fn residue_degree_one_and_proper_subfield_cases_agree_with_prime_field_operators() {
    let source = r"
        rvhyperelliptic Base[5]{0,2,-3,1};
        rvunram U[5,3]{0,1};
        rvhyperelliptic A[U]{rvUnramElem[U]{0},rvUnramElem[U]{2},rvUnramElem[U]{-3},rvUnramElem[U]{1}};
        assert rvHWWDescend(rvHWW[A,3]) == rvHW[Base,3,3];
        assert rvHWWDescend(rvHWWRatio[A]) == rvHWUnit[Base,3];
        assert rvHWWDescend(rvHWWUnit[A]) == rvHWUnit[Base,3];
        rvHWWWeil[A];
        rvunram V[5,3]{2,0,1};
        rvhyperelliptic B[V]{rvUnramElem[V]{0},rvUnramElem[V]{2},rvUnramElem[V]{-3},rvUnramElem[V]{1}};
        let u = rvHWUnit[Base,3];
        assert rvHWWDescend(rvHWWRatio[B]) == u;
        assert rvHWWDescend(rvHWWUnit[B]) == rvHWCompose(u,u);
        rvHWWWeil[B];
    ";
    assert_eq!(
        stdout(interpret(source)),
        "rvHWWWeil[A] = [1, 2, 5]\nrvHWWWeil[B] = [1, 6, 25]\n"
    );
}

#[test]
fn the_genus_one_unit_agrees_with_the_signed_closed_legendre_norm() {
    let source = r"
        rvunram V[5,3]{2,0,1};
        let z = rvUnramElem[V]{0};
        let u = rvUnramElem[V]{1};
        let lambda = rvUnramElem[V]{0,1};
        rvhyperelliptic Legendre[V]{z,lambda,neg(add(lambda,u)),u};
        assert rvHWWEntry[1,1](rvHWWUnit[Legendre]) == rvLegendreUnit(lambda);
        assert rvHWWChar[1](rvHWWUnit[Legendre]) == neg(rvLegendreUnit(lambda));
        assert rvHWWEntry[1,1](rvHWWUnit[Legendre]) == rvUnramElem[V]{73};
        rvHWWVerify[Legendre];

        rvunram W[3,3]{1,-1,0,1};
        let z3 = rvUnramElem[W]{0};
        let u3 = rvUnramElem[W]{1};
        let lambda3 = rvUnramElem[W]{0,1};
        rvhyperelliptic Signed[W]{z3,lambda3,neg(add(lambda3,u3)),u3};
        assert rvHWWEntry[1,1](rvHWWUnit[Signed]) == rvLegendreUnit(lambda3);
        assert rvHWWEntry[1,1](rvHWWUnit[Signed]) == rvUnramElem[W]{4};
        rvHWWVerify[Signed];
    ";
    assert_eq!(
        stdout(interpret(source)),
        "rvHWWVerify[Legendre] mod 125 = unit=[1, 52], quotient=[1, 75]\nrvHWWVerify[Signed] mod 27 = unit=[1, 23], quotient=[1, 0]\n"
    );
}

#[test]
fn precision_reduction_is_coherent_in_the_same_witt_presentation() {
    let source = format!(
        "{CURVE} {}",
        r"
        rvunram V[3,2]{1,-1,0,1};
        rvhyperelliptic D[V]{rvUnramElem[V]{0,1,2},rvUnramElem[V]{1,2},rvUnramElem[V]{1,1},
            rvUnramElem[V]{2,2,1},rvUnramElem[V]{2},rvUnramElem[V]{1}};
        assert rvHWWReduce[V](rvHWW[C,2]) == rvHWW[D,2];
        assert rvHWWReduce[V](rvHWWRatio[C]) == rvHWWRatio[D];
        assert rvHWWReduce[V](rvHWWUnit[C]) == rvHWWUnit[D];
    "
    );
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn degree_even_models_count_the_quadratic_fiber_at_infinity() {
    let source = r"
        rvunram V[5,2]{0,1};
        let z=rvUnramElem[V]{0};
        let u=rvUnramElem[V]{1};
        rvhyperelliptic E[V]{z,u,z,u,rvUnramElem[V]{2}};
        rvHWWPoints[E,1];
        rvHWWWeil[E];
        rvHWWVerify[E];
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "rvHWWPoints[E, 1] = 4, trace=2\n",
            "rvHWWWeil[E] = [1, -2, 5]\n",
            "rvHWWVerify[E] mod 25 = unit=[1, 13], quotient=[1, 10]\n",
        )
    );
}

#[test]
fn nonordinary_curves_keep_raw_matrices_and_point_counts() {
    let setup = r"
        rvunram V[3,2]{0,1};
        let z=rvUnramElem[V]{0};
        let u=rvUnramElem[V]{1};
        rvhyperelliptic E[V]{z,neg(u),z,u};
    ";
    let source = format!("{setup} assert rvHWW[E,1] == zero[RVHWWMatrix(V,1)]; rvHWWWeil[E];");
    assert_eq!(stdout(interpret(&source)), "rvHWWWeil[E] = [1, 0, 3]\n");
    for query in [
        "show rvHWWRatio[E];",
        "show rvHWWUnit[E];",
        "rvHWWVerify[E];",
    ] {
        let error = diagnostic(interpret(&format!("{setup} {query}")));
        assert!(error.contains("nonordinary"), "{error}");
    }
}

#[test]
fn matrix_values_have_normal_type_and_finite_domain_checks() {
    let setup = r"
        rvunram V[3,1]{0,1};
        let z=rvUnramElem[V]{0};
        let u=rvUnramElem[V]{1};
        rvhyperelliptic E[V]{z,neg(u),z,u};
    ";
    let values = format!(
        "{setup} {}",
        r"
        let i=rvHWW[E,0];
        assert rvHWWInv(i) == i;
        assert add(i,neg(i)) == zero[RVHWWMatrix(V,1)];
        card(\x:RVHWWMatrix(V,1)=>S.hi);
    "
    );
    assert_eq!(stdout(interpret(&values)), "card = 3 / 3\n");
    for (source, expected) in [
        (
            "show rvHWWEntry[0,1](rvHWW[E,0]);",
            "matrix entry is out of range",
        ),
        (
            "show rvHWWChar[2](rvHWW[E,0]);",
            "characteristic coefficient is out of range",
        ),
        ("rvhyperelliptic Bad[V]{z,u,u,u};", "squarefree"),
        (
            "rvunram W[3,2]{0,1}; show rvHWWReduce[W](rvHWW[E,0]);",
            "no greater precision",
        ),
        (
            "rvunram W[3,1]{0,1}; rvhyperelliptic Bad[W]{z,neg(u),z,u};",
            "expected type RVUnram(W)",
        ),
        (
            "rvunram T[2,1]{0,1}; show zero[RVHWWMatrix(T,1)];",
            "odd residue characteristic",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup} {source}")));
        assert!(error.contains(expected), "{source}\n{error}");
    }
}
