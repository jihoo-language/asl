mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn miller_tate_returns_nontrivial_roots_in_genera_one_and_two() {
    assert_eq!(
        stdout(example("hyperelliptic-miller-tate")),
        "RVUnram(F)[4]\nRVUnram(W)[4]\nRVUnram(W)[3]\n"
    );
}

#[test]
fn miller_quadratic_pairing_agrees_at_closed_ramified_and_repeated_places() {
    assert_eq!(
        stdout(interpret(
            r"
        rvunram F[7,1]{0,1};
        let z=rvUnramElem[F]{};
        let o=rvUnramElem[F]{1};
        let t=rvUnramElem[F]{3};
        let f=rvUnramElem[F]{4};
        rvhyperelliptic C[F]{z,t,o,f,o,o};
        let closed=rvjac[C]{[1],[],[1]}{};
        let r=rvJacPoint[C](z,z);
        assert rvJacTate[2](closed,closed)==rvJacTate2(closed,closed);
        assert rvJacTate[2](r,closed)==rvJacTate2(r,closed);
        assert rvJacTate[2](r,r)==rvJacTate2(r,r);
        assert rvJacTate[2](r,zero[RVJac(C)])==o;
        assert rvJacTate[2](zero[RVJac(C)],r)==o;

        rvhyperelliptic D[F]{z,t,rvUnramElem[F]{-1},z,f,o};
        let a=rvJacPoint[D](z,z);
        let p=rvJacPoint[D](rvUnramElem[F]{5},o);
        assert rvJacTate[2](a,add(p,p))==rvJacTate2(a,add(p,p));
        assert rvJacTate[2](a,add(p,a))==rvJacTate2(a,add(p,a));
    "
        )),
        ""
    );
}

#[test]
fn miller_uses_nonmonic_infinity_units_and_weights_mixed_support() {
    assert_eq!(
        stdout(interpret(
            r"
        rvunram F[7,1]{0,1};
        let z=rvUnramElem[F]{};
        let o=rvUnramElem[F]{1};
        let f=rvUnramElem[F]{4};
        // With s=1+2x²+2x³, f_C=s²−4x⁶ and y−s has leading unit −2 at infinity.
        rvhyperelliptic C[F]{o,z,f,f,f,o};
        let p=rvJacPoint[C](z,o);
        let t=add(p,p);
        let q=rvJacPoint[C](rvUnramElem[F]{6},rvUnramElem[F]{2});
        assert rvJacTate[3](t,q)==rvUnramElem[F]{2};
        assert rvJacTate[3](t,neg(q))==f;
        assert rvJacTate[3](t,rvJacPoint[C](o,z))==o;

        rvunram W[29,1]{0,1};
        let zeroW=rvUnramElem[W]{};
        let oneW=rvUnramElem[W]{1};
        let fourW=rvUnramElem[W]{4};
        // The leading unit of y−2x−1 is 1/2 at a repeated point and 5 at the other point.
        rvhyperelliptic D[W]{oneW,fourW,fourW,zeroW,zeroW,zeroW,zeroW,oneW};
        let r=rvJacPoint[D](zeroW,oneW);
        let other=rvJacPoint[D](fourW,rvUnramElem[W]{14});
        let mixed=add(add(r,r),other);
        assert rvJacU[3](mixed)==oneW;
        assert rvJacTate[7](r,mixed)==rvUnramElem[W]{20};
    "
        )),
        ""
    );
}

#[test]
fn miller_uses_hasse_terms_in_small_characteristic_and_wide_field_orders() {
    assert_eq!(
        stdout(interpret(
            r"
        rvunram W[3,1]{1,0,1};
        let z=rvUnramElem[W]{};
        let o=rvUnramElem[W]{1};
        let d=rvUnramElem[W]{1,1};
        let dd=rvUnramElem[W]{0,2};
        // f_C=(x⁴+2x+d)²−x⁸. The shared local unit after normalization is 1/(2d).
        rvhyperelliptic C[W]{dd,d,o,z,neg(d),o};
        let t=rvJacPoint[C](z,d);
        let expected=rvUnramElem[W]{1,2};
        assert rvJacTate[8](t,t)==expected;
        assert rvJacTate[8](t,add(t,t))==rvUnramMul(expected,expected);
        assert rvJacTate[8](t,rvJacPoint[C](d,z))==neg(o);

        rvunram L[65537,1]{-3,0,0,0,1};
        let a=rvUnramElem[L]{0,1};
        let zeroL=rvUnramElem[L]{};
        let oneL=rvUnramElem[L]{1};
        rvhyperelliptic G[L]{zeroL,rvUnramElem[L]{-3},zeroL,zeroL,rvUnramElem[L]{21846},oneL};
        let r=rvJacPoint[G](zeroL,zeroL);
        let q=rvJacPoint[G](a,oneL);
        assert rvJacTate[2](r,q)==rvUnramElem[L]{-1};
    "
        )),
        ""
    );
}

#[test]
fn miller_requires_valid_order_torsion_and_named_curve() {
    let source = r"
        rvunram F[7,1]{0,1};
        let z=rvUnramElem[F]{};
        let o=rvUnramElem[F]{1};
        rvhyperelliptic C[F]{o,z,z,o};
        rvhyperelliptic D[F]{o,z,z,o};
        let t=rvJacPoint[C](z,o);
        let r=zero[RVJac(C)];
    ";
    for (code, expected) in [
        ("show rvJacTate[0](r,r);", "at least 2"),
        ("show rvJacTate[1](r,r);", "at least 2"),
        ("show rvJacTate[4](r,r);", "must divide q−1"),
        ("show rvJacTate[2](t,r);", "must have order dividing 2"),
        ("show rvJacTate[2](t,t);", "must have order dividing 2"),
        ("show rvJacTate[18446744073709551616](r,r);", "exceeds u64"),
        (
            "show rvJacTate[2](r,zero[RVJac(D)]);",
            "expected type RVJac(C)",
        ),
    ] {
        let problem = diagnostic(interpret(&format!("{source}{code}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
