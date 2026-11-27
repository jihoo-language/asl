mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn miller_weil_matches_explicit_functions_in_genera_one_and_two() {
    assert_eq!(
        stdout(example("hyperelliptic-miller-weil")),
        "RVUnram(F)[4]\nRVUnram(W)[5]\n"
    );
}

#[test]
fn miller_weil_agrees_with_the_quadratic_pairing_at_shared_branch_places() {
    assert_eq!(
        stdout(interpret(
            r"
        rvunram F[7,1]{0,1};
        let z=rvUnramElem[F]{};
        let o=rvUnramElem[F]{1};
        let t=rvUnramElem[F]{3};
        let f=rvUnramElem[F]{4};
        rvhyperelliptic C[F]{z,t,rvUnramElem[F]{-1},z,f,o};
        let a=rvJacPoint[C](z,z);
        let b=rvJacPoint[C](o,z);
        let c=rvJacPoint[C](rvUnramElem[F]{2},z);
        let ab=add(a,b);
        let bc=add(b,c);
        assert rvJacWeil[2](a,b)==rvJacWeil2(a,b);
        assert rvJacWeil[2](ab,bc)==rvJacWeil2(ab,bc);
        assert rvJacWeil[2](ab,ab)==o;
        assert rvJacWeil[2](bc,zero[RVJac(C)])==o;
        assert rvJacWeil[4](ab,bc)==o;

        rvhyperelliptic D[F]{z,t,o,f,o,o};
        let branch=rvjac[D]{[1],[],[1]}{};
        let r=rvJacPoint[D](z,z);
        assert rvJacWeil[2](branch,branch)==rvJacWeil2(branch,branch);
        assert rvJacWeil[2](branch,r)==rvJacWeil2(branch,r);
    "
        )),
        ""
    );
}

#[test]
fn miller_weil_accepts_orders_not_dividing_the_field_size_minus_one() {
    assert_eq!(
        stdout(interpret(
            r"
        rvunram F[5,1]{0,1};
        let z=rvUnramElem[F]{};
        let o=rvUnramElem[F]{1};
        rvhyperelliptic C[F]{o,z,o,rvUnramElem[F]{2},rvUnramElem[F]{4},o};
        let p=rvJacPoint[C](z,o);
        let t=add(p,p);
        assert rvJacScale[3](t)==zero[RVJac(C)];
        assert rvJacWeil[3](t,neg(t))==o;

        // Over F25 with a²=2, the functions y−(1+3x²+x³) and
        // y−a(1+3x+3x²+2x³) have divisors 3T and 3U. Their normalized
        // evaluations are 3 and 1+a, respectively.
        rvunram W[5,1]{-2,0,1};
        rvjacbase L[C,W];
        let first=rvJacExtend[L](t);
        let q=rvJacPoint[L](rvUnramElem[W]{4,2},rvUnramElem[W]{1});
        let r=rvJacPoint[L](rvUnramElem[W]{4,3},rvUnramElem[W]{-1});
        let second=add(q,r);
        let value=rvJacWeil[3](first,second);
        assert value==rvUnramElem[W]{2,3};
        assert rvJacFrob[C,1](second)==neg(second);
        assert rvJacWeil[3](rvJacFrob[C,1](first),rvJacFrob[C,1](second))==rvUnramFrob[1](value);
    "
        )),
        ""
    );
}

#[test]
fn miller_weil_checks_both_torsion_arguments_and_the_order() {
    let source = r"
        rvunram F[7,1]{0,1};
        let z=rvUnramElem[F]{};
        let o=rvUnramElem[F]{1};
        let two=rvUnramElem[F]{2};
        rvhyperelliptic C[F]{two,z,z,o};
        rvhyperelliptic D[F]{two,z,z,o};
        let r=zero[RVJac(C)];
        let t=rvJacPoint[C](z,rvUnramElem[F]{3});
    ";
    for (code, expected) in [
        ("show rvJacWeil[1](r,r);", "at least 2"),
        (
            "show rvJacWeil[14](r,r);",
            "coprime to the residue characteristic",
        ),
        (
            "show rvJacWeil[2](t,r);",
            "first argument must have order dividing 2",
        ),
        (
            "show rvJacWeil[2](r,t);",
            "second argument must have order dividing 2",
        ),
        (
            "show rvJacWeil[2](r,zero[RVJac(D)]);",
            "expected type RVJac(C)",
        ),
    ] {
        let problem = diagnostic(interpret(&format!("{source}{code}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
