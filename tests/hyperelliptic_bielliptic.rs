mod common;

use common::{diagnostic, example, interpret, stdout};

const F7: &str = r"
    rvunram F[7,1]{0,1};
    let z=rvUnramElem[F]{0};
    let o=rvUnramElem[F]{1};
    let two=rvUnramElem[F]{2};
    let four=rvUnramElem[F]{4};
    let six=rvUnramElem[F]{6};
    rvhyperelliptic C[F]{six,z,four,z,six,z,o};
";

const F9: &str = r"
    rvunram F[3,1]{1,0,1};
    let z=rvUnramElem[F]{0};
    let o=rvUnramElem[F]{1};
    let two=rvUnramElem[F]{2};
    let i=rvUnramElem[F]{0,1};
    let i2=neg(i);
    let a=rvUnramElem[F]{1,1};
    let b=rvUnramElem[F]{1,2};
    let c=rvUnramElem[F]{2,1};
    // x=i/(t-1) gives H(u)=i+u+(1+i)u²+(1+i)u³.
    rvhyperelliptic C[F]{a,z,neg(a),z,o,o};
    rvbielliptic C[Plus,Minus](z,i,o,two);
";

fn assert_weil_product(text: &str, plus: [i128; 3], minus: [i128; 3]) {
    let mut product = [0; 5];
    for (i, a) in plus.into_iter().enumerate() {
        for (j, b) in minus.into_iter().enumerate() {
            product[i + j] += a * b;
        }
    }
    for (name, polynomial) in [
        ("C", product.as_slice()),
        ("Plus", plus.as_slice()),
        ("Minus", minus.as_slice()),
    ] {
        let expected = format!("rvHWWWeil[{name}] = {polynomial:?}");
        assert!(
            text.lines().any(|line| line == expected),
            "{expected}\n{text}"
        );
    }
}

#[test]
fn even_sextic_splits_with_two_quotient_maps() {
    let text = stdout(example("hyperelliptic-bielliptic"));
    assert!(text.starts_with(concat!(
        "RVJac(Plus){u=[[5], [1]], v=[[2]]}\n",
        "RVJac(Minus){u=[[4], [1]], v=[[2]]}\n",
    )));
    assert_weil_product(&text, [1, 0, 7], [1, 0, 7]);
}

#[test]
fn quintic_mobius_chart_over_f9_preserves_twists_and_both_involutions() {
    let text = stdout(interpret(&format!(
        r"{F9}
        let p=rvBiellipticPoint[Plus](o,i);
        let q=rvBiellipticPoint[Minus](o,i);
        assert p==rvJacPoint[Plus](b,a);
        assert q==rvJacPoint[Minus](two,b);
        assert rvBiellipticPoint[Plus](b,c)==p;
        assert rvBiellipticPoint[Minus](b,c)==neg(q);
        assert rvBiellipticPoint[Plus](b,neg(c))==neg(p);
        assert rvBiellipticPoint[Minus](b,neg(c))==q;
        assert rvBiellipticPoint[Plus](o,neg(i))==neg(p);
        assert rvBiellipticPoint[Minus](o,neg(i))==neg(q);

        // t=0 is the minus pole. Inverting t exchanges the two factors.
        assert rvBiellipticPoint[Minus](i2,c)==zero[RVJac(Minus)];
        assert rvBiellipticPoint[Plus](i2,c)==rvJacPoint[Plus](z,two);
        rvbielliptic C[ReversePlus,ReverseMinus](i,z,two,o);
        assert rvBiellipticPoint[ReversePlus](i2,c)==zero[RVJac(ReversePlus)];
        assert rvBiellipticPoint[ReverseMinus](i2,c)==rvJacPoint[ReverseMinus](z,two);

        rvHWWPoints[C,1];
        rvHWWPoints[C,2];
        rvHWWPoints[Plus,1];
        rvHWWPoints[Minus,1];
        rvHWWVerify[Plus];
        rvHWWVerify[Minus];
        rvHWWWeil[C];
        rvHWWWeil[Plus];
        rvHWWWeil[Minus];
    "
    )));
    for expected in [
        "rvHWWPoints[C, 1] = 10, trace=0",
        "rvHWWPoints[C, 2] = 110, trace=-28",
        "rvHWWPoints[Plus, 1] = 8, trace=2",
        "rvHWWPoints[Minus, 1] = 12, trace=-2",
        "rvHWWVerify[Plus] mod 3 = unit=[1, 1], quotient=[1, 0]",
        "rvHWWVerify[Minus] mod 3 = unit=[1, 2], quotient=[1, 0]",
    ] {
        assert!(
            text.lines().any(|line| line == expected),
            "{expected}\n{text}"
        );
    }
    assert_weil_product(&text, [1, -2, 9], [1, 2, 9]);
}

#[test]
fn rejects_invalid_charts_and_conflicting_or_wrong_field_names() {
    for (input, expected) in [
        ("rvbielliptic C[Plus,Minus](o,o,o,o);", "determinant must be nonzero"),
        ("rvbielliptic C[Plus,Minus](o,o,z,o);", "must give an even sextic"),
        ("rvbielliptic C[Plus,Minus](o,two,z,o);", "nonzero constant and sextic coefficients"),
        ("rvbielliptic C[Plus,Minus](two,z,o,o);", "nonzero constant and sextic coefficients"),
        ("rvbielliptic C[Plus,Plus](o,z,z,o);", "target names must be distinct"),
        ("rvbielliptic C[Plus,C](o,z,z,o);", "curve is already defined"),
        (
            "rvhyperelliptic Legacy[7]{1,1,0,1}; rvbielliptic C[Legacy,Minus](o,z,z,o);",
            "curve is already defined",
        ),
        (
            "rvhyperelliptic E[F]{o,o,z,o}; rvbielliptic E[Plus,Minus](o,z,z,o);",
            "smooth genus-two curve",
        ),
        (
            "rvunram W[7,2]{0,1}; let w=rvUnramElem[W]{1}; rvhyperelliptic D[W]{w,w,w,w,w,w}; rvbielliptic D[Plus,Minus](w,w,w,w);",
            "precision-one field",
        ),
        (
            "rvunram G[7,1]{0,1}; rvbielliptic C[Plus,Minus](rvUnramElem[G]{1},z,z,o);",
            "bielliptic chart coefficient: expected type RVUnram(F), found RVUnram(G)",
        ),
        (
            "rvrichelot D[C]{[1],[],[1]}{[2],[],[1]}{[3],[],[1]};",
            "product of elliptic curves over the algebraic closure",
        ),
    ] {
        let problem = diagnostic(interpret(&format!("{F7}{input}")));
        assert!(problem.contains(expected), "{problem}");
    }
}

#[test]
fn validates_source_points_before_mapping_poles_and_keeps_nominal_targets() {
    for (input, expected) in [
        (
            "show rvBiellipticPoint[Minus](i2,z);",
            "point does not lie on the bielliptic source curve",
        ),
        (
            "show rvBiellipticPoint[Plus](o,z);",
            "point does not lie on the bielliptic source curve",
        ),
        (
            "show rvBiellipticPoint[C](o,i);",
            "has no stored bielliptic quotient",
        ),
        (
            "rvunram G[3,1]{1,0,1}; show rvBiellipticPoint[Plus](o,rvUnramElem[G]{0,1});",
            "expected type RVUnram(F), found RVUnram(G)",
        ),
        (
            "show add(rvBiellipticPoint[Plus](o,i),rvBiellipticPoint[Minus](o,i));",
            "expected type RVJac(Plus)",
        ),
    ] {
        let problem = diagnostic(interpret(&format!("{F9}{input}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
