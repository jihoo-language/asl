mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const MARKED: &str = r"
    rvunram F[3,1]{0,1};
    rvunram W[3,1]{1,0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
    rvjacbase L[C,W];
    let zeroW=rvUnramElem[W]{};
    let oneW=rvUnramElem[W]{1};
    let alpha=rvUnramElem[W]{0,1};
    let e1=rvJacPoint[L](zeroW,zeroW);
    let f1=rvJacPoint[L](oneW,zeroW);
    let e2=add(add(e1,f1),rvJacPoint[L](neg(oneW),zeroW));
    let f2=add(add(e1,f1),rvJacPoint[L](alpha,zeroW));
    rvjactorsionbasis B[L,2](neg(oneW)){e1,e2,f1,f2};
";

#[test]
fn anticommutators_recover_the_full_integral_gram_matrix() {
    let report = stdout(interpret(&format!("{MARKED} rvK3LatticeGram[B];")));
    let gram: Vec<i128> = report
        .split_once("gram=")
        .unwrap()
        .1
        .split(|c: char| c != '-' && !c.is_ascii_digit())
        .filter(|entry| !entry.is_empty())
        .map(|entry| entry.parse().unwrap())
        .collect();
    assert_eq!(gram.len(), 22 * 22);
    assert!((0..22).any(|i| (i + 1..22).any(|j| gram[22 * i + j] != 0)));
    let mut source = format!("{MARKED} let unit=rvCliffordOne[B];");
    for i in 0..22 {
        writeln!(
            source,
            "let g{i}=rvCliffordVector(rvK3LatticeBasis[B,{i}]);"
        )
        .unwrap();
    }
    for i in 0..22 {
        for j in i..22 {
            let entry = gram[22 * i + j];
            writeln!(
                source,
                "assert add(rvCliffordMul(g{i},g{j}),rvCliffordMul(g{j},g{i})) == rvCliffordScale[{entry},1](unit);"
            ).unwrap();
        }
    }
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn reflections_volume_and_half_spin_idempotents() {
    let output = stdout(example("hyperelliptic-kummer-clifford"));
    assert_eq!(
        output,
        concat!(
            "rvCliffordCoeff[0] = -1\n",
            "rvCliffordCoeff[4194303] = 1/8\n",
            "rvCliffordCoeff[4194303] = 1/16\n",
            "rvCliffordCoeff[4194303] = -1/16\n",
        ),
    );
}

#[test]
fn sparse_associativity_reversion_and_exact_rational_cancellation() {
    let huge = i128::MAX - 1;
    let source = format!(
        "{MARKED}
         let a=rvCliffordVector(rvK3LatticeBasis[B,0]);
         let b=rvCliffordVector(rvK3LatticeBasis[B,1]);
         let c=rvCliffordVector(rvK3LatticeBasis[B,2]);
         let unit=rvCliffordOne[B];
         let v=rvCliffordVolume[B];
         let x=add(unit,rvCliffordMul(a,b));
         let y=add(v,rvCliffordMul(b,c));
         let odd=add(a,c);
         assert rvCliffordMul(rvCliffordMul(x,y),odd)==rvCliffordMul(x,rvCliffordMul(y,odd));
         assert rvCliffordReverse(rvCliffordMul(x,y))==rvCliffordMul(rvCliffordReverse(y),rvCliffordReverse(x));
         assert rvCliffordReverse(rvCliffordReverse(x))==x;
         assert rvCliffordMul(x,add(a,c))==add(rvCliffordMul(x,a),rvCliffordMul(x,c));
         assert add(odd,neg(odd))==zero[RVK3Clifford(B,1)];
         assert rvCliffordMul(x,zero[RVK3Clifford(B,1)])==zero[RVK3Clifford(B,1)];
         assert rvCliffordScale[1,2](add(rvCliffordScale[1,3](unit),rvCliffordScale[5,3](unit)))==unit;
         let wide=rvK3Lattice[B]{{[{huge},1,0,0,0,{huge}],[0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,1]}};
         let actor=rvCliffordScale[-1,{denominator}](unit);
         assert rvCliffordAct(actor,wide)==wide;
         let tiny=rvCliffordMul(actor,actor);
         assert rvCliffordAct(tiny,wide)==wide;
         rvCliffordCoeff[1](rvCliffordVector(wide));
         rvCliffordCoeff[2](rvCliffordVector(wide));
         rvCliffordCoeff[0](tiny);",
        denominator = u128::MAX,
    );
    assert_eq!(
        stdout(interpret(&source)),
        format!(
            "rvCliffordCoeff[1] = {}\nrvCliffordCoeff[2] = 0\nrvCliffordCoeff[0] = 1/115792089237316195423570985008687907852589419931798687112530834793049593217025\n",
            huge / 2
        )
    );
}

#[test]
fn checks_normalizer_conditions_and_integral_gluing() {
    let mut setup = format!("{MARKED} let unit=rvCliffordOne[B];");
    for i in 0..6 {
        writeln!(
            setup,
            "let node{i}=rvK3FromKummerDiv(rvKummerNodeDiv[B,{i}]); let v{i}=rvCliffordVector(node{i});"
        )
        .unwrap();
    }
    setup.push_str(
        "let v4form=rvCliffordMul(rvCliffordMul(v0,v1),rvCliffordMul(v2,v3));
         let v6form=rvCliffordMul(v4form,rvCliffordMul(v4,v5));
         let glue=rvK3Lattice[B]{[0,0,0,0,0,0],[0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]};",
    );
    for (command, expected) in [
        (
            "show rvCliffordAct(zero[RVK3Clifford(B,0)],node0);",
            "times its reverse to be a nonzero scalar",
        ),
        (
            "show rvCliffordAct(add(unit,v4form),node0);",
            "times its reverse to be a nonzero scalar",
        ),
        (
            "show rvCliffordAct(add(unit,v6form),node0);",
            "conjugation did not return a vector",
        ),
        (
            "show rvCliffordAct(add(v0,rvCliffordScale[2,1](v1)),node0);",
            "image has a nonintegral K3 lattice coordinate",
        ),
        (
            "show rvCliffordAct(add(v0,neg(v1)),glue);",
            "affine gluing condition",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}{command}")));
        assert!(error.contains(expected), "{error}");
    }
}

#[test]
fn checks_parity_nominal_types_and_numeric_bounds() {
    let max = i128::MAX;
    let overflow = format!(
        "show rvCliffordScale[{max},1](rvCliffordScale[{max},1](rvCliffordScale[{max},1](unit)));"
    );
    let setup = format!(
        "{MARKED}
         rvjactorsionbasis Twin[L,2](neg(oneW)){{e1,e2,f1,f2}};
         let unit=rvCliffordOne[B];
         let x=rvK3LatticeBasis[B,0];
         let g=rvCliffordVector(x);"
    );
    for (command, expected) in [
        ("show zero[RVK3Clifford(B,2)];", "parity must be 0 or 1"),
        ("show add(unit,g);", "RVK3Clifford(B,0)"),
        (
            "show rvCliffordMul(g,rvCliffordOne[Twin]);",
            "RVK3Clifford(B,0)",
        ),
        (
            "show rvCliffordAct(g,rvK3LatticeBasis[Twin,0]);",
            "RVK3Lattice(B)",
        ),
        (
            "show rvCliffordScale[1,0](g);",
            "denominator must be positive",
        ),
        ("rvCliffordCoeff[4194304](g);", "22 bits"),
        (
            &overflow,
            "Clifford rational coefficient exceeds the supported exact integer range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}{command}")));
        assert!(error.contains(expected), "{error}");
    }
}

#[test]
fn large_sparse_product_stops_at_the_statement_work_limit() {
    let mut source = format!(
        "{MARKED}
         let unit=rvCliffordOne[B];
         let u=rvCliffordVector(rvK3Lattice[B]{{[2,0,0,0,0,2],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}});
         let v=rvCliffordVector(rvK3Lattice[B]{{[2,0,0,0,0,-2],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}});
         let g0=add(unit,rvCliffordMul(u,v));"
    );
    for i in 0..16 {
        writeln!(
            source,
            "let exc{i}=rvCliffordVector(rvK3FromKummerDiv(rvKummerNodeDiv[B,{i}]));"
        )
        .unwrap();
    }
    for i in 0..8 {
        let left = 2 * i;
        let right = left + 1;
        let next = i + 1;
        writeln!(
            source,
            "let g{next}=rvCliffordMul(g{i},add(unit,rvCliffordMul(exc{left},exc{right})));"
        )
        .unwrap();
    }
    source.push_str("show rvCliffordMul(g8,g8);");
    let error = diagnostic(interpret(&source));
    assert!(
        error.contains("Clifford operation exceeds the sparse work or storage limit"),
        "{error}"
    );
}
