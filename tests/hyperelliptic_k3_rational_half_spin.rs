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
fn the_volume_projectors_and_glued_root_actions_are_exact() {
    let output = stdout(example("hyperelliptic-k3-rational-half-spin"));
    assert_eq!(
        output,
        concat!(
            "rvK3SpinCoeff[0] = 2/3\n",
            "rvK3SpinCoeff[1023] = -5/7\n",
            "rvK3SpinCoeff[1023] = 1\n",
        )
    );
}

#[test]
fn signed_permutation_coordinates_and_wide_rational_cancellation() {
    let (n, d) = (i128::MAX, i128::MAX - 2);
    let source = format!(
        "{MARKED}
         let s=rvK3SpinorBasis[B,0,0];
         let unit=rvCliffordOne[B];
         let wedge=rvCliffordVector(rvK3Lattice[B]{{[2,0,0,0,0,2],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}});
         let first=rvCliffordVector(rvK3FromKummerDiv(rvKummerNodeDiv[B,0]));
         let second=rvCliffordVector(rvK3FromKummerDiv(rvKummerNodeDiv[B,8]));
         assert rvK3SpinAct(wedge,s)==rvK3SpinorBasis[B,1,1];
         assert rvK3SpinAct(wedge,rvK3SpinorBasis[B,1,1])==rvK3SpinorScale[-2,1](s);
         assert rvK3SpinAct(first,s)==rvK3SpinorBasis[B,1,72];
         assert rvK3SpinAct(second,s)==rvK3SpinorBasis[B,1,128];
         assert rvK3SpinAct(zero[RVK3Clifford(B,1)],s)==zero[RVK3Spinor(B,1)];
         assert rvK3SpinAct(wedge,zero[RVK3Spinor(B,1)])==zero[RVK3Spinor(B,0)];
         let a=rvCliffordScale[{n},{d}](unit);
         let actor=rvCliffordMul(a,a);
         let wide=rvK3SpinorScale[{d},{n}](rvK3SpinorScale[{d},{n}](s));
         assert rvK3SpinAct(actor,wide)==s;
         assert add(wide,neg(wide))==zero[RVK3Spinor(B,0)];
         rvK3SpinCoeff[10](wide);
         rvK3SpinCoeff[0](rvK3SpinAct(actor,wide));"
    );
    assert_eq!(
        stdout(interpret(&source)),
        "rvK3SpinCoeff[10] = 0\nrvK3SpinCoeff[0] = 1\n"
    );
}

#[test]
fn bounds_and_nominal_type_errors() {
    let setup = format!(
        "{MARKED}
         rvjactorsionbasis Twin[L,2](neg(oneW)){{e1,e2,f1,f2}};
         let s=rvK3SpinorBasis[B,0,0];
         let t=rvK3SpinorBasis[B,1,0];"
    );
    for (command, expected) in [
        (
            "show zero[RVK3Spinor(B,2)];",
            "K3 spinor chirality must be 0 or 1",
        ),
        (
            "show rvK3SpinorBasis[B,0,1024];",
            "K3 spinor index must be between 0 and 1023",
        ),
        (
            "rvK3SpinCoeff[1024](s);",
            "K3 spinor index must be between 0 and 1023",
        ),
        ("show add(s,t);", "expected type RVK3Spinor(B,0)"),
        (
            "show rvK3SpinAct(rvCliffordOne[Twin],s);",
            "expected type RVK3Spinor(Twin,0)",
        ),
        ("show rvK3SpinAct(s,t);", "an RVK3Clifford element"),
        (
            "show rvK3SpinorScale[1,0](s);",
            "K3 spinor scalar denominator must be positive",
        ),
        (
            "rvK3SpinCoeff[1](rvCliffordOne[B]);",
            "an RVK3Spinor vector",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}{command}")));
        assert!(error.contains(expected), "{error}");
    }
}

#[test]
fn sparse_action_respects_the_statement_work_limit() {
    let mut source = format!("{MARKED} let unit=rvCliffordOne[B];");
    for i in 0..16 {
        writeln!(
            source,
            "let exc{i}=rvCliffordVector(rvK3FromKummerDiv(rvKummerNodeDiv[B,{i}]));"
        )
        .unwrap();
    }
    source.push_str("let g0=unit;");
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
    for i in 0..1024 {
        writeln!(source, "let s0x{i}=rvK3SpinorBasis[B,0,{i}];").unwrap();
    }
    for layer in 1..=10 {
        let child = layer - 1;
        for i in 0..(1024 >> layer) {
            let left = 2 * i;
            let right = left + 1;
            writeln!(
                source,
                "let s{layer}x{i}=add(s{child}x{left},s{child}x{right});"
            )
            .unwrap();
        }
    }
    source.push_str("show rvK3SpinAct(g8,s10x0);");
    let error = diagnostic(interpret(&source));
    assert!(
        error.contains("K3 spinor operation exceeds the sparse work or storage limit"),
        "{error}"
    );
}
