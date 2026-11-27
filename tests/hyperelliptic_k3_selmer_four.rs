mod common;

use common::{diagnostic, example, interpret, stdout};

const F23: &str = include_str!("../examples/hyperelliptic-k3-cassels-pairing.asl");
const THIRD: &str = r"
rvk3cassels Third[Sum]{
  0{[0],[3],[9]}{[9],[4],[22]}{[1],[1]},
  1{[2]}{[1],[17]}{[1]},
  2{[10]}{}{[1]}
};
";

fn f23() -> &'static str {
    F23.split_once("rvk3selmer4 R[M]").unwrap().0
}

#[test]
fn f23_pairing_identifies_exact_four_selmer_image_and_sha_four() {
    let output = stdout(example("hyperelliptic-k3-cassels-pairing"));
    let report = output
        .lines()
        .find(|line| line.starts_with("rvK3Selmer4[R]"))
        .unwrap();
    assert!(report.contains("dim Sel2=6, dim E(K)/2=4, quotient source masks=[16962, 33441], Cassels matrix=[[0, 1], [1, 0]], rank=2"), "{report}");
    assert!(report.contains("symplectic pairs=[(16962, 33441)], quotient radical lifts=[], dim image Sel4->Sel2=4, image basis masks=[55489, 18754, 18528, 14848]"), "{report}");
    assert!(report.contains("dim Sha[2] intersect 2Sha=0, Sha[4]=(Z/2)^2, Sel4 order=256, torsion=(Z/2)^2 certified by I0* at t=[0] and I4 at t=[21]"));
    assert!(output.lines().any(|line| line == "PrimeScalar(2){1}"));
    assert!(output
        .lines()
        .last()
        .unwrap()
        .starts_with("RVK3PencilKummer2(P)"));
}

#[test]
fn actual_lifts_and_pairing_are_independent_of_quotient_basis() {
    let source = format!(
        r"{}{THIRD}
        rvk3selmer4 R[M]{{T,Second}};
        rvk3selmer4 Changed[M]{{Third,T}};
        let zeroClass=zero[RVK3PencilKummer2(P)];
        let c=add(a,b);
        let mwq=rvK3PencilDescent2(q);
        let torsion=rvK3PencilDescent2(rvK3PencilTorsion2At[P,0]);
        let e=rvK3Selmer4QuotientBasis[Changed,0];
        let f=rvK3Selmer4QuotientBasis[Changed,1];
        assert rvK3Selmer4Pair[Changed](e,f)==rvK3Selmer4Pair[R](a,b);
        assert rvK3Selmer4Pair[Changed](add(a,mwq),add(b,torsion))==rvK3Selmer4Pair[R](a,b);
        assert rvK3Selmer4Pair[Changed](c,b)==rvK3Selmer4Pair[R](a,b);
        assert rvK3Selmer4Pair[Changed](a,a)==rvK3Selmer4Pair[R](zeroClass,b);
        assert rvK3Selmer4Pair[R](a,add(b,mwq))==add(rvK3Selmer4Pair[R](a,b),rvK3Selmer4Pair[R](a,mwq));
        assert rvK3Selmer4Contains[R](torsion)==S.hi;
        assert rvK3Selmer4Contains[R](mwq)==S.hi;
        assert rvK3Selmer4Contains[Changed](delta)==S.hi;
        assert rvK3Selmer4Contains[Changed](c)==S.lo;
        assert rvK3Selmer4Contains[Changed](add(a,mwq))==S.lo;
        assert rvK3Selmer4Contains[Changed](rvK3PencilKummer2[P]{{[5]}}{{[5]}}{{[1]}})==S.lo;
        assert rvK3Selmer4ImageBasis[R,0]==rvK3Selmer4ImageBasis[Changed,0];
        assert rvK3Selmer4ImageBasis[R,1]==rvK3Selmer4ImageBasis[Changed,1];
        assert rvK3Selmer4ImageBasis[R,2]==rvK3Selmer4ImageBasis[Changed,2];
        assert rvK3Selmer4ImageBasis[R,3]==rvK3Selmer4ImageBasis[Changed,3];
        show rvK3Selmer4Pair[Changed](a,mwq);
        show rvK3Selmer4Pair[Changed](a,b);
        ",
        f23()
    );
    let output = stdout(interpret(&source));
    assert!(output.ends_with("PrimeScalar(2){0}\nPrimeScalar(2){1}\n"));
}

#[test]
fn zero_sha_needs_no_tangent_witnesses_and_the_image_is_all_selmer_two() {
    let base = include_str!("../examples/hyperelliptic-k3-pencil-regulator.asl")
        .split_once("rvK3PencilSectionHeight[P]")
        .unwrap()
        .0;
    let source = format!(
        r"{base}
        rvk3pencilmw2 M[P]{{p,q}};
        rvk3selmer4 R[M]{{}};
        rvK3Selmer4[R];
        let a=rvK3PencilSelmerBasis[P,0];
        let b=rvK3PencilSelmerBasis[P,1];
        assert rvK3Selmer4Contains[R](a)==S.hi;
        assert rvK3Selmer4Contains[R](add(a,b))==S.hi;
        assert rvK3Selmer4Contains[R](rvK3Selmer4ImageBasis[R,3])==S.hi;
        show rvK3Selmer4Pair[R](a,b);
        "
    );
    let output = stdout(interpret(&source));
    assert!(output
        .contains("dim Sel2=4, dim E(K)/2=4, quotient source masks=[], Cassels matrix=[], rank=0"));
    assert!(output.contains("dim image Sel4->Sel2=4"));
    assert!(output.contains("dim Sha[2] intersect 2Sha=0, Sha[4]=0, Sel4 order=64"));
    assert!(output.ends_with("PrimeScalar(2){0}\n"));
}

#[test]
fn incomplete_dependent_or_wrong_typed_inputs_are_rejected() {
    let base = format!(
        r"{}
        rvk3selmer4 R[M]{{T,Second}};
        rvk3pencil Other[B,14,0];
        let pother=rvK3PencilSectionPoint[Other]{{[0],[18],[2],[5]}}{{[1]}}{{[0],[0],[18],[11],[10],[5],[6]}}{{[1]}};
        let qother=rvK3PencilSectionPoint[Other]{{[0],[18],[4],[5]}}{{[1]}}{{[0],[0],[14],[7],[18],[5],[6]}}{{[1]}};
        rvk3pencilmw2 Mother[Other]{{pother,qother}};
        let bad=rvK3PencilKummer2[P]{{[5]}}{{[5]}}{{[1]}};
        ",
        f23()
    );
    for (operation, expected) in [
        ("rvk3selmer4 X[M]{T};", "requires 2 Cassels witnesses"),
        (
            "rvk3selmer4 X[M]{T,T};",
            "dependent modulo the global Mordell-Weil image",
        ),
        (
            "rvk3selmer4 X[Mother]{T,Second};",
            "must use the certificate's named pencil",
        ),
        ("rvk3selmer4 R[M]{T,Second};", "already defined"),
        ("rvK3Selmer4[Missing];", "unknown four-Selmer certificate"),
        (
            "show rvK3Selmer4Pair[R](a,bad);",
            "requires two everywhere locally soluble classes",
        ),
        (
            "show rvK3Selmer4Contains[R](zero[RVK3PencilKummer2(Other)]);",
            "expected type RVK3PencilKummer2(P)",
        ),
        ("show rvK3Selmer4QuotientBasis[R,2];", "outside dimension 2"),
        ("show rvK3Selmer4ImageBasis[R,4];", "outside dimension 4"),
    ] {
        let error = diagnostic(interpret(&format!("{base}\n{operation}")));
        assert!(error.contains(expected), "{operation}: {error}");
    }
}
