mod common;

use common::{diagnostic, example, interpret, stdout};

const EXAMPLE: &str = include_str!("../examples/sheaf-cone-gluing.asl");

#[test]
fn sphere_link_changes_the_costalk_and_the_explicit_homotopy_is_invertible() {
    let output = stdout(example("sheaf-cone-gluing"));
    assert_eq!(
        output,
        concat!(
            "shhyperdim[Split__fiber, 1] = 1
",
            "shhyperdim[Split__fiber, 3] = 1
",
            "shhyperdim[Split__fiber, 4] = 1
",
            "shhyperdim[Attached__fiber, 1] = 1
",
            "shhyperdim[Attached__fiber, 3] = 0
",
            "shhyperdim[Attached__fiber, 4] = 0
",
            "shlocaldim[SplitAtTip, pt, relative, 1] = 1
",
            "shlocaldim[SplitAtTip, pt, relative, 3] = 1
",
            "shlocaldim[SplitAtTip, pt, relative, 4] = 1
",
            "shlocaldim[AttachedAtTip, pt, relative, 1] = 1
",
            "shlocaldim[AttachedAtTip, pt, relative, 3] = 0
",
            "shlocaldim[AttachedAtTip, pt, relative, 4] = 0
",
        )
    );
}

#[test]
fn three_apex_columns_kill_exactly_one_direction() {
    let tail = r"
        shconeglue Mixed[j,tip,F,3] {zero[ShC(F,3)],fundamental,shScale[2](fundamental)};
        let m0 = shGlueApex[Mixed,0];
        let m1 = shGlueApex[Mixed,1];
        let m2 = shGlueApex[Mixed,2];
        assert shD(m0) == zero[ShTot(Mixed__fiber,4)];
        assert shD(add(shScale[2](m1),shScale[-1](m2))) == zero[ShTot(Mixed__fiber,4)];
        assert shClass(shD(m1)) == zero[ShHyper(Mixed__fiber,4)];
        shhyperdim[Mixed__fiber,3];
        shhyperdim[Mixed__fiber,4];
    ";
    let output = stdout(interpret(&format!("{EXAMPLE}{tail}")));
    assert!(output.ends_with(
        "shhyperdim[Mixed__fiber, 3] = 2
shhyperdim[Mixed__fiber, 4] = 0
"
    ));
}

#[test]
fn genuinely_nonzero_vertical_complex_differential_and_prime_two() {
    let input = r"
        space U {a<b};
        space X {tip<a,a<b};
        let j = tab U -> X {a=>X.a,b=>X.b};
        sheaf A[U,2] {a:2,b:2} {a<b:[1,0,0,1]};
        sheaf B[U,2] {a:2,b:2} {a<b:[1,0,0,1]};
        shmap d[A,B] {a:[1,0,0,0],b:[1,0,0,0]};
        shcomplex K {A,d,B};
        let h = shTotEmbed[K,0](shTerm[A,[a],0]);
        let closed = add(shTotEmbed[K,1](shTerm[B,[a],1]),shTotEmbed[K,1](shTerm[B,[b],1]));
        shconegluecomplex Original[j,tip,K,1] {closed};
        shconegluecomplex Changed[j,tip,K,1] {add(closed,shD(h))};
        shconehomotopy H[Original,Changed]{h};
        let a = shGlueApex[Original,0];
        assert shGlueTail[Original](shD(a)) == closed;
        assert shD(shConeTransfer[H](a)) == shConeTransfer[H](shD(a));
        assert shConeUntransfer[H](shConeTransfer[H](a)) == a;
        shhyperdim[Original__fiber,1];
        shhyperdim[Original__fiber,2];
    ";
    assert_eq!(
        stdout(interpret(input)),
        "shhyperdim[Original__fiber, 1] = 1
shhyperdim[Original__fiber, 2] = 0
"
    );
}

#[test]
fn cone_and_homotopy_errors_certify_actual_input_conditions() {
    for (tail, reason) in [
        (
            "shconeglue N[j,tip,F,0]{shTerm[F,[a0],0]};",
            "column 0 is not closed",
        ),
        (
            "shconeglue N[j,tip,F,1]{fundamental};",
            "requires cochain degree 1, found 3",
        ),
        (
            "shconeglue N[j,u7,F,3]{fundamental};",
            "strictly below exactly",
        ),
        (
            "shconehomotopy Bad[Attached,Adjusted]{zero[ShC(F,2)]};",
            "does not satisfy D h",
        ),
        (
            "shconehomotopy Bad[Attached,Adjusted]{};",
            "requires 1 correcting cochains",
        ),
        ("show shGlueApex[Attached,1];", "outside its basis"),
        (
            "show shConeTransfer[H](shGlueApex[Adjusted,0]);",
            "source fiber",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{EXAMPLE}{tail}")));
        assert!(error.contains(reason), "{error}");
    }
}

#[test]
fn the_empty_link_and_degree_zero_homotopies_have_the_expected_contract() {
    let code = r"
        space Empty {};
        let open = tab Empty -> One {};
        sheaf F[Empty,3] {} {};
        let z = zero[ShC(F,0)];
        shconeglue G[open,pt,F,0] {z};
        shconeglue H[open,pt,F,0] {z};
        shconehomotopy Zero[G,H] {};
        let a = shGlueApex[G,0];
        assert shConeTransfer[Zero](a) == shGlueApex[H,0];
        assert shConeUntransfer[Zero](shConeTransfer[Zero](a)) == a;
        shhyperdim[G__fiber,0];
        shhyperdim[G__fiber,1];
    ";
    assert_eq!(
        stdout(interpret(code)),
        "shhyperdim[G__fiber, 0] = 1\nshhyperdim[G__fiber, 1] = 0\n"
    );
}
