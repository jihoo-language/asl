mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn genuine_class_groups_genera_and_oriented_period_bridge() {
    assert_eq!(
        stdout(example("hyperelliptic-k3-cm-class-group")),
        concat!(
            "rvCMClassOrder = 3\n",
            "rvCMSameGenus = true\n",
            "rvCMClassGroup[O] = discriminant=-23, order=3, invariant_factors=[3], genera=1\n",
            "rvCMClassGroup[C2] = discriminant=-20, order=2, invariant_factors=[2], genera=2\n",
            "rvCMClassGroup[V] = discriminant=-84, order=4, invariant_factors=[2, 2], genera=4\n",
            "rvCMClassGroup[N] = discriminant=-63, order=4, invariant_factors=[4], genera=2\n",
            "rvCMSameGenus = false\n",
            "rvCMSameGenus = true\n",
            "RVCMClass(O)[2, 1, 3]\n",
        ),
    );
    let source = include_str!("../examples/hyperelliptic-k3-cm-class-group.asl");
    let error = diagnostic(interpret(&format!(
        "{source} rvcmorder Raw[-92]; show rvK3CMClass[Raw,P];"
    )));
    assert!(
        error.contains("primitive form discriminant differs from the named order"),
        "{error}"
    );
}

#[test]
fn nonfundamental_composite_conductors_and_class_orders() {
    let output = stdout(interpret(
        "rvcmorder O[-368];
         rvcmorder G[-576];
         rvcmorder Unit[-3];
         rvCMClassGroup[O]; rvCMClassGroup[G]; rvCMClassGroup[Unit];
         let a=rvCMClass[O,3,2,31];
         let b=rvCMClass[G,5,2,29];
         rvCMClassOrder(a); rvCMClassOrder(b);
         assert rvCMPow[6](a)==rvCMPrincipal[O];
         assert rvCMPow[4](b)==rvCMPrincipal[G];",
    ));
    assert_eq!(
        output,
        concat!(
            "rvCMClassGroup[O] = discriminant=-368, order=6, invariant_factors=[6], genera=2\n",
            "rvCMClassGroup[G] = discriminant=-576, order=8, invariant_factors=[2, 4], genera=4\n",
            "rvCMClassGroup[Unit] = discriminant=-3, order=1, invariant_factors=[], genera=1\n",
            "rvCMClassOrder = 6\nrvCMClassOrder = 4\n",
        ),
    );
}

#[test]
fn large_coefficients_cancel_and_shared_leading_coefficients_compose() {
    let min = i128::MIN;
    let k = 1i128 << 62;
    let b = 1 + 4 * k;
    let c = 3 + k + 2 * k * k;
    let source = format!(
        "rvcmorder Huge[{min}];
         let one=rvCMPrincipal[Huge];
         assert one==rvCMClass[Huge,1,0,42535295865117307932921825928971026432];
         let x=rvCMClass[Huge,3,2,14178431955039102644307275309657008811];
         let xx=rvCMClass[Huge,9,-4,4726143985013034214769091769885669604];
         assert rvCMCompose(x,x)==xx;
         assert rvCMCompose(x,rvCMInverse(x))==one;
         rvcmorder O[-23];
         let alpha=rvCMClass[O,2,1,3];
         assert rvCMClass[O,2,{b},{c}]==alpha;
         assert rvCMPow[{min}](alpha)==alpha;
         rvCMClassOrder(alpha);"
    );
    assert_eq!(stdout(interpret(&source)), "rvCMClassOrder = 3\n");
}

#[test]
fn rejects_invalid_discriminants_imprimitive_and_nominal_mixed_classes() {
    for (command, expected) in [
        ("rvcmorder Invalid[0];", "must be negative"),
        ("rvcmorder Invalid[16];", "must be negative"),
        ("rvcmorder Invalid[-2];", "congruent to 0 or 1 modulo 4"),
        ("rvcmorder O[-23];", "already defined"),
        ("show rvCMPrincipal[Missing];", "unknown CM ring order"),
        ("show rvCMClass[O,-1,1,-6];", "positive definite"),
        ("show rvCMClass[O,1,0,5];", "discriminant differs"),
        ("show rvCMClass[U,4,2,6];", "requires a primitive form"),
        (
            "show rvCMCompose(rvCMPrincipal[O],rvCMPrincipal[Twin]);",
            "expected type RVCMClass(O)",
        ),
        (
            "rvCMSameGenus(rvCMPrincipal[O],rvCMPrincipal[Twin]);",
            "expected type RVCMClass(O)",
        ),
        (
            "rvCMClassGroup[Huge];",
            "class number enumeration exceeds the computation limit",
        ),
    ] {
        let error = diagnostic(interpret(&format!(
            "rvcmorder O[-23]; rvcmorder Twin[-23]; rvcmorder U[-92]; rvcmorder Huge[{}]; {command}",
            i128::MIN,
        )));
        assert!(error.contains(expected), "{error}");
    }
}

#[test]
fn canonical_class_type_enumerates_finitely() {
    let output = stdout(interpret(
        r"rvcmorder O[-23]; show (\x:RVCMClass(O)=>rvCMInverse(x));",
    ));
    assert!(output.contains("RVCMClass(O)[1, 1, 6]"), "{output}");
    assert!(output.contains("RVCMClass(O)[2, -1, 3]"), "{output}");
    assert!(output.contains("RVCMClass(O)[2, 1, 3]"), "{output}");
}
