mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn odd_prime_tensor_shortens_the_long_bar_and_keeps_the_secondary_class() {
    assert_eq!(
        stdout(example("mixed-tensor")),
        concat!(
            "mixbarrank[T, 10, 8] = 1\n",
            "mixbarrank[T, 10, 6] = 0\n",
            "mixbarrank[T, 5, 3] = 1\n",
            "mixbarrank[T, 5, 1] = 0\n",
            "mixbars[T] = free [], torsion [(5, 2), (10, 2)]\n",
        )
    );
}

#[test]
fn free_finite_and_free_free_cycles_have_the_expected_orders() {
    let fixture = include_str!("../examples/mixed-tensor.asl");
    let source = format!(
        r"{fixture}
        mixedcomplex G[5] {{g:-1}} {{}} {{}};
        mixedcomplex H[5] {{h:-2}} {{}} {{}};
        mixtensor GR[G,Right];
        mixtensor GH[G,H];
        let g = mixBarRep[G,0];
        let h = mixBarRep[H,0];
        assert mixTensorCycle[GH,0,0] == mixTensor[GH](g,h);
        assert mixD(mixTensorCycle[GH,0,0]) == zero[MixedU(GH,-4)];
        assert mixTensorNull[GR,0,0] == pscale[4](mixTensor[GR](g,Y));
        assert mixD(mixTensorNull[GR,0,0]) == mixU[3](mixTensorCycle[GR,0,0]);
        mixbars[GR];
        mixbars[GH];
        mixbarperiodicdim[GH,-99];
        "
    );
    assert!(stdout(interpret(&source)).ends_with(concat!(
        "mixbars[GR] = free [], torsion [(5, 3)]\n",
        "mixbars[GH] = free [-3], torsion []\n",
        "mixbarperiodicdim[GH, -99] = 1\n",
    )));
    for (expr, expected) in [
        ("mixTensorNull[GH,0,0]", "requires at least one torsion bar"),
        ("mixTensorTor[GR,0,0]", "requires two torsion bars"),
        ("mixTensorTorNull[GR,0,0]", "requires two torsion bars"),
    ] {
        let e = diagnostic(interpret(&format!("{source} show {expr};")));
        assert!(e.contains(expected), "{expr}: {e}");
    }
}

#[test]
fn kunneth_witnesses_on_actual_imported_exterior_hochschild_fibers() {
    let source = r"
        rvideal J[[0]] {rvmon[[0]]{1}};
        rvHHMixed Odd[J,[1]];
        rvHHMixed Even[J,[2]];
        mixtensor OO[Odd,Odd];
        mixtensor OE[Odd,Even];
        assert mixD(mixTensorTor[OO,0,0]) == zero[MixedU(OO,2)];
        assert mixD(mixTensorTorNull[OO,0,0]) == mixU[1](mixTensorTor[OO,0,0]);
        assert mixD(mixTensorNull[OO,0,0]) == mixU[1](mixTensorCycle[OO,0,0]);
        assert mixD(mixTensorNull[OE,0,1]) == mixU[1](mixTensorCycle[OE,0,1]);
        mixbars[OO];
        mixbars[OE];
    ";
    assert_eq!(
        stdout(interpret(source)),
        "mixbars[OO] = free [], torsion [(3, 1), (4, 1)]\nmixbars[OE] = free [], torsion [(5, 1), (6, 1)]\n",
    );
}

#[test]
fn tensor_types_check_field_endpoint_order_sizes_and_signed_sums() {
    let source = include_str!("../examples/mixed-tensor.asl");
    for (statement, expected) in [
        ("show mixTensor[T](y,x);", "factor complexes in order"),
        ("show mixTensor[T](mixU[0](x),y);", "same form"),
        ("show mixSwap[T](raw);", "factors reversed"),
        ("show mixSwap[Left](raw);", "defined by mixtensor"),
        ("show mixTensorCycle[T,1,0];", "index is out of range"),
        (
            "mixedcomplex W[3] {x:0} {} {}; mixtensor U[Left,W];",
            "same prime field",
        ),
        (
            "mixedcomplex Q[5] {x:9223372036854775807} {} {}; mixtensor U[Left,Q];",
            "degree is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{source} {statement}")));
        assert!(error.contains(expected), "{statement}: {error}");
    }
    let gens = (0..300)
        .map(|i| format!("x{i}:0"))
        .collect::<Vec<_>>()
        .join(",");
    let oversized = format!("mixedcomplex V[2] {{{gens}}} {{}} {{}}; mixtensor X[V,V];");
    assert!(diagnostic(interpret(&oversized))
        .contains("mixed tensor complex exceeds the storage limit"));
}
