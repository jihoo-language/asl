mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn a_signed_scrambled_complex_produces_free_and_finite_bars_with_actual_witnesses() {
    assert_eq!(
        stdout(example("mixed-bars")),
        concat!(
            "mixbardim[M, -2] = 2\n",
            "mixbarrank[M, 2, -2] = 1\n",
            "mixbarrank[M, 2, -4] = 0\n",
            "mixbarperiodicdim[M, 0] = 1\n",
            "mixbars[M] = free [-2, 1], torsion [(-1, 4), (2, 3)]\n",
        )
    );
}

#[test]
fn actual_shape_of_a_torsion_witness_over_small_and_large_primes() {
    for p in [2, 3, 5, 4_294_967_291u64] {
        let source = format!(
            r"mixedcomplex C[{p}] {{ q:-3, x:-1, z:-2, y:0 }}
            {{ x=>[(z,-1)] }}
            {{ q=>[(z,1)], x=>[(y,-1)] }};
            let rep = mixBarRep[C,0];
            let witness = mixBarNull[C,0];
            "
        );
        let source = format!(
            r"{source}
            let q = mixU[0](mixchain[C,-3][(q,-1)]);
            let x = mixU[1](mixchain[C,-1][(x,-1)]);
            assert witness == add(q,x);
            assert rep == mixU[0](mixchain[C,0][(y,1)]);
            assert mixD(witness) == mixU[2](rep);
            assert mixD(rep) == zero[MixedU(C,-1)];
            assert mixD(mixPageClass[2,0](witness)) == mixPageClass[2,2](mixU[2](rep));
            mixbarrank[C,0,-2];
            mixbarrank[C,0,-4];
            mixbarperiodicdim[C,0];
            mixbars[C];
            "
        );
        assert_eq!(
            stdout(interpret(&source)),
            "mixbarrank[C, 0, -2] = 1\nmixbarrank[C, 0, -4] = 0\nmixbarperiodicdim[C, 0] = 0\nmixbars[C] = free [], torsion [(0, 2)]\n"
        );
    }
}

#[test]
fn extremes_free_multiplicity_and_bad_bar_queries() {
    let program = r"
        mixedcomplex C[2] { a:9223372036854775806, b:-9223372036854775808,
                            c:-9223372036854775808, d:0 } {} {};
        mixbars[C];
        mixbardim[C,-9223372036854775808];
        mixbarrank[C,9223372036854775806,-9223372036854775808];
        mixbarperiodicdim[C,9223372036854775806];
        mixbarperiodicdim[C,-9223372036854775807];
        ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "mixbars[C] = free [-9223372036854775808, -9223372036854775808, 0, 9223372036854775806], torsion []\n",
            "mixbardim[C, -9223372036854775808] = 4\n",
            "mixbarrank[C, 9223372036854775806, -9223372036854775808] = 1\n",
            "mixbarperiodicdim[C, 9223372036854775806] = 4\n",
            "mixbarperiodicdim[C, -9223372036854775807] = 0\n",
        )
    );
    let example = include_str!("../examples/mixed-bars.asl");
    for (extra, expected) in [
        ("show mixBarNull[M,0];", "requires a torsion bar"),
        ("show mixBarRep[M,4];", "bar index is out of range"),
        ("mixbarrank[M,2,-1];", "nonnegative power of u"),
        ("mixbarrank[M,2,4];", "nonnegative power of u"),
    ] {
        let error = diagnostic(interpret(&format!("{example} {extra}")));
        assert!(error.contains(expected), "{extra}: {error}");
    }
}
