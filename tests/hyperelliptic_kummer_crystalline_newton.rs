mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const EXAMPLE: &str = include_str!("../examples/hyperelliptic-crystalline-newton.asl");

#[test]
fn example_distinguishes_heights_and_the_odd_generalized_kummer_polygons() {
    assert_eq!(
        stdout(example("hyperelliptic-crystalline-newton")),
        concat!(
            "rvKummerK3Height[O] = 1\n",
            "rvKummerK3Height[P] = 2\n",
            "rvKummerK3Height[S] = infinite\n",
            "rvKummerHilbNewton[O, 2, 4] = [0: 1, 1: 21, 2: 232, 3: 21, 4: 1]\n",
            "rvKummerHilbNewton[P, 2, 4] = [1: 3, 3/2: 38, 2: 194, 5/2: 38, 3: 3]\n",
            "rvKummerHilbNewton[S, 2, 4] = [2: 276]\n",
            "rvJacHilbNewton[O, 2, 3] = [0: 2, 1: 14, 2: 14, 3: 2]\n",
            "rvGenKummerNewton[A, 3, 3] = [1: 4, 2: 4]\n",
            "rvGenKummerNewton[B, 3, 3] = [3/2: 8]\n",
        )
    );
}

#[test]
fn normalized_slopes_survive_base_extension_and_duality() {
    let out = stdout(interpret(&format!(
        r"{EXAMPLE}
        rvunram V[3,1]{{1,0,1}};
        rvjacbase P2[P,V];
        rvKummerK3Height[P2];
        rvKummerHilbNewton[P2,2,4];
        rvKummerHilbNewton[P,2,2];
        rvKummerHilbNewton[P,2,6];
        rvJacHilbNewton[O,2,5];
        rvGenKummerNewton[A,3,5];
        rvGenKummerNewton[O,2,2];
        rvKummerHilbNewton[O,1,2];
        "
    )));
    assert_eq!(
        out.lines().skip(9).collect::<Vec<_>>().join("\n"),
        concat!(
            "rvKummerK3Height[P2] = 2\n",
            "rvKummerHilbNewton[P2, 2, 4] = [1: 3, 3/2: 38, 2: 194, 5/2: 38, 3: 3]\n",
            "rvKummerHilbNewton[P, 2, 2] = [1/2: 2, 1: 19, 3/2: 2]\n",
            "rvKummerHilbNewton[P, 2, 6] = [5/2: 2, 3: 19, 7/2: 2]\n",
            "rvJacHilbNewton[O, 2, 5] = [1: 2, 2: 14, 3: 14, 4: 2]\n",
            "rvGenKummerNewton[A, 3, 5] = [2: 4, 3: 4]\n",
            "rvGenKummerNewton[O, 2, 2] = [0: 1, 1: 20, 2: 1]\n",
            "rvKummerHilbNewton[O, 1, 2] = [0: 1, 1: 20, 2: 1]",
        )
    );
}

#[test]
fn scalar_slopes_have_the_known_betti_multiplicities_in_every_degree() {
    let mut source = EXAMPLE.to_string();
    let jacobian = [1, 4, 13, 32, 44, 32, 13, 4, 1];
    let k3 = [1, 0, 23, 0, 276, 0, 23, 0, 1];
    let generalized = [1, 0, 7, 8, 108, 8, 7, 0, 1];
    let mut expected = Vec::new();
    for (command, curve, n, bettis) in [
        ("rvJacHilbNewton", "S", 2, jacobian),
        ("rvKummerHilbNewton", "S", 2, k3),
        ("rvGenKummerNewton", "B", 3, generalized),
    ] {
        for (i, betti) in bettis.into_iter().enumerate() {
            writeln!(source, "{command}[{curve},{n},{i}];").unwrap();
            let slope = if i % 2 == 0 {
                (i / 2).to_string()
            } else {
                format!("{i}/2")
            };
            let bins = if betti == 0 {
                String::new()
            } else {
                format!("{slope}: {betti}")
            };
            expected.push(format!("{command}[{curve}, {n}, {i}] = [{bins}]"));
        }
    }
    let out = stdout(interpret(&source));
    assert_eq!(out.lines().skip(9).collect::<Vec<_>>(), expected);
}

#[test]
fn composite_length_and_stable_low_degrees() {
    let out = stdout(interpret(&format!(
        r"{EXAMPLE}
        rvGenKummerNewton[O,4,6];
        rvGenKummerNewton[O,340282366920938463463374607431768211454,2];
        rvJacHilbNewton[O,340282366920938463463374607431768211455,1];
        rvJacHilbNewton[O,0,0];
        rvKummerHilbNewton[O,2,9];
        rvGenKummerNewton[O,1,0];
        rvGenKummerNewton[O,1,1];
        "
    )));
    assert_eq!(
        out.lines().skip(9).collect::<Vec<_>>(),
        [
            "rvGenKummerNewton[O, 4, 6] = [0: 1, 1: 5, 2: 37, 3: 372, 4: 37, 5: 5, 6: 1]",
            "rvGenKummerNewton[O, 340282366920938463463374607431768211454, 2] = [0: 1, 1: 5, 2: 1]",
            "rvJacHilbNewton[O, 340282366920938463463374607431768211455, 1] = [0: 2, 1: 2]",
            "rvJacHilbNewton[O, 0, 0] = [0: 1]",
            "rvKummerHilbNewton[O, 2, 9] = []",
            "rvGenKummerNewton[O, 1, 0] = [0: 1]",
            "rvGenKummerNewton[O, 1, 1] = []",
        ],
    );
}

#[test]
fn rejects_bad_characteristic_bad_models_and_large_coefficient_tables() {
    for (statement, expected) in [
        ("rvGenKummerNewton[O,0,0];", "length must be positive"),
        (
            "rvGenKummerNewton[O,3,0];",
            "length prime to the residue characteristic",
        ),
        (
            "rvGenKummerNewton[O,3,999];",
            "length prime to the residue characteristic",
        ),
        ("rvJacHilbNewton[O,30,60];", "coefficient storage limit"),
        (
            "rvKummerK3Height[Missing];",
            "unknown residue hyperelliptic curve",
        ),
        (
            "rvhyperelliptic E[F]{one,one,z,one}; rvJacHilbNewton[E,1,0];",
            "monic quintic",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{EXAMPLE}\n{statement}")));
        assert!(error.contains(expected), "{error}");
    }
}
