mod common;

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
fn discriminant_twenty_three_detects_conjugation_and_saturates() {
    assert_eq!(
        stdout(example("hyperelliptic-k3-shioda-inose")),
        concat!(
            "rvK3CMForm[P] = form=[4, 2, 6], content=2, discriminant=-92, order_discriminant=-23\n",
            "rvK3CMForm[Conjugate] = form=[4, -2, 6], content=2, discriminant=-92, order_discriminant=-23\n",
            "rvK3CMEquivalent[P,Conjugate] = false\n",
            "rvK3CMEquivalent[P,Shear] = true\n",
            "rvK3CMEquivalent[P,Nonprimitive] = true\n",
            "rvK3CMPeriods[P] = tau1=(-2 + sqrt(-92))/8, tau2=(2 + sqrt(-92))/2, branch=positive_imaginary\n",
            "rvK3CMClassNumber[P] = 3\n",
            "rvK3LatticePair = 8\n",
            "rvK3LatticePair = 2\n",
            "rvK3LatticePair = 12\n",
        ),
    );
}

#[test]
fn proper_boundary_reduction_and_orders() {
    let mut source = MARKED.to_string();
    for (name, a, b, c) in [
        ("A", 1, 1, 4),
        ("C", 2, 1, 2),
        ("D", 1, 0, 1),
        ("E", 1, 1, 1),
        ("Ff", 1, 0, 3),
    ] {
        let (u, v) = wedge_pair(a, b, c);
        source.push_str(&format!(
            "let u{name}={u}; let v{name}={v};
             rvk3period {name}(u{name},v{name});
             rvk3period reverse{name}(u{name},neg(v{name}));
             rvK3CMEquivalent[{name},reverse{name}];
             rvK3CMClassNumber[{name}];"
        ));
    }
    source.push_str("rvK3CMEquivalent[A,C]; rvK3CMEquivalent[D,E];");
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "rvK3CMEquivalent[A,reverseA] = true\n",
            "rvK3CMClassNumber[A] = 2\n",
            "rvK3CMEquivalent[C,reverseC] = true\n",
            "rvK3CMClassNumber[C] = 2\n",
            "rvK3CMEquivalent[D,reverseD] = true\n",
            "rvK3CMClassNumber[D] = 1\n",
            "rvK3CMEquivalent[E,reverseE] = true\n",
            "rvK3CMClassNumber[E] = 1\n",
            "rvK3CMEquivalent[Ff,reverseFf] = true\n",
            "rvK3CMClassNumber[Ff] = 1\n",
            "rvK3CMEquivalent[A,C] = false\n",
            "rvK3CMEquivalent[D,E] = false\n",
        ),
    );
}

#[test]
fn wide_unimodular_shear_and_swap_reduce_to_the_same_actual_basis() {
    let (u, v) = wedge_pair(2, 1, 3);
    let k = (1i128 << 124) - 1;
    let huge_v = wedge([2, 2 * k, 0, 0, 2 + 4 * k, -6]);
    let output = stdout(interpret(&format!(
        "{MARKED} let u={u}; let v={v}; let hv={huge_v};
         rvk3period P(u,v); rvk3period Large(u,hv); rvk3period Swapped(hv,neg(u));
         rvK3CMEquivalent[P,Large]; rvK3CMEquivalent[P,Swapped];
         assert rvK3CMBasis[Large,0]==u; assert rvK3CMBasis[Large,1]==v;
         rvK3CMForm[Swapped];"
    )));
    assert_eq!(
        output,
        concat!(
            "rvK3CMEquivalent[P,Large] = true\n",
            "rvK3CMEquivalent[P,Swapped] = true\n",
            "rvK3CMForm[Swapped] = form=[4, 2, 6], content=2, discriminant=-92, order_discriminant=-23\n",
        )
    );
}

#[test]
fn glued_basis_and_equivalence_across_markings() {
    let source = format!(
        "{MARKED}
         rvjactorsionbasis Twin[L,2](neg(oneW)){{e1,e2,f1,f2}};
         let u=rvK3Lattice[B]{{[0,2,0,0,2,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}};
         let v=rvK3Lattice[B]{{[3,0,0,0,0,-4],[0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1]}};
         let otherU=rvK3Lattice[Twin]{{[0,2,0,0,2,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}};
         let otherV=rvK3Lattice[Twin]{{[3,0,0,0,0,-4],[0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1]}};
         rvk3period P(u,v); rvk3period Q(otherU,otherV);
         assert rvK3CMBasis[P,0]==u; assert rvK3CMBasis[Q,1]==otherV;
         rvK3CMEquivalent[P,Q]; rvK3CMForm[P];"
    );
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "rvK3CMEquivalent[P,Q] = true\n",
            "rvK3CMForm[P] = form=[2, 0, 5], content=1, discriminant=-40, order_discriminant=-40\n",
        ),
    );
}

#[test]
fn saturation_keeps_orientation_with_an_index_beyond_u128() {
    let (u, v) = wedge_pair(2, 1, 3);
    let k = 1i128 << 122;
    let combination = |a, b| wedge([2 * b, 2 * a, 0, 0, 4 * a + 2 * b, -6 * b]);
    let first = combination(k, k + 1);
    let second = combination(-k - 2, k);
    // The determinant k²+(k+1)(k+2) is positive, and the first
    // coefficient column is primitive by gcd(k,k+1)=1.
    let output = stdout(interpret(&format!(
        "{MARKED} rvk3period P({u},{v}); rvk3period Wide({first},{second});
         rvK3CMEquivalent[P,Wide]; rvK3CMForm[Wide];"
    )));
    assert_eq!(
        output,
        concat!(
            "rvK3CMEquivalent[P,Wide] = true\n",
            "rvK3CMForm[Wide] = form=[4, 2, 6], content=2, discriminant=-92, order_discriminant=-23\n",
        )
    );
}

#[test]
fn large_forms_remain_available_when_class_number_enumeration_exceeds_budget() {
    let (u, v) = wedge_pair(2, 1, 1i128 << 125);
    let source = format!("{MARKED} rvk3period P({u},{v});");
    let output = stdout(interpret(&format!("{source} rvK3CMForm[P];")));
    assert!(
        output.starts_with("rvK3CMForm[P] = form=[4, 2, "),
        "{output}"
    );
    for (command, expected) in [
        (
            "rvK3CMClassNumber[P];",
            "K3 CM class number enumeration exceeds the computation limit",
        ),
        ("rvK3CMForm[Missing];", "unknown K3 period 'Missing'"),
        ("show rvK3CMBasis[P,2];", "basis index must be 0 or 1"),
    ] {
        let error = diagnostic(interpret(&format!("{source}{command}")));
        assert!(error.contains(expected), "{error}");
    }
}

fn wedge_pair(a: i128, b: i128, c: i128) -> (String, String) {
    (
        wedge([0, 2, 0, 0, 2 * a, 0]),
        wedge([2, 0, 0, 0, 2 * b, -2 * c]),
    )
}

fn wedge(w: [i128; 6]) -> String {
    format!("rvK3Lattice[B]{{{w:?},[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}}")
}
