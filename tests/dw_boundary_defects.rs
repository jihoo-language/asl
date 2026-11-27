mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const EXAMPLE: &str = include_str!("../examples/twisted-double-boundary-defects.asl");
const KLEIN: &str = "
    rvgroup K{{0,1,2,3},{1,0,3,2},{2,3,0,1},{3,2,1,0}};
    rvdwcoeff O[K,2];
";

fn type_iii() -> &'static str {
    EXAMPLE.split_once("// The plane").unwrap().0
}

fn list<T: ToString>(values: impl IntoIterator<Item = T>) -> String {
    values
        .into_iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

#[test]
fn plane_and_klein_torsion_reproduce_the_independent_defect_counts() {
    let text = stdout(example("twisted-double-boundary-defects"));
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines[0], "rvDoubleDefectData[C] = boundaries=[Plane,Plane], intersection=[0, 1, 2, 3], cosets=2, conductor=4, simples=5, sum_squared_fiber_dimensions=8, annulus_verified=true");
    assert!(
        lines[1].ends_with("radical=[0, 1, 2, 3], first_simple=0, simples=4, fiber_dimension=1")
    );
    assert!(lines[2].ends_with("radical=[0], first_simple=4, simples=1, fiber_dimension=2"));
    assert!(lines[3].ends_with("isotropic=[0, 1], coset_representatives=[0, 2], fp_dimension=2"));
    assert_eq!(lines[4], "RVQZ(1/2)");
    assert_eq!(lines[5..9],[
        "rvDoubleDefectRep[1,1] = target=1, phase=2",
        "rvDoubleDefectRep[2,0] = target=1, phase=0",
        "rvDoubleDefectFuse[4,4] = coset=0, fiber_dimension=4, sectors=[(0, 1), (1, 1), (2, 1), (3, 1)]",
        "rvDoubleDefectFusionProjector[0,0,3] = RVCycloQ(4){[1] / 2}",
    ]);
    assert_eq!(lines[9], "rvDoubleDefectData[AB] = boundaries=[A,B], intersection=[0, 1, 2, 3], cosets=1, conductor=4, simples=1, sum_squared_fiber_dimensions=4, annulus_verified=true");
}

#[test]
fn type_iii_projectors_are_rank_one_orthogonal_and_resolve_the_balanced_fiber() {
    let mut source = format!(
        "{}rvdoubleboundary A[D]{{0,1,2,3}};rvdoubledefects C[A,A];
        let m=rvDoubleDefectBasis[C,4];
        let one=rvCycloRat[1](rvCyclo[4,0]);let minus=neg(one);
        assert rvDoubleDefectEntry[1,0,0](m)==one;
        assert rvDoubleDefectEntry[1,1,1](m)==minus;
        assert rvDoubleDefectEntry[2,1,0](m)==one;
        assert rvDoubleDefectEntry[2,0,1](m)==one;
        assert rvDoubleDefectEntry[3,1,0](m)==one;
        assert rvDoubleDefectEntry[3,0,1](m)==minus;
        assert rvDoubleDefectCharacter[1](m)==zero[RVCycloQ(4)];
        assert rvDoubleDefectCharacter[0](m)==add(one,one);",
        type_iii()
    );
    for c in 0..4 {
        for r in 0..4 {
            for s in 0..4 {
                writeln!(source, "rvDoubleDefectFusionProjector[{c},{r},{s}](m,m);").unwrap();
            }
        }
    }
    let data = stdout(interpret(&source));
    let coefficients = data
        .lines()
        .map(|line| match line.split_once(" = ").unwrap().1 {
            "RVCycloQ(4){[] / 1}" => 0,
            "RVCycloQ(4){[1] / 2}" => 1,
            "RVCycloQ(4){[-1] / 2}" => -1,
            other => panic!("unexpected projector coefficient {other}"),
        })
        .collect::<Vec<i32>>();
    assert_eq!(coefficients.len(), 64);
    let projectors = coefficients.chunks_exact(16).collect::<Vec<_>>();
    for (i, p) in projectors.iter().enumerate() {
        assert_eq!((0..4).map(|j| p[5 * j]).sum::<i32>(), 2);
        for (j, q) in projectors.iter().enumerate() {
            for r in 0..4 {
                for s in 0..4 {
                    let product = (0..4).map(|u| p[4 * r + u] * q[4 * u + s]).sum::<i32>();
                    assert_eq!(product, if i == j { 2 * p[4 * r + s] } else { 0 });
                }
            }
        }
    }
    for r in 0..4 {
        for s in 0..4 {
            assert_eq!(
                projectors.iter().map(|p| p[4 * r + s]).sum::<i32>(),
                2 * i32::from(r == s)
            );
        }
    }
}

#[test]
fn mixed_intersections_have_actual_distinct_double_cosets() {
    let source = format!(
        "{}rvdoubleboundary Plane[D]{{0,1,2,3}};rvdoubleboundary Line[D]{{0,1}};
         rvdoubleboundary Unit[D]{{0}};
         rvdoubledefects LP[Line,Plane];rvdoubledefects PL[Plane,Line];
         rvdoubledefects LL[Line,Line];rvdoubledefects EP[Unit,Plane];
         rvDoubleDefectData[LP];rvDoubleDefectData[LL];rvDoubleDefectData[EP];
         rvDoubleDefectCoset[LL,3];
         let f=\\x : RVDoubleDefect(LP) => rvDoubleDefectOpposite[LP](rvDoubleDefectOpposite[PL](x));
         assert f == (\\x : RVDoubleDefect(LP) => x);
         assert rvDoubleDefectCharacter[0](rvDoubleDefectBasis[EP,0])==rvCycloRat[1](rvCyclo[1,0]);",
        type_iii()
    );
    let lines = stdout(interpret(&source));
    let lines: Vec<_> = lines.lines().collect();
    assert!(lines[0].ends_with("intersection=[0, 1], cosets=2, conductor=2, simples=4, sum_squared_fiber_dimensions=4, annulus_verified=true"));
    assert!(lines[1].ends_with("intersection=[0, 1], cosets=4, conductor=2, simples=8, sum_squared_fiber_dimensions=8, annulus_verified=true"));
    assert!(lines[2].ends_with("intersection=[0], cosets=2, conductor=1, simples=2, sum_squared_fiber_dimensions=2, annulus_verified=true"));
    assert!(lines[3].contains("representative=6, support=[6, 7]"));
}

fn odd_mixed_carry() -> String {
    let mut group = vec![(0, 0), (1, 1), (0, 2)];
    for a in 0..3 {
        for b in 0..3 {
            if !group.contains(&(a, b)) {
                group.push((a, b));
            }
        }
    }
    let table = list(group.iter().map(|&(a, b)| {
        format!(
            "{{{}}}",
            list(group.iter().map(|&(c, d)| {
                group
                    .iter()
                    .position(|&e| e == ((a + c) % 3, (b + d) % 3))
                    .unwrap()
            }))
        )
    }));
    let k = list((1..9).flat_map(|a| {
        (1..9).flat_map({
            let group = &group;
            move |b| (1..9).map(move |c| group[a].0 * ((group[b].1 + group[c].1) / 3) % 3)
        })
    }));
    let subgroup = list(
        group
            .iter()
            .enumerate()
            .filter_map(|(i, &(_, b))| (b == 0).then_some(i)),
    );
    format!(
        "rvgroup G{{{table}}};rvdwcoeff P[G,3];let k=rvPostnikovCochain[P,3]{{{k}}};rvdouble D[P](k);
         rvdoubleboundary H[D]{{{subgroup}}};rvdoubledefects C[H,H];"
    )
}

#[test]
fn opposite_rephases_an_odd_complex_simple_at_a_noncanonical_inverse() {
    let source = format!(
        "{}
        let x=rvDoubleDefectBasis[C,3];
        let y=rvDoubleDefectBasis[C,6];
        let xd=rvDoubleDefectOpposite[C](x);
        assert xd==rvDoubleDefectBasis[C,8];
        assert rvDoubleDefectOpposite[C](xd)==x;
        assert rvDoubleDefectCharacter[4](xd)==rvCycloRat[1](rvCyclo[3,2]);
        rvDoubleDefectFuse(x,x);rvDoubleDefectFuse(x,y);rvDoubleDefectFuse(x,xd);",
        odd_mixed_carry()
    );
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "rvDoubleDefectFuse[3,3] = coset=2, fiber_dimension=1, sectors=[(6, 1)]\n",
            "rvDoubleDefectFuse[3,6] = coset=0, fiber_dimension=1, sectors=[(1, 1)]\n",
            "rvDoubleDefectFuse[3,8] = coset=0, fiber_dimension=1, sectors=[(0, 1)]\n"
        )
    );
}

#[test]
fn opposite_grade_translation_orders_the_left_action_before_the_fiber() {
    let mut group = vec![[0, 0, 0], [0, 1, 1], [0, 0, 2]];
    for a in 0..3 {
        for b in 0..3 {
            for c in 0..3 {
                if !group.contains(&[a, b, c]) {
                    group.push([a, b, c]);
                }
            }
        }
    }
    let table = list(group.iter().map(|a| {
        format!(
            "{{{}}}",
            list(group.iter().map(|b| {
                let product = [0, 1, 2].map(|i| (a[i] + b[i]) % 3);
                group.iter().position(|g| *g == product).unwrap()
            }))
        )
    }));
    let k = list((1..27).flat_map(|a| {
        (1..27).flat_map({
            let group = &group;
            move |b| (1..27).map(move |c| group[a][0] * group[b][1] * group[c][2] % 3)
        })
    }));
    let line = list(
        group
            .iter()
            .enumerate()
            .filter_map(|(i, &[_, b, c])| (b == 0 && c == 0).then_some(i)),
    );
    let plane = list(
        group
            .iter()
            .enumerate()
            .filter_map(|(i, &[_, _, c])| (c == 0).then_some(i)),
    );
    let source = format!(
        "rvgroup G{{{table}}};rvdwcoeff P[G,3];let k=rvPostnikovCochain[P,3]{{{k}}};rvdouble D[P](k);
         rvdoubleboundary H[D]{{{line}}};rvdoubleboundary K[D]{{{plane}}};
         rvdoubledefects C[H,K];rvdoubledefects E[K,H];
         let x=rvDoubleDefectBasis[C,3];let dual=rvDoubleDefectOpposite[E](x);
         assert dual==rvDoubleDefectBasis[E,6];
         assert rvDoubleDefectOpposite[C](dual)==x;"
    );
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn quarter_bockstein_boundaries_produce_conductor_eight() {
    let lift = |a: i32, b: i32| (a & 1) * ((b >> 1) & 1);
    let k = list((1..4).flat_map(|a| {
        (1..4).flat_map(move |b| {
            (1..4).map(move |c| {
                ((lift(b, c) - lift(a ^ b, c) + lift(a, b ^ c) - lift(a, b)) / 2).rem_euclid(2)
            })
        })
    }));
    let supplied = list((1..4).flat_map(|a| (1..4).map(move |b| format!("({},4)", lift(a, b)))));
    let source = format!(
        "{KLEIN}let k=rvPostnikovCochain[O,3]{{{k}}};rvdouble D[O](k);
         rvdoubleboundary A[D]{{0,1,2,3}};
         rvdoubleboundary B[D]{{0,1,2,3}}{{{supplied}}};
         rvdoubledefects C[A,B];rvdoubledefects E[B,A];rvDoubleDefectData[C];
         assert rvDoubleDefectTheta[C,0,1,1]==rvQZ[3,4];
         assert rvDoubleDefectTheta[E,0,1,1]==rvQZ[1,4];
         let m=rvDoubleDefectBasis[C,0];
         assert rvDoubleDefectEntry[1,0,0](m)==rvCycloRat[1](rvCyclo[8,3]);
         assert rvDoubleDefectOpposite[E](m)==rvDoubleDefectBasis[E,0];"
    );
    assert_eq!(
        stdout(interpret(&source)),
        "rvDoubleDefectData[C] = boundaries=[A,B], intersection=[0, 1, 2, 3], cosets=1, conductor=8, simples=1, sum_squared_fiber_dimensions=4, annulus_verified=true\n"
    );
}

#[test]
fn defects_reject_invalid_boundaries_domains_and_self_fusion() {
    let data = format!(
        "{KLEIN}rvdouble D[O](rvPostnikovK[O]);rvdoubleboundary A[D]{{0,1,2,3}};
         rvdoubleboundary B[D]{{0,1}};
         rvdoubledefects C[A,B];rvdoubledefects E[B,A];"
    );
    for (operation, expected) in [
        (
            "show rvDoubleDefectBasis[C,4];",
            "simple index is out of range",
        ),
        ("rvDoubleDefectCoset[C,2];", "coset index is out of range"),
        (
            "show rvDoubleDefectEntry[2,0,0](rvDoubleDefectBasis[C,0]);",
            "boundary intersection",
        ),
        (
            "show rvDoubleDefectOpposite[C](rvDoubleDefectBasis[C,0]);",
            "reverse the same two nominal",
        ),
        (
            "rvDoubleDefectFuse(rvDoubleDefectBasis[C,0],rvDoubleDefectBasis[C,0]);",
            "requires one nominal boundary",
        ),
        (
            "rvDoubleDefectFuse(rvDoubleDefectBasis[C,0],rvDoubleDefectBasis[E,0]);",
            "expected type RVDoubleDefect(C)",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{data}{operation}")));
        assert!(error.contains(expected), "{expected} in {error}");
    }
    let source = format!(
        "{}rvdoubleboundary A[D]{{0}};rvdoubleboundary O[D]{{0,1,2,3,4,5,6,7}};rvdoubledefects F[A,O];",
        type_iii()
    );
    assert!(diagnostic(interpret(&source)).contains("has an obstructed subgroup"));
}
