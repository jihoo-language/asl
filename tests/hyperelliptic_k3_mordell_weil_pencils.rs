mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const BASIS: &str = r"
    let r0=rvJacPoint[C](z,z);
    let r1=rvJacPoint[C](o,z);
    let r2=rvJacPoint[C](two,z);
    let r3=rvJacPoint[C](three,z);
    rvjactorsionbasis B[C,2](neg(o)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),r3)};
";

const F5: &str = r"
    rvunram F[5,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let two=rvUnramElem[F]{2};
    let three=rvUnramElem[F]{3};
    rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
";

const F7: &str = r"
    rvunram F[7,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let two=rvUnramElem[F]{2};
    let three=rvUnramElem[F]{3};
    let four=rvUnramElem[F]{4};
    rvhyperelliptic C[F]{z,three,neg(o),z,four,o};
";

#[test]
fn all_zero_sections_and_trope_coordinates_match_the_branch_incidence() {
    let mut program = format!("{F5}{BASIS} rvk3pencil P[B,0,1];");
    let section_masks = [4u8, 5, 6, 7, 12, 13, 14, 15];
    for zero in section_masks {
        writeln!(
            program,
            "rvk3mordell M{zero}[P,{zero}]; rvK3MordellLattice[M{zero}];"
        )
        .unwrap();
        for trope in 0..16 {
            writeln!(
                program,
                "rvK3MordellCoord[M{zero}](rvKummerTropeDiv[B,{trope}]);"
            )
            .unwrap();
        }
    }
    let output = stdout(interpret(&program));
    let mut lines = output.lines();
    let branch = [0, 1, 4, 7, 13, 15];
    for zero in section_masks {
        let report = lines.next().unwrap();
        assert!(report.contains(
            "D4 stars=[(0, [4, 7, 13, 15]), (1, [5, 6, 12, 14])], A1 nodes=[2, 3, 8, 9, 10, 11]"
        ));
        assert!(report.contains("trivial rank=16, trivial |det|=1024, marked rank=17, marked |det|=64, quotient=Z+(Z/2)^2"));
        assert!(report.ends_with("shift free=-2, shift height=4"));
        let on_left = branch.contains(&zero);
        let (torsion, free) = if on_left {
            ("[4, 7, 13, 15]", "[(5, -1), (6, -1), (12, -1), (14, -1)]")
        } else {
            ("[5, 6, 12, 14]", "[(4, 1), (7, 1), (13, 1), (15, 1)]")
        };
        assert!(report.contains(&format!("torsion sections={torsion}, free sections={free}")));
        for trope in 0..16 {
            let left = i32::from(branch.contains(&trope));
            let right = i32::from(branch.contains(&(trope ^ 1)));
            let degree = 2 - left - right;
            let origin_sign = if on_left { 1 } else { -1 };
            let coordinate = (left - right - degree * origin_sign) / 2;
            assert_eq!(
                lines.next().unwrap(),
                format!("rvK3MordellCoord[M{zero}] = degree={degree}, free={coordinate}")
            );
        }
    }
    assert!(lines.next().is_none());
}

#[test]
fn divisor_translation_matches_both_orders_of_node_projections_on_generators() {
    let mut program = format!(
        r"{F5}{BASIS}
        rvk3pencil P[B,1,0];
        rvk3mordell M[P,4];
        rvk3autom Forward[B]{{project(0),project(1)}};
        rvk3autom Back[B]{{project(1),project(0)}};
        let oClass=rvKummerTropeDiv[B,4];
        let gClass=rvKummerTropeDiv[B,5];
        let sClass=add(rvKummerNodeDiv[B,0],neg(rvKummerNodeDiv[B,1]));
        let fClass=rvK3PencilDiv[P];
        rvK3MordellCoord[M](fClass);
        rvK3MordellCoord[M](oClass);
        rvK3MordellCoord[M](rvK3MordellShiftDiv[M,1](oClass));
        rvK3MordellCoord[M](rvK3MordellShiftDiv[M,-3](oClass));
        rvK3MordellHeight[M](sClass,gClass);
        rvK3MordellHeight[M](sClass,sClass);
        rvK3MordellHeight[M](sClass,fClass);
        rvK3MordellHeight[M](oClass,gClass);
        "
    );
    let generators = ["rvKummerHyperplane[B]".to_string()]
        .into_iter()
        .chain((0..16).flat_map(|mask| {
            [
                format!("rvKummerNodeDiv[B,{mask}]"),
                format!("rvKummerTropeDiv[B,{mask}]"),
            ]
        }));
    for (index, generator) in generators.enumerate() {
        writeln!(
            program,
            "let d{index}={generator};
            assert rvK3MordellShiftDiv[M,1](d{index})==rvK3AutoDiv[Forward](d{index});
            assert rvK3MordellShiftDiv[M,-1](d{index})==rvK3AutoDiv[Back](d{index});
            assert rvK3MordellShiftDiv[M,2](d{index})==rvK3AutoDiv[Forward](rvK3AutoDiv[Forward](d{index}));
            assert rvK3MordellShiftDiv[M,-7](rvK3MordellShiftDiv[M,7](d{index}))==d{index};"
        )
        .unwrap();
    }
    program.push_str(
        "rvKummerDivPair(oClass,gClass);
         rvKummerDivPair(rvK3MordellShiftDiv[M,9](oClass),rvK3MordellShiftDiv[M,9](gClass));",
    );
    let output = stdout(interpret(&program));
    let mut lines = output.lines();
    for expected in [
        "rvK3MordellCoord[M] = degree=0, free=0",
        "rvK3MordellCoord[M] = degree=1, free=0",
        "rvK3MordellCoord[M] = degree=1, free=-2",
        "rvK3MordellCoord[M] = degree=1, free=6",
        "rvK3MordellHeight[M] = 2",
        "rvK3MordellHeight[M] = 4",
        "rvK3MordellHeight[M] = 0",
        "rvK3MordellHeight[M] = 0",
        "rvKummerDivPair = 0",
        "rvKummerDivPair = 0",
    ] {
        assert_eq!(lines.next().unwrap(), expected);
    }
    assert!(lines.next().is_none());
}

#[test]
fn actual_shifts_preserve_every_fiber_including_exceptional_and_twist_points() {
    for (field, i, j, zero) in [(F5, 0, 1, 4), (F7, 4, 7, 4)] {
        let i_points = [
            ("zero[RVKummer(C)]", "rvJacKummerPoint(r0)"),
            ("rvJacKummerPoint(r1)", "rvJacKummerPoint(r2)"),
        ];
        let (a, b) = i_points[usize::from(i != 0)];
        let mut program = format!(
            r"{field}{BASIS}
            rvk3pencil P[B,{j},{i}];
            rvk3mordell M[P,{zero}];
            rvk3autom W[B]{{project({i}),project({j})}};
            let iNode={a};
            let jNode={b};
            let expected=\x:RVKummerK3(C)=>rvKummerNodeInvolution(jNode,rvKummerNodeInvolution(iNode,x));
            let reversed=\x:RVKummerK3(C)=>rvKummerNodeInvolution(iNode,rvKummerNodeInvolution(jNode,x));
            let pencil=\x:RVKummerK3(C)=>rvK3PencilMap[P](x);
            assert (\x:RVKummerK3(C)=>rvK3MordellShift[M,1](x))==expected;
            assert (\x:RVKummerK3(C)=>rvK3MordellShift[M,-1](x))==reversed;
            assert (\x:RVKummerK3(C)=>rvK3MordellShift[M,-1](expected(x)))==(\x:RVKummerK3(C)=>x);
            assert (\x:RVKummerK3(C)=>rvK3MordellShift[M,0](x))==(\x:RVKummerK3(C)=>x);
            assert (\x:RVKummerK3(C)=>rvK3MordellShift[M,2](x))==(\x:RVKummerK3(C)=>expected(expected(x)));
            assert (\x:RVKummerK3(C)=>rvK3PencilMap[P](expected(x)))==pencil;
            assert (\x:RVKummerK3(C)=>rvK3PencilMap[P](reversed(x)))==pencil;
            let p=rvK3PencilParam[P](o,two);
            let x=rvK3PencilSection[P,{zero}](p);
            assert rvK3PencilMap[P](rvK3MordellShift[M,-170141183460469231731687303715884105728](x))==p;
            assert rvK3MordellShift[M,-170141183460469231731687303715884105728](x)==rvK3AutoIter[W,-170141183460469231731687303715884105728](x);
            rvK3MordellLattice[M];
            rvK3AutoCycles[W];
            "
        );
        if i != 0 {
            program.push_str(
                "let twist=rvKummer[C](z,o,neg(o),o);
                 assert rvKummerLiftable(twist)==S.lo;
                 let tx=rvKummerSmooth(twist);
                 assert rvK3PencilMap[P](rvK3MordellShift[M,1](tx))==rvK3PencilMap[P](tx);",
            );
        }
        let output = stdout(interpret(&program));
        assert!(output.contains(&format!("nodes=[{i}, {j}], zero={zero}")));
        let cycles = output.lines().last().unwrap();
        assert!(cycles.contains("cycles=[("));
    }
}

#[test]
fn wide_translation_cancels_before_narrowing_and_checks_the_final_range() {
    let k = 1i128 << 124;
    let big = i128::MAX - 3;
    let left = big + 2;
    let right = big - 2;
    let min = i128::MIN;
    let program = format!(
        r"{F5}{BASIS}
        rvk3pencil P[B,0,1];
        rvk3mordell M[P,4];
        let x=rvKummerDiv[B]{{-{big},{left},{right},0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
        let y=rvKummerDiv[B]{{4,-2,-6,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
        let f=rvK3PencilDiv[P];
        assert rvK3MordellShiftDiv[M,{k}](x)==y;
        assert rvK3MordellShiftDiv[M,-{k}](y)==x;
        assert rvK3MordellShiftDiv[M,{min}](f)==f;
        rvK3MordellHeight[M](rvKummerDiv[B]{{{big},2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}},f);
        "
    );
    assert_eq!(stdout(interpret(&program)), "rvK3MordellHeight[M] = 0\n");
    let error = diagnostic(interpret(&format!(
        "{program} show rvK3MordellShiftDiv[M,{min}](rvKummerTropeDiv[B,4]);"
    )));
    assert!(
        error.contains("result exceeds signed 128-bit range"),
        "{error}"
    );
}

#[test]
fn wrong_zero_tropes_curves_and_markings_are_rejected() {
    let setup = format!("{F5}{BASIS} rvk3pencil P[B,0,1]; rvk3mordell M[P,4];");
    for mask in [0, 1, 2, 3, 8, 9, 10, 11] {
        let error = diagnostic(interpret(&format!("{setup} rvk3mordell Bad[P,{mask}];")));
        assert!(error.contains("must contain exactly one"), "{error}");
    }
    for (operation, expected) in [
        ("rvk3mordell Bad[Unknown,4];", "unknown K3 pencil"),
        ("rvK3MordellLattice[Unknown];", "unknown K3 Mordell–Weil marking"),
        ("rvk3mordell M[P,4];", "already defined"),
        ("rvk3mordell N[P,16];", "zero trope mask must be between 0 and 15"),
        ("rvK3MordellHeight[M](zero[RVKummer(C)],rvK3PencilDiv[P]);", "expected type RVKummerDiv(B)"),
        ("rvjactorsionbasis B2[C,2](neg(o)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),r3)}; rvK3MordellCoord[M](rvKummerHyperplane[B2]);", "expected type RVKummerDiv(B)"),
        ("show rvK3MordellShift[M,1](rvK3PencilDiv[P]);", "expected type RVKummerK3(C)"),
        ("rvhyperelliptic Twin[F]{z,neg(o),z,z,z,o}; show rvK3MordellShift[M,1](rvKummerSmooth(rvKummer[Twin](z,o,two,rvUnramMul(two,two))));", "expected type RVKummerK3(C), found RVKummerK3(Twin)"),
    ] {
        let error = diagnostic(interpret(&format!("{setup} {operation}")));
        assert!(error.contains(expected), "{error}");
    }
}

#[test]
fn runnable_example_reports_the_free_shift_and_its_height() {
    let output = stdout(example("hyperelliptic-k3-mordell-weil-pencil"));
    assert!(output.contains("shift free=-2, shift height=4"));
    assert!(output.contains("rvK3MordellCoord[M] = degree=1, free=-2\n"));
    assert!(output.contains("rvK3MordellHeight[M] = 4\n"));
    assert!(output.ends_with("rvK3AutoPeriod[ProjectPair] = 3\n"));
}
