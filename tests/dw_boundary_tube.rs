mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const EXAMPLE: &str = include_str!("../examples/twisted-double-boundary-tube.asl");

fn type_iii() -> &'static str {
    EXAMPLE.split_once("// The plane boundary").unwrap().0
}

#[test]
fn example_builds_both_type_iii_tubes_and_a_complex_carry_tube() {
    let output = stdout(example("twisted-double-boundary-tube"));
    assert!(output.contains(
        "dimension=40, center=22, boundary_simples=5, conductor=4, blocks={1: 16, 2: 6}"
    ));
    assert!(output.contains(
        "dimension=64, center=22, boundary_simples=8, conductor=4, blocks={1: 8, 2: 14}"
    ));
    assert!(output.contains("rvDoubleTubeAnnular[U,4,4] = dimension=8"));
    assert!(output.contains("rvDoubleTubeBlockData[T,8] = dimension=2, channels=[[0, 0], [1, 0]]"));
    assert!(output.contains("rvDoubleTubeDecompose = conductor=4, blocks=[(8,0,0,[1]/1)]"));
    assert!(output.contains("rvDoubleTubeDecompose = conductor=4, blocks=[(8,1,1,[1]/1)]"));
    assert!(output.contains("rvDoubleTubeStack[36,37] = terms=4, verified=true"));
    assert!(
        output.contains("dimension=9, center=9, boundary_simples=3, conductor=9, blocks={1: 9}")
    );
    assert!(output.ends_with("RVCycloQ(9){[0, 0, 0, 1] / 1}\n"));
}

#[test]
fn all_type_iii_annular_dimensions_match_independent_fusion_counts() {
    let mut source = format!("{}rvdoubletube T[PP];rvdoubletube U[QQ];", type_iii());
    let expected = [
        [5, 1, 1, 1, 0],
        [1, 5, 1, 1, 0],
        [1, 1, 5, 1, 0],
        [1, 1, 1, 5, 0],
        [0, 0, 0, 0, 8],
    ];
    let mut answer = String::new();
    for (j, row) in expected.into_iter().enumerate() {
        for (i, n) in row.into_iter().enumerate() {
            writeln!(source, "rvDoubleTubeAnnular[T,{j},{i}];").unwrap();
            writeln!(answer, "rvDoubleTubeAnnular[T,{j},{i}] = dimension={n}").unwrap();
        }
    }
    for j in 0..8 {
        for i in 0..8 {
            writeln!(source, "rvDoubleTubeAnnular[U,{j},{i}];").unwrap();
            let n = if i == j { 8 } else { 0 };
            writeln!(answer, "rvDoubleTubeAnnular[U,{j},{i}] = dimension={n}").unwrap();
        }
    }
    assert_eq!(stdout(interpret(&source)), answer);
}

#[test]
fn primitive_central_idempotents_sum_to_the_geometric_unit_and_split_matrix_units() {
    let mut source = format!(
        "{}rvdoubletube T[PP];let zeroT=rvDoubleTubeZero[T];let unit=rvDoubleTubeUnit[T];",
        type_iii()
    );
    for a in 0..22 {
        writeln!(
            source,
            "let p{a}=rvDoubleTubeCentral[T,{a}];
            assert rvDoubleTubeProduct(p{a},p{a})==p{a};
            assert rvDoubleTubeProduct(unit,p{a})==p{a};"
        )
        .unwrap();
    }
    for a in 0..22 {
        for b in 0..a {
            writeln!(source, "assert rvDoubleTubeProduct(p{a},p{b})==zeroT;").unwrap();
        }
    }
    let sum = (0..22)
        .map(|a| format!("p{a}"))
        .reduce(|a, b| format!("rvDoubleTubeAdd({a},{b})"))
        .unwrap();
    writeln!(source, "assert {sum}==unit;").unwrap();
    source.push_str(
        "
        let x=rvDoubleTubeMatrixUnit[T,8,0,1];let y=rvDoubleTubeMatrixUnit[T,8,1,0];
        let half=rvCycloRat[2](rvCyclo[4,0]);let quarter=rvCycloRat[4](rvCyclo[4,0]);
        assert x==rvDoubleTubeScale(half,rvDoubleTubeBasis[T,5]);
        assert y==rvDoubleTubeScale(half,rvDoubleTubeBasis[T,8]);
        assert rvDoubleTubeCoefficient[0](rvDoubleTubeProduct(x,y))==quarter;
        assert rvDoubleTubeCoefficient[2](rvDoubleTubeProduct(x,y))==neg(quarter);
        assert rvDoubleTubeCoefficient[9](rvDoubleTubeProduct(y,x))==quarter;
        assert rvDoubleTubeProduct(x,y)==rvDoubleTubeMatrixUnit[T,8,0,0];
        assert rvDoubleTubeProduct(y,x)==rvDoubleTubeMatrixUnit[T,8,1,1];
        assert rvDoubleTubeProduct(x,x)==zeroT;assert rvDoubleTubeProduct(x,p8)==x;
        assert rvDoubleTubeBlockEntry[8,0,1](x)==rvCycloRat[1](rvCyclo[4,0]);
        assert rvDoubleTubeBlockEntry[8,1,0](x)==zero[RVCycloQ(4)];
        assert rvDoubleTubeAdd(rvDoubleTubeMatrixUnit[T,8,0,0],rvDoubleTubeMatrixUnit[T,8,1,1])==p8;
        rvDoubleTubeDecompose(x);",
    );
    assert_eq!(
        stdout(interpret(&source)),
        "rvDoubleTubeDecompose = conductor=4, blocks=[(8,0,1,[1]/1)]\n"
    );
}

#[test]
fn annular_stacking_checks_noninvertible_loop_resolutions_and_pauli_signs() {
    let source = format!(
        "{}rvdoubletube T[PP];rvdoubletube U[QQ];
        rvDoubleTubeStack[T,36,37];rvDoubleTubeStack[T,39,39];
        rvDoubleTubeStack[T,5,8];rvDoubleTubeStack[T,8,5];
        rvDoubleTubeStack[T,5,5];
        rvDoubleTubeStack[U,33,34];rvDoubleTubeStack[U,34,33];
        let f=rvDoubleTubeBasis[U,33];let g=rvDoubleTubeBasis[U,34];
        let h=rvDoubleTubeBasis[U,35];
        let minus=neg(rvCycloRat[1](rvCyclo[4,0]));
        assert rvDoubleTubeProduct(f,g)==rvDoubleTubeScale(minus,h);
        assert rvDoubleTubeProduct(g,f)==h;",
        type_iii()
    );
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "rvDoubleTubeStack[36,37] = terms=4, verified=true\n",
            "rvDoubleTubeStack[39,39] = terms=4, verified=true\n",
            "rvDoubleTubeStack[5,8] = terms=4, verified=true\n",
            "rvDoubleTubeStack[8,5] = terms=4, verified=true\n",
            "rvDoubleTubeStack[5,5] = terms=0, verified=true\n",
            "rvDoubleTubeStack[33,34] = terms=1, verified=true\n",
            "rvDoubleTubeStack[34,33] = terms=1, verified=true\n",
        )
    );
}

#[test]
fn all_cyclic_three_stacks_give_the_dpr_multiplier_and_fourier_inverses() {
    let mut source =
        "rvgroup G[3];rvdwcoeff P[G,3];rvdouble D[P](rvPostnikovCochain[P,3]{0,1,1,1,0,2,2,2});
        rvdoubleboundary H[D]{0};rvdoubledefects C[H,H];rvdoubletube T[C];"
            .to_owned();
    for g in 0..3 {
        for x in 0..3 {
            let a = 3 * g + x;
            for y in 0..3 {
                let b = 3 * g + y;
                let c = 3 * g + (x + y) % 3;
                let phase = 3 * g * ((x + y) / 3);
                writeln!(source,"rvDoubleTubeStack[T,{a},{b}];
                    assert rvDoubleTubeProduct(rvDoubleTubeBasis[T,{a}],rvDoubleTubeBasis[T,{b}])
                    ==rvDoubleTubeScale(rvCycloRat[1](rvCyclo[9,{phase}]),rvDoubleTubeBasis[T,{c}]);").unwrap();
            }
        }
    }
    source.push_str(
        "assert rvDoubleTubeBlockEntry[3,0,0](rvDoubleTubeBasis[T,4])==rvCycloRat[1](rvCyclo[9,1]);
        let e=rvDoubleTubeCentral[T,3];
        assert rvDoubleTubeCoefficient[3](e)==rvCycloRat[3](rvCyclo[9,0]);
        assert rvDoubleTubeCoefficient[4](e)==rvCycloRat[3](rvCyclo[9,8]);
        rvDoubleTubeDecompose(e);",
    );
    let output = stdout(interpret(&source));
    assert_eq!(output.matches("terms=1, verified=true").count(), 27);
    assert!(output.ends_with("rvDoubleTubeDecompose = conductor=9, blocks=[(3,0,0,[1]/1)]\n"));
}

#[test]
fn tube_declarations_reject_mixed_nominal_operands_and_bad_indices() {
    let setup = format!("{}rvdoubletube T[PP];rvdoubletube U[QQ];", type_iii());
    for (op, want) in [
        ("rvdoubletube X[Y];", "nominal self-boundary"),
        (
            "show rvDoubleTubeBasis[T,40];",
            "basis index is out of range",
        ),
        (
            "show rvDoubleTubeMatrixUnit[T,8,2,1];",
            "outside its bulk block",
        ),
        (
            "show rvDoubleTubeProduct(rvDoubleTubeUnit[T],rvDoubleTubeUnit[U]);",
            "expected type RVDoubleTube(T)",
        ),
        (
            "show rvDoubleTubeScale(rvCycloRat[1](rvCyclo[9,0]),rvDoubleTubeUnit[T]);",
            "expected type RVCycloQ(4)",
        ),
        ("rvDoubleTubeStack[T,40,0];", "basis index is out of range"),
    ] {
        let error = diagnostic(interpret(&format!("{setup}{op}")));
        assert!(error.contains(want), "{want}: {error}");
    }
}

#[test]
fn a_boundary_coboundary_keeps_its_cochain_conductor() {
    let b = |i| i32::from(i == 1);
    let mu = (1..4)
        .flat_map(|x| (1..4).map(move |y| format!("({},3)", b(x) + b(y) - b(x ^ y))))
        .collect::<Vec<_>>()
        .join(",");
    let source = format!(
        "rvgroup G{{{{0,1,2,3}},{{1,0,3,2}},{{2,3,0,1}},{{3,2,1,0}}}};
         rvdwcoeff P[G,2];rvdouble D[P](rvPostnikovK[P]);
         rvdoubleboundary H[D]{{0,1,2,3}}{{{mu}}};
         rvdoubledefects C[H,H];rvdoubletube T[C];
         rvDoubleTubeData[T];rvDoubleTubeStack[T,5,6];"
    );
    assert_eq!(stdout(interpret(&source)),concat!(
        "rvDoubleTubeData[T] = dimension=16, center=16, boundary_simples=4, conductor=12, blocks={1: 16}, fourier_verified=true\n",
        "rvDoubleTubeStack[5,6] = terms=1, verified=true\n",
    ));
}

#[test]
fn an_oversized_exact_tube_has_a_specific_bound() {
    let source = "rvgroup G[12];rvdwcoeff P[G,2];rvdouble D[P](rvPostnikovK[P]);
        rvdoubleboundary Z[D]{0};rvdoubledefects C[Z,Z];rvdoubletube T[C];";
    let error = diagnostic(interpret(source));
    assert!(
        error.contains("supports at most 128 annular basis elements"),
        "{error}"
    );
}
