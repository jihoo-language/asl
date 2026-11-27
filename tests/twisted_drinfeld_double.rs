mod common;

use common::{diagnostic, example, interpret, stdout};

const CARRY: &str = "
    rvgroup G[3];rvdwcoeff P[G,3];let k=rvPostnikovCochain[P,3]{0,1,1,1,0,2,2,2};
    rvdouble D[P](k);
";

#[test]
fn cyclic_carry_flux_has_order_nine_and_the_lens_amplitude_matches_dw() {
    let actual = stdout(example("twisted-drinfeld-double"));
    let lines = actual.lines().collect::<Vec<_>>();
    assert_eq!(lines, [
        "rvDoubleData[D] = group_order=3, conductor=9, simples=9, dimension_counts={1: 9}, twist_counts={0: 3, 1: 2, 4: 2, 7: 2}, total_squared_dimension=9, modular_relations=true",
        "rvDoubleSector[D,3] = flux=1, dimension=1, radical_phases=[(0, 0), (1, 1), (2, 2)], dual=7, twist=1",
        "RVDouble(D){[(2, 1)]}",
        "RVCycloQ(9){[1, 0, 0, 2] / 3}",
    ]);
}

fn list<T: ToString>(entries: impl IntoIterator<Item = T>) -> String {
    entries
        .into_iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn table(order: usize, mul: impl Fn(usize, usize) -> usize) -> String {
    list((0..order).map(|i| format!("{{{}}}", list((0..order).map(|j| mul(i, j))))))
}

#[test]
fn type_three_cocycle_has_nonpointed_fusion_and_exact_modular_data() {
    let xor = table(8, |a, b| a ^ b);
    let k = list((1..8).flat_map(|a| {
        (1..8).flat_map(move |b| (1..8).map(move |c| (a & 1) * ((b >> 1) & 1) * ((c >> 2) & 1)))
    }));
    let mut source = format!(
        "rvgroup G{{{xor}}};rvdwcoeff P[G,2];
        let k=rvPostnikovCochain[P,3]{{{k}}};rvdouble D[P](k);"
    );
    source.push_str(
        r"
        rvDoubleData[D];
        rvDoubleSector[D,8];rvDoubleSector[D,20];rvDoubleSector[D,21];
        rvDoublePhase[D,8,2];
        let x=rvDoubleSimple[D,1]{0,0};
        let x2=rvDoubleSimple[D,1]{0,2};
        let y=rvDoubleSimple[D,2]{0,0};
        let z=rvDoubleSimple[D,4]{0,0};
        let xx=rvDoubleFuse(x,x);
        let answer=rvDoubleLinear[D]{1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        assert x==rvDoubleBasis[D,8];
        assert xx==answer;
        assert x==rvDoubleDual(x);
        assert rvDoubleFuse(rvDoubleFuse(x,y),z)==rvDoubleFuse(x,rvDoubleFuse(y,z));
        assert rvDoubleFuse(x,y)==rvDoubleFuse(y,x);
        assert rvDoubleFuse(add(x,neg(x2)),y)==add(rvDoubleFuse(x,y),neg(rvDoubleFuse(x2,y)));
        assert rvDoubleScale[2](x)==add(x,x);
        assert add(x,neg(x))==zero[RVDouble(D)];
        let one=rvCyclo[4,0];
        assert rvDoubleCharacter[D,8,0]==rvCycloRat[1](rvCycloScale[2](one));
        assert rvDoubleCharacter[D,8,2]==zero[RVCycloQ(4)];
        assert rvDoubleS[D,8,10]==zero[RVCycloQ(4)];
        assert rvDoubleS[D,0,8]==rvCycloRat[4](one);
        assert rvDoubleT[D,20]==rvCycloRat[1](rvCyclo[4,1]);
        assert rvDoubleT[D,21]==rvCycloRat[1](rvCyclo[4,3]);
        assert rvDoubleLens[D,0]==rvCycloRat[1](one);
        assert rvDoubleLens[D,1]==rvCycloRat[8](one);
        rvDoubleCoeff[3](xx);rvDoubleCoeff[4](xx);
        ",
    );
    let actual = stdout(interpret(&source));
    let lines = actual.lines().collect::<Vec<_>>();
    assert_eq!(lines[0], "rvDoubleData[D] = group_order=8, conductor=4, simples=22, dimension_counts={1: 8, 2: 14}, twist_counts={0: 14, 1: 1, 2: 6, 3: 1}, total_squared_dimension=64, modular_relations=true");
    assert!(lines[1].contains("flux=1, dimension=2, radical_phases=[(0, 0), (1, 0)]"));
    assert!(lines[2].contains("flux=7, dimension=2, radical_phases=[(0, 0), (7, 1)]"));
    assert!(lines[3].contains("flux=7, dimension=2, radical_phases=[(0, 0), (7, 3)]"));
    assert_eq!(lines[4], "rvDoublePhase[D,8,2] = outside radical");
    assert_eq!(lines[5..], ["rvDoubleCoeff[3] = 1", "rvDoubleCoeff[4] = 0"]);
}

#[test]
fn relabeled_cyclic_group_exercises_nonsplit_power_relations_in_the_radical_chain() {
    let relabel = [0, 2, 1, 3];
    let table = table(4, |a, b| {
        relabel
            .iter()
            .position(|&x| x == (relabel[a] + relabel[b]) % 4)
            .unwrap()
    });
    let source = format!(
        "rvgroup G{{{table}}};rvdwcoeff P[G,2];rvdouble U[P](rvPostnikovK[P]);
        rvDoubleData[U];
        let x=rvDoubleSimple[U,0]{{0,4,2,6}};
        let x2=rvDoubleSimple[U,0]{{0,0,4,4}};
        assert rvDoubleFuse(x,x)==x2;
        assert rvDoubleFuse(x2,x2)==rvDoubleBasis[U,0];
        assert rvDoubleFuse(x,rvDoubleDual(x))==rvDoubleBasis[U,0];
        assert rvDoubleLens[U,1]==rvCycloRat[4](rvCyclo[8,0]);"
    );
    let actual = stdout(interpret(&source));
    assert!(actual.contains("group_order=4, conductor=8, simples=16, dimension_counts={1: 16}"));
}

#[test]
fn gauge_transform_preserves_modular_data_with_rephased_characters() {
    let source = format!(
        "{CARRY}
        let shifted=rvPostnikovGauge(k,rvPostnikovCochain[P,2]{{0,1,0,0}});
        rvdouble Q[P](shifted);
        let x=rvDoubleSimple[D,1]{{0,1,2}};
        let xa=rvDoubleSimple[Q,1]{{0,1,8}};
        assert xa==rvDoubleBasis[Q,3];
        assert rvDoubleT[D,3]==rvDoubleT[Q,3];
        assert rvDoubleLens[D,3]==rvDoubleLens[Q,3];
        assert rvDoubleFuse(x,rvDoubleFuse(x,x))==rvDoubleSimple[D,0]{{0,6,3}};
        assert rvDoubleFuse(xa,rvDoubleFuse(xa,xa))==rvDoubleSimple[Q,0]{{0,6,3}};"
    );
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn scope_invalid_characters_and_nominal_operations_have_diagnostics() {
    for (suffix, expected) in [
        ("show rvDoubleSimple[D,1]{0,0,0};","no simple at this flux"),
        ("show rvDoubleBasis[D,9];","simple index is out of range"),
        ("show rvDoubleCharacter[D,3,4];","evaluation element is out of range"),
        ("rvdouble Q[P](k);show rvDoubleFuse(rvDoubleBasis[D,0],rvDoubleBasis[Q,0]);","expected type RVDouble(D)"),
        ("rvdouble Q[P](rvPostnikovCochain[P,3]{1,0,0,0,0,0,0,0});","fails the twisted cocycle"),
        ("show rvDoubleScale[2](rvDoubleScale[170141183460469231731687303715884105727](rvDoubleBasis[D,0]));","signed 128-bit coefficient range"),
        ("show rvDoubleLinear[D]{0,1};","requires 9 coefficients"),
    ] {
        let error = diagnostic(interpret(&format!("{CARRY}{suffix}")));
        assert!(error.contains(expected), "{expected} in {error}");
    }
    let s3 = table(6, |x, y| {
        let (a, b, c, d) = (x % 3, x / 3, y % 3, y / 3);
        (a + if b == 0 { c } else { 3 - c }) % 3 + 3 * ((b + d) % 2)
    });
    let source = format!("rvgroup G{{{s3}}};rvdwcoeff P[G,2];rvdouble D[P](rvPostnikovK[P]);");
    assert!(diagnostic(interpret(&source)).contains("requires an abelian group"));
    let twisted = "rvgroup H[3];rvgroup E[2];
        rvcrossed C[H,E]{0,0,0}{{0,1,2},{0,2,1}};
        rvpostnikov P[C,1];rvdouble D[P](rvPostnikovK[P]);";
    assert!(diagnostic(interpret(twisted)).contains("requires trivial coefficient action"));
}
