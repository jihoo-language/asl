mod common;

use common::{diagnostic, example, interpret, stdout};

const EXAMPLE: &str = "dwork-fricke-cyclic-isogeny-towers";

#[test]
fn primitive_level_nine_has_the_reference_coefficients_and_resultant_identity() {
    let out = stdout(example(EXAMPLE));
    assert!(out.contains("rvDworkHeckeEquation[9] = bidegree=(12, 12), terms=139"));
    assert!(out.contains("certified Res(Phi_3,Phi_3)=(X-Y)^4 Phi_9 over Z"));
    assert!(out.contains("rvDworkHeckeCoeff[9, 0, 0] = 4571271210813549960910454993738206238929240990423040610826738238226432"));
}

#[test]
fn local_towers_preserve_special_and_nonsplit_multiplicities() {
    let out = stdout(example(EXAMPLE));
    let p41 = out
        .lines()
        .find(|row| row.starts_with("rvDworkHeckeTower[3, 3] = H=[0]"))
        .unwrap();
    assert!(p41.contains("resultant degree=16, cyclic degree=12, backtracking divisor degree=4, cyclic diagonal multiplicity=0"));
    assert_eq!(p41.matches("(2, 1,").count(), 2);
    assert!(p41.contains("(6, 1,"));
    assert!(p41.ends_with("rational roots (value, multiplicity)=[([36], 2)]"));
    assert!(out.contains("rvDworkHeckeTower[3, 5] = H=[0], resultant degree=24, cyclic degree=24, coprime composition agrees in both orders"));

    let extra = stdout(interpret(
        "rvunram K[23,1]{0,1}; rvDworkHeckeTower[3](rvUnramElem[K]{4});
         rvunram F[7,1]{0,1}; rvDworkHeckeNeighbors[9](rvUnramElem[F]{1});",
    ));
    let diagonal = extra.lines().next().unwrap();
    assert!(diagonal.contains("cyclic diagonal multiplicity=2"));
    assert!(diagonal.contains("(2, 2,"));
    assert!(diagonal.contains("(6, 1,"));
    assert!(extra
        .trim_end()
        .ends_with("rational roots (value, multiplicity)=[([1], 4), ([5], 8)]"));

    let extension = stdout(interpret(
        "rvunram W[7,1]{1,0,1}; let x=rvUnramElem[W]{0,1}; rvDworkHeckeTower[5](x);",
    ));
    assert!(
        extension.contains("resultant degree=36, cyclic degree=30, backtracking divisor degree=6")
    );
    assert!(extension.contains("(20, 1,"));
    assert!(extension
        .trim_end()
        .ends_with("rational roots (value, multiplicity)=[([0, 6], 2)]"));
}

#[test]
fn cyclic_and_backtracking_factors_in_characteristic_three_are_schematic() {
    let out = stdout(example(EXAMPLE));
    assert!(out.contains("backtracking divisor degree=4, cyclic diagonal multiplicity=2, level characteristic with schematic multiplicities"));
    assert!(out.contains(
        "rational roots (value, multiplicity)=[([0, 1, 0], 2), ([1, 1, 0], 9), ([2, 1, 0], 1)]"
    ));
    assert!(out.contains("cyclic level in characteristic three, schematic multiplicities"));
}

#[test]
fn constructs_a_cyclic_degree_nine_map_between_distinct_modular_parameters() {
    let out = stdout(example(EXAMPLE));
    assert!(out.contains("rvECMapDegree = 9"));
    assert!(out.contains("rvECMapCertify = degree=9, cleared identity=yes"));
}

#[test]
fn tower_rejects_unimplemented_global_levels_and_wrong_domains() {
    for (source, fragment) in [
        ("rvDworkHeckeEquation[25];", "support good levels 3 and 5"),
        ("rvunram F[41,1]{0,1}; rvDworkHeckeTower[2,3](rvUnramElem[F]{1});", "support good levels 3 and 5"),
        ("rvunram F[3,2]{0,1}; show rvDworkHeckeCyclicCoeff[3,3,1](rvUnramElem[F]{1});", "odd residue field at precision one"),
        ("rvunram F[41,1]{0,1}; rvunram G[41,1]{0,1}; show rvDworkHeckeCompose[3](rvUnramElem[F]{1},rvUnramElem[G]{1});", "rvDworkHeckeCompose"),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(fragment), "{error}");
    }
}
