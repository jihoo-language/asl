mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn example_mutates_twice_and_checks_a_nonsimple_brick_in_both_dual_bases() {
    let output = stdout(example("deligne-support-tau-tilting"));
    assert_eq!(output, concat!(
        "rvPartitionTauTilting[2,0] = faces=[1, 8, 18, 12], chambers=12, torsionClasses=12, covers=18\n",
        "rvPartitionTauSummands = [M(1,1), M(2,1), P(0)[1]]\n",
        "rvPartitionTauVectors = g=[[0, 1, 0], [-1, 0, 1], [-1, 0, 0]], c=[[0, 1, 0], [0, 0, 1], [-1, 0, -1]]\n",
        "rvPartitionTauVectors = g=[[0, 1, 0], [-1, 0, 1], [0, 0, 1]], c=[[0, 1, 0], [-1, 0, 0], [1, 0, 1]]\n",
        "rvPartitionTauWall = cover=(3, 2), brick=M(2,2), normal=[1, 0, 1], signs=[-1, 1]\n",
        "rvPartitionTauWall = cover=(7, 3), brick=M(1,1), normal=[0, 1, 0], signs=[-1, 1]\n",
        "rvPartitionTauInclusion = true\n",
        "rvPartitionTauGreen[2,0] = length=4, chambers=[11, 7, 3, 2, 1], bricks=[[2, 1], [1, 1], [2, 2], [0, 1]]\n",
    ));
}

#[test]
fn all_fixtures_match_face_counts_unique_torsions_and_shortest_chains() {
    let output = stdout(interpret(concat!(
        "rvPartitionTauTilting[2,1];rvPartitionTauTilting[2,2];rvPartitionTauTilting[2,3];",
        "rvPartitionTauTilting[1,0];rvPartitionTauTilting[0,0];",
        "rvPartitionTauGreen[2,1];rvPartitionTauGreen[2,2];rvPartitionTauGreen[1,0];",
    )));
    assert_eq!(output,concat!(
        "rvPartitionTauTilting[2,1] = faces=[1, 10, 34, 48, 24], chambers=24, torsionClasses=24, covers=48\n",
        "rvPartitionTauTilting[2,2] = faces=[1, 10, 34, 48, 24], chambers=24, torsionClasses=24, covers=48\n",
        "rvPartitionTauTilting[2,3] = faces=[1, 8, 24, 32, 16], chambers=16, torsionClasses=16, covers=32\n",
        "rvPartitionTauTilting[1,0] = faces=[1, 2], chambers=2, torsionClasses=2, covers=1\n",
        "rvPartitionTauTilting[0,0] = faces=[1, 2], chambers=2, torsionClasses=2, covers=1\n",
        "rvPartitionTauGreen[2,1] = length=5, chambers=[23, 15, 11, 5, 2, 1], bricks=[[3, 1], [3, 2], [1, 1], [0, 1], [2, 1]]\n",
        "rvPartitionTauGreen[2,2] = length=5, chambers=[23, 11, 7, 5, 2, 1], bricks=[[3, 1], [2, 1], [2, 2], [0, 1], [1, 1]]\n",
        "rvPartitionTauGreen[1,0] = length=1, chambers=[1, 0], bricks=[[0, 1]]\n",
    ));
}

#[test]
fn fac_and_torsion_lattice_distinguish_quotients_extensions_and_isolated_blocks() {
    let output = stdout(interpret(concat!(
        "let a=rvPartitionTauFac(rvPartitionTauModule[2,0,0,1]);\n",
        "let b=rvPartitionTauFac(rvPartitionTauModule[2,0,2,1]);\n",
        "let isolated=rvPartitionTauFac(rvPartitionTauModule[2,0,1,1]);\n",
        "let linked=rvPartitionTauJoin(a,b);\n",
        "let bottom=rvPartitionTauClass[2,0,11];\n",
        "assert rvPartitionTauMeet(a,b)==bottom;\n",
        "assert rvPartitionTauMeet(linked,isolated)==bottom;\n",
        "assert rvPartitionTauJoin(linked,isolated)==rvPartitionTauClass[2,0,1];\n",
        "assert linked==rvPartitionTauClass[2,0,5];\n",
        "assert rvPartitionTauRecover(linked)==rvPartitionTauChamber[2,0,5];\n",
        "rvPartitionTauMembers(a);\n",
        "rvPartitionTauMembers(linked);\n",
        "rvPartitionTauContains[0,2](a);\n",
        "rvPartitionTauContains[0,2](linked);\n",
        "rvPartitionTauContains[1,1](linked);\n",
        "rvPartitionTauInclusion(a,linked);\n",
        "rvPartitionTauInclusion(linked,b);\n",
    )));
    assert_eq!(
        output,
        concat!(
            "rvPartitionTauMembers = [[0, 1]]\n",
            "rvPartitionTauMembers = [[0, 1], [0, 2], [0, 3], [2, 1], [2, 2], [2, 3]]\n",
            "rvPartitionTauContains[0,2] = false\n",
            "rvPartitionTauContains[0,2] = true\n",
            "rvPartitionTauContains[1,1] = false\n",
            "rvPartitionTauInclusion = true\n",
            "rvPartitionTauInclusion = false\n",
        )
    );
}

#[test]
fn actual_finite_types_enumerate_faces_and_chambers_without_duplicates() {
    let output = stdout(interpret(
        r"
        card \x : RVPartitionTauPair(2,0) => S.hi;
        card \x : RVPartitionTauTorsion(2,0) => S.hi;
        card \x : RVPartitionTauPair(1,0) => S.hi;
        let x=rvPartitionTauCombine(rvPartitionTauModule[2,0,2,1],rvPartitionTauShift[2,0,0]);
        let y=rvPartitionTauCombine(x,rvPartitionTauModule[2,0,1,1]);
        assert y == rvPartitionTauChamber[2,0,3];
        assert rvPartitionTauCombine(y,y)==y;
        assert rvPartitionTauFace[2,0,0,0]==rvPartitionTauEmpty[2,0];
        assert rvPartitionTauMutate[2](rvPartitionTauMutate[2](y))==y;
        rvPartitionTauCompatible(rvPartitionTauModule[2,0,0,1],rvPartitionTauModule[2,0,2,1]);
        rvPartitionTauCompatible(x,rvPartitionTauModule[2,0,1,1]);
        rvPartitionTauVectors(rvPartitionTauShift[2,0,0]);
    ",
    ));
    assert_eq!(
        output,
        concat!(
            "card = 39 / 39\n",
            "card = 12 / 12\n",
            "card = 3 / 3\n",
            "rvPartitionTauCompatible = false\n",
            "rvPartitionTauCompatible = true\n",
            "rvPartitionTauVectors = g=[[-1, 0, 0]]\n",
        )
    );
}

#[test]
fn every_mutation_factors_through_the_reported_codimension_one_face() {
    let mut source = String::new();
    for (n, rank, chambers) in [(0, 3, 12), (1, 4, 24), (2, 4, 24)] {
        for j in 0..chambers {
            for k in 0..rank {
                let chamber = format!("rvPartitionTauChamber[2,{n},{j}]");
                writeln!(
                    source,
                    "let chosen{n}{j}_{k}=rvPartitionTauMutate[{k}]({chamber});"
                )
                .unwrap();
                writeln!(source, "rvPartitionTauWall({chamber},chosen{n}{j}_{k});").unwrap();
                writeln!(source,"rvPartitionTauCompatible(rvPartitionTauDrop[{k}]({chamber}),chosen{n}{j}_{k});").unwrap();
            }
        }
    }
    let output = stdout(interpret(&source));
    assert_eq!(output.lines().count(), 2 * (36 + 96 + 96));
    for (k, line) in output.lines().enumerate() {
        if k.is_multiple_of(2) {
            assert!(
                line.contains("brick=M(") && line.ends_with("signs=[-1, 1]"),
                "{line}"
            );
        } else {
            assert_eq!(line, "rvPartitionTauCompatible = true");
        }
    }
}

#[test]
fn wrong_fibers_nonrigid_summands_and_nonadjacent_chambers_are_checked() {
    for (source,expected) in [
        ("rvPartitionTauTilting[3,0];","limited to rank 2"),
        ("show rvPartitionTauModule[2,0,0,2];","module is not tau-rigid"),
        ("show rvPartitionTauShift[2,0,3];","vertex is out of range"),
        ("show rvPartitionTauFace[2,0,3,12];","face index is out of range"),
        ("show rvPartitionTauMutate[0](rvPartitionTauEmpty[2,0]);","requires a maximal tau-tilting chamber"),
        ("show rvPartitionTauDrop[3](rvPartitionTauChamber[2,0,0]);","summand index is out of range"),
        ("show rvPartitionTauCombine(rvPartitionTauModule[2,0,0,1],rvPartitionTauModule[2,0,2,1]);","incompatible summands"),
        ("show rvPartitionTauCombine(rvPartitionTauEmpty[2,0],rvPartitionTauEmpty[2,1]);","same rank and integral parameter"),
        ("rvPartitionTauWall(rvPartitionTauChamber[2,0,1],rvPartitionTauChamber[2,0,11]);","do not share a wall"),
        ("rvPartitionTauInclusion(rvPartitionTauEmpty[2,0],rvPartitionTauClass[2,0,11]);","requires an RVPartitionTauTorsion"),
        ("rvPartitionTauContains[2,4](rvPartitionTauClass[2,0,1]);","module top or length is out of range"),
    ] {
        let message=diagnostic(interpret(source));
        assert!(message.contains(expected),"{source}: {message}");
    }
}
