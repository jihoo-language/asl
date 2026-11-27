mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn ambient_trace_and_character_are_nontrivial_off_the_critical_locus() {
    let output = stdout(example("deligne-quiver-potential-exponentials"));
    for line in [
        "RVPartitionQuiver(C,[2])[[1, 6, 0, 3]]",
        "rvPartitionQuiverCritical[C] = false",
        "rvPartitionQuiverTrace[C] = 5",
        "rvPartitionQuiverTrace[C] = 2",
        "rvPartitionQuiverTrace[C] = 0",
        "rvPartitionQuiverDistribution[C] = q=7, framed=[1, 0, 3, 0, 0, 3, 0], total=7, groupOrder=6, stack=[1/6, 0/1, 1/2, 0/1, 0/1, 1/2, 0/1]",
        "rvPartitionQuiverFourierCheck[C] = q=7, recovered=[1, 0, 3, 0, 0, 3, 0]",
        "rvPartitionQuiverDistribution[B] = q=3, framed=[1953, 1440, 3168], total=6561, groupOrder=2304, stack=[217/256, 5/8, 11/8]",
        "rvPartitionQuiverDistribution[E] = q=4, framed=[1, 3, 0, 0], total=4, groupOrder=3, stack=[1/3, 1/1, 0/1, 0/1]",
    ] {assert!(output.contains(line),"missing {line}:\n{output}");}
}

#[test]
fn exact_primary_type_distributions_agree_with_raw_matrices() {
    let text = stdout(interpret(concat!(
        "rvpartitionginzburg A[1,0,5]; rvpartitionginzburg B[1,0,2];\n",
        "rvPartitionQuiverDistribution[A,structure,[2]];\n",
        "rvPartitionQuiverDistribution[A,direct,[2]];\n",
        "rvPartitionQuiverDistribution[B,structure,[3]];\n",
        "rvPartitionQuiverDistribution[B,direct,[3]];\n",
        "rvunram F[2,1]{1,1,1}; rvpartitionginzburg E[1,0,F];\n",
        "rvPartitionQuiverDistribution[E,structure,[2]];\n",
        "rvPartitionQuiverDistribution[E,direct,[2]];\n",
        "rvPartitionQuiverFourierCheck[E,structure,[2]];\n",
        "assert rvPartitionQuiverFourier[A,auto,1,[2]]==rvCycloScale[100](rvCyclo[5,0]);\n",
        "assert rvPartitionQuiverFourier[A,structure,0,[2]]==rvCycloScale[625](rvCyclo[5,0]);\n",
    )));
    let rows: Vec<_> = text.lines().collect();
    assert_eq!(rows[0], rows[1]);
    assert_eq!(rows[2], rows[3]);
    assert_eq!(rows[4], rows[5]);
    assert!(rows[0].contains("framed=[205, 105, 105, 105, 105], total=625, groupOrder=480, stack=[41/96, 7/32, 7/32, 7/32, 7/32]"));
    assert!(rows[6].starts_with("rvPartitionQuiverFourierCheck[E] = q=4, recovered="));
}

#[test]
fn rank_fibres_agree_with_raw_two_cycle_matrices_in_prime_and_named_fields() {
    let text=stdout(interpret(concat!(
        "rvpartitionginzburg A[2,0,3]; rvpartitionginzburg B[2,0,5];\n",
        "rvPartitionQuiverDistribution[A,structure,[2,0,2]];\n",
        "rvPartitionQuiverDistribution[A,direct,[2,0,2]];\n",
        "rvPartitionQuiverDistribution[B,structure,[2,0,2]];\n",
        "rvPartitionQuiverDistribution[B,direct,[2,0,2]];\n",
        "rvunram F[3,1]{1,0,1}; rvpartitionginzburg C[2,0,F];\n",
        "rvPartitionQuiverDistribution[C,structure,[1,0,1]];\n",
        "rvPartitionQuiverDistribution[C,direct,[1,0,1]];\n",
        "rvPartitionQuiverFourierCheck[C,structure,[1,0,1]];\n",
        "let u=rvPartitionQuiverRep[C,[1,0,1]]{[3],[1]};\n",
        "assert rvPartitionQuiverCharacter[1](u)==rvCyclo[3,2];\n",
        "assert rvPartitionQuiverCharacter[3](u)==rvCyclo[3,0];\n",
        "let x=rvPartitionQuiverRep[A,[1,0,1]]{[1],[1]};\n",
        "let y=rvPartitionQuiverRep[A,[2,0,1]]{[1,0],[1,1]};\n",
        "let z=rvPartitionQuiverTriangular(x,y){[1,2],[2]};\n",
        "assert rvPartitionQuiverCharacter[1](z)==rvCyclo[3,1];\n",
        "assert rvPartitionQuiverFourier[A,structure,1,[2,0,2]]==neg(add(rvCycloScale[1215](rvCyclo[3,0]),rvCycloScale[1728](rvCyclo[3,1])));\n",
    )));
    let rows: Vec<_> = text.lines().collect();
    assert_eq!(rows[0], rows[1]);
    assert_eq!(rows[2], rows[3]);
    assert_eq!(rows[4], rows[5]);
    assert!(rows[2]
        .contains("framed=[88225, 57600, 93600, 93600, 57600], total=390625, groupOrder=230400"));
    assert!(rows[4].contains("framed=[17, 16, 16, 16, 0, 0, 16, 0, 0], total=81, groupOrder=64"));
    assert!(rows[6].ends_with("recovered=[17, 16, 16, 16, 0, 0, 16, 0, 0]"));
}

#[test]
fn structured_counts_scale_beyond_ambient_enumeration_and_include_isolated_vertices() {
    let output = stdout(interpret(concat!(
        "rvpartitionginzburg A[2,0,3]; rvpartitionginzburg L[1,0,5];\n",
        "rvPartitionQuiverDistribution[A,structure,[5,0,5]];\n",
        "rvPartitionQuiverDistribution[A,structure,[0,2,5]];\n",
        "rvPartitionQuiverDistribution[L,structure,[4]];\n",
        "rvPartitionQuiverDistribution[L,structure,[0]];\n",
        "rvpartitionginzburg S[2,3,2]; rvPartitionQuiverDistribution[S,structure,[0,2,0,0]];\n",
    )));
    assert!(output.contains("framed=[239325499508521555281285, 239343322457739366979200, 239229165725591666509764], total=717897987691852588770249"));
    assert!(output.contains("framed=[1, 0, 0], total=1, groupOrder=22827190763520"));
    assert!(output.contains("framed=[30672578125, 30478828125, 30478828125, 30478828125, 30478828125], total=152587890625"));
    assert!(output.contains("framed=[1, 0, 0, 0, 0], total=1, groupOrder=1"));
    assert!(output.contains("framed=[1, 0], total=1, groupOrder=6"));
}

#[test]
fn invalid_potential_domain_nominality_and_budget_errors_are_specific() {
    let fixture="rvpartitionginzburg D[1,0,7]; rvpartitionginzburg E[1,0,7]; let x=rvPartitionQuiverRep[D,[1]]{[1]}; let y=rvPartitionQuiverRep[E,[1]]{[1]};";
    for (suffix,expected) in [
        ("show rvPartitionQuiverRep[D,[1,0]]{[0]};","one dimension per vertex"),
        ("show rvPartitionQuiverRep[D,[65]]{[0]};","total dimension exceeds 64"),
        ("show rvPartitionQuiverRep[D,[1]]{};","one matrix entry list per original arrow"),
        ("show rvPartitionQuiverRep[D,[2]]{[0]};","entry count does not match"),
        ("show rvPartitionQuiverTriangular(x,y){[0]};","same nominal Ginzburg quiver"),
        ("show rvPartitionQuiverTriangular(x,x){[0,1]};","entry count does not match"),
        ("show rvPartitionQuiverCharacter[7](x);","outside the finite field codes"),
        ("show rvPartitionQuiverFourier[D,auto,7,[1]];","outside the finite field codes"),
        ("rvPartitionQuiverDistribution[D,anything,[1]];","count method must be auto, structure, or direct"),
        ("rvPartitionQuiverDistribution[D,direct,[4]];","500,000 matrix assignments"),
        ("rvPartitionQuiverDistribution[D,structure,[16]];","exceed unsigned 128-bit arithmetic"),
        ("rvpartitionginzburg W[1,0,257];rvPartitionQuiverDistribution[W,auto,[1]];","field of order at most 256"),
        ("rvpartitionginzburg W[1,0,3];","potential denominator 3 vanishes"),
        ("rvpartitionginzburggradient W[1,0,3];rvPartitionQuiverTrace(rvPartitionQuiverRep[W,[1]]{[0]});","integral differentiated data only"),
        ("rvpartitionginzburggradient W[2,0,2];rvPartitionQuiverDistribution[W,structure,[1,0,1]];","integral differentiated data only"),
    ] {
        let error=diagnostic(interpret(&format!("{fixture} {suffix}")));
        assert!(error.contains(expected),"{suffix}: {error}");
    }
}
