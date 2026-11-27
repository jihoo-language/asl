mod common;

use common::{diagnostic, example, interpret, stdout};

const SINGLE: &str = r"
    rvunram F[3,1]{0,1};
    rvunram W[3,1]{1,0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
    rvjacbase L[C,W];
    let zeroW=rvUnramElem[W]{};
    let oneW=rvUnramElem[W]{1};
    let twoW=neg(oneW);
    let a=rvUnramElem[W]{0,1};
    let r0=rvJacPoint[L](zeroW,zeroW);
    let r1=rvJacPoint[L](oneW,zeroW);
    let r2=rvJacPoint[L](twoW,zeroW);
    let r3=rvJacPoint[L](a,zeroW);
    let e1=r0;
    let e2=add(add(r0,r1),r2);
    let f1=r1;
    let f2=add(add(r0,r1),r3);
    rvjactorsionbasis B[L,2](twoW){e1,e2,f1,f2};
";

#[test]
fn transposition_fixes_three_planes_and_constructs_the_actual_smooth_quotient() {
    assert_eq!(
        stdout(example("hyperelliptic-richelot-hecke")),
        concat!(
            "rvJacRichelotKernels[B, C] = orbits=[1, 1, 1, 2, 2, 2, 2, 2, 2], stable=[",
            "{masks=[1, 2, 3], rational=4, determinant_zero=false}, ",
            "{masks=[2, 4, 6], rational=4, determinant_zero=true}, ",
            "{masks=[2, 5, 7], rational=4, determinant_zero=true}]\n",
            "rvRichelotData[D] = source=L, determinant=[1, 0], f=[[0, 0], [2, 0], [0, 0], [0, 0], [0, 0], [1, 0]], ",
            "factors=[[[0, 0], [2, 0]], [[1, 0], [0, 0], [1, 0]], [[1, 0], [0, 0], [2, 0]]]\n"
        )
    );
}

#[test]
fn identity_frobenius_lists_all_fifteen_planes_in_lexicographic_order() {
    assert_eq!(
        stdout(interpret(&format!("{SINGLE} rvJacRichelotKernels[B,L];"))),
        concat!(
            "rvJacRichelotKernels[B, L] = orbits=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], stable=[",
            "{masks=[1, 2, 3], rational=4, determinant_zero=false}, ",
            "{masks=[1, 8, 9], rational=4, determinant_zero=true}, ",
            "{masks=[1, 10, 11], rational=4, determinant_zero=true}, ",
            "{masks=[2, 4, 6], rational=4, determinant_zero=true}, ",
            "{masks=[2, 5, 7], rational=4, determinant_zero=true}, ",
            "{masks=[3, 12, 15], rational=4, determinant_zero=true}, ",
            "{masks=[3, 13, 14], rational=4, determinant_zero=true}, ",
            "{masks=[4, 8, 12], rational=4, determinant_zero=false}, ",
            "{masks=[4, 10, 14], rational=4, determinant_zero=false}, ",
            "{masks=[5, 8, 13], rational=4, determinant_zero=false}, ",
            "{masks=[5, 10, 15], rational=4, determinant_zero=false}, ",
            "{masks=[6, 9, 15], rational=4, determinant_zero=false}, ",
            "{masks=[6, 11, 13], rational=4, determinant_zero=false}, ",
            "{masks=[7, 9, 14], rational=4, determinant_zero=false}, ",
            "{masks=[7, 11, 12], rational=4, determinant_zero=false}]\n"
        )
    );
}

#[test]
fn two_disjoint_transpositions_can_fix_a_kernel_without_fixing_its_generators() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram F[3,1]{0,1};
            rvunram W[3,1]{1,0,1};
            let z=rvUnramElem[F]{};
            let o=rvUnramElem[F]{1};
            rvhyperelliptic C[F]{z,neg(o),o,z,o,o};
            rvjacbase L[C,W];
            let zeroW=rvUnramElem[W]{};
            let oneW=rvUnramElem[W]{1};
            let a=rvUnramElem[W]{0,1};
            let r0=rvJacPoint[L](zeroW,zeroW);
            let r1=rvJacPoint[L](a,zeroW);
            let r2=rvJacPoint[L](neg(a),zeroW);
            let r3=rvJacPoint[L](add(oneW,a),zeroW);
            rvjactorsionbasis B[L,2](neg(oneW)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),r3)};
            rvJacRichelotKernels[B,C];
            let k0=rvJacRichelotKernel[B,C,1,0];
            let k1=rvJacRichelotKernel[B,C,1,1];
            let k2=rvJacRichelotKernel[B,C,1,2];
            assert k0==r0;
            assert rvJacFrob[C,1](k0)==k0;
            assert rvJacFrob[C,1](k1)==k2;
            assert rvJacFrob[C,1](k2)==k1;
            assert k2==add(k0,k1);
            rvrichelotfromlevel D[B,C,1];
            rvjacchart N[D](zeroW);
            assert rvRichelotPush[N](k0)==zero[RVJac(N)];
            assert rvRichelotPush[N](k1)==zero[RVJac(N)];
            assert rvRichelotPush[N](k2)==zero[RVJac(N)];
            "
        )),
        concat!(
            "rvJacRichelotKernels[B, C] = orbits=[1, 1, 1, 2, 2, 2, 2, 2, 2], stable=[",
            "{masks=[1, 2, 3], rational=4, determinant_zero=false}, ",
            "{masks=[1, 8, 9], rational=2, determinant_zero=false}, ",
            "{masks=[1, 10, 11], rational=2, determinant_zero=false}]\n"
        )
    );
}

#[test]
fn a_level_target_descends_when_frobenius_exchanges_its_branch_pairs() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram F[3,1]{0,1};
            rvunram W[3,1]{1,0,1};
            let z=rvUnramElem[F]{};
            let o=rvUnramElem[F]{1};
            let t=neg(o);
            rvhyperelliptic C[F]{z,t,t,z,t,o};
            rvjacbase L[C,W];
            let zeroW=rvUnramElem[W]{};
            let twoW=rvUnramElem[W]{2};
            let a=rvUnramElem[W]{0,1};
            let r0=rvJacPoint[L](zeroW,zeroW);
            let r1=rvJacPoint[L](neg(a),zeroW);
            let r2=rvJacPoint[L](a,zeroW);
            let r3=rvJacPoint[L](add(twoW,a),zeroW);
            rvjactorsionbasis B[L,2](twoW){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),r3)};
            rvrichelotfromlevel DL[B,C,1];
            rvrichelotdescend D[C,DL];
            rvRichelotDescent[D];
            rvhyperelliptic Known[F]{o,t,t,o,t,o};
            assert rvJacExtend[Known](zero[RVJac(D)])==zero[RVJac(Known)];
            assert (\p:RVJac(C)=>rvJacExtend[DL](rvRichelotDescPush[D](p)))
                ==(\p:RVJac(C)=>rvRichelotPush[DL](rvJacExtend[L](p)));
            assert rvRichelotDescPush[D](rvJacPoint[C](z,z))==zero[RVJac(D)];
            rvJacOrder(rvRichelotDescPush[D](rvJacPoint[C](o,o)));
            "
        )),
        concat!(
            "rvRichelotDescent[D] = relative degree=2, permutation=[0, 2, 1], rational kernel order=2\n",
            "rvJacOrder = 10\n"
        )
    );
}

#[test]
fn a_single_three_cycle_has_no_stable_perfect_matching() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram F[3,1]{0,1};
            rvunram W[3,1]{1,2,0,1};
            let z=rvUnramElem[F]{};
            let o=rvUnramElem[F]{1};
            let t=neg(o);
            rvhyperelliptic C[F]{z,t,t,t,t,o};
            rvjacbase L[C,W];
            let zeroW=rvUnramElem[W]{};
            let oneW=rvUnramElem[W]{1};
            let a=rvUnramElem[W]{0,1};
            let r0=rvJacPoint[L](zeroW,zeroW);
            let r1=rvJacPoint[L](oneW,zeroW);
            let r2=rvJacPoint[L](a,zeroW);
            let r3=rvJacPoint[L](rvUnramFrob[1](a),zeroW);
            rvjactorsionbasis B[L,2](neg(oneW)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),r3)};
            rvJacRichelotKernels[B,C];
            "
        )),
        "rvJacRichelotKernels[B, C] = orbits=[3, 3, 3, 3, 3], stable=[]\n"
    );
}

#[test]
fn reports_invalid_level_data_split_quotients_and_bad_indices() {
    for (input, expected) in [
        ("rvrichelotfromlevel D[B,C,1];", "product of elliptic curves over the algebraic closure"),
        ("rvrichelotfromlevel D[B,C,3];", "stable Richelot plane index is out of range"),
        ("show rvJacRichelotKernel[B,C,3,0];", "stable Richelot plane index is out of range"),
        ("show rvJacRichelotKernel[B,C,0,3];", "kernel generator index is out of range"),
        ("rvhyperelliptic Wrong[F]{z,o,z,z,z,o}; rvJacRichelotKernels[B,Wrong];", "not the declared base extension"),
        ("rvJacRichelotKernels[Wrong,C];", "unknown Jacobian torsion basis"),
        ("show add(rvJacRichelotKernel[B,C,0,0],zero[RVJac(C)]);", "expected type RVJac(L)"),
        ("rvhyperelliptic E[F]{z,neg(o),z,o}; rvjactorsionbasis BE[E,2](neg(o)){rvJacPoint[E](z,z),rvJacPoint[E](o,z)}; rvJacRichelotKernels[BE,E];", "require a certified genus-two level-2 basis"),
        ("rvunram K[5,1]{3,0,1}; let zeroK=rvUnramElem[K]{}; let oneK=rvUnramElem[K]{1}; let alpha=rvUnramElem[K]{0,1}; rvhyperelliptic E[K]{oneK,zeroK,zeroK,oneK}; rvjactorsionbasis BE[E,3](rvUnramElem[K]{2,2}){rvJacPoint[E](zeroK,oneK),rvJacPoint[E](oneK,alpha)}; rvJacRichelotKernels[BE,E];", "require a certified genus-two level-2 basis"),
    ] {
        let problem = diagnostic(interpret(&format!("{SINGLE}{input}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
