mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn integral_split_and_filtered_kappa_have_typed_witnesses() {
    let output = stdout(example("annular-khovanov"));
    for line in [
        "khAnnularCheck[Mixed] = verified (36 enhanced states)",
        "khEuler[U] = [(-1, 1), (1, 1)]",
        "khEuler[P] = [(-1, 1), (1, 1)]",
        "khAnnularEuler[U] = [((-1, -1), 1), ((1, 1), 1)]",
        "khAnnularGroup[P, 0, 5, 3] = chain=1, free=1, torsion=[]",
        "khAnnularGroup[U, 0, 5, 3] = chain=0, free=0, torsion=[]",
        "khKappa[U] = infinity",
        "khKappa[P] = infinity",
        "khKappa[N] = 2",
        "khKappa[A] = 4",
        "khKappa[B] = 2",
        "khKappa[BPlus] = 4",
    ] {
        assert!(
            output.lines().any(|actual| actual == line),
            "{line}\n{output}"
        );
    }
}

#[test]
fn annular_smith_representatives_and_binary_duals_work_with_their_types() {
    let source = r"
        khbraid U[1,[]];
        khbraid N[2,[-1]];
        khbraid T[2,[1,1,1]];
        let class = khAGen[U,0,1,1,0];
        let chain = khABasis[U,0,1,1,0];
        assert khClass(chain) == class;
        assert khRep(class) == chain;
        assert add(class,neg(class)) == zero[KhAH(U,0,1,1)];
        assert khForgetK(chain) == khTerm[U,[],[]];
        assert khAnnular[-1](khTerm[U,[],[]]) == zero[KhAC(U,0,1,-1)];
        let cycle = khPsi[U];
        let dual = khKappaDual[U];
        assert khD(cycle) == zero[KhFC(U,1,-1,-1)];
        assert add(cycle,cycle) == zero[KhFC(U,0,-1,-1)];
        assert khD(dual) == zero[KhFD(U,-1,-1,1)];
        assert khPair(dual,khToF2[1](cycle)) == S.hi;
        assert khPair(dual,khToF2[1](khScale[2](cycle))) == S.lo;
        assert khPsi[N] == khToF2[-2](khTerm[N,[0],[0,1]]);
        khAnnularGroup[U,0,1,-9223372036854775808];
        let torsion = khAGen[T,3,7,0,0];
        assert khD0(khANull[T,3,7,0,0]) == khScale[2](khRep(torsion));
        assert add(torsion,torsion) == zero[KhAH(T,3,7,0)];
        khAnnularGroup[T,3,7,0];
        let paired = \x : KhFC(U,0,-1,-1) => khPair(dual,khToF2[1](x));
        card paired;
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "khAnnularGroup[U, 0, 1, -9223372036854775808] = chain=0, free=0, torsion=[]\n",
            "khAnnularGroup[T, 3, 7, 0] = chain=3, free=0, torsion=[2]\n",
            "card = 1 / 2\n",
        )
    );
}

#[test]
fn annular_frontend_checks_degrees_cutoffs_and_exhaustive_work() {
    for (source, expected) in [
        (
            "khbraid U[1,[]]; show khKappaNull[U];",
            "transverse cycle has no filtered null",
        ),
        (
            "khbraid U[1,[]]; show khPair(khKappaDual[U],khPsi[U]);",
            "requires the same braid, bidegree, and filtration cutoff",
        ),
        (
            "khbraid U[1,[]]; show khD0(khTerm[U,[],[]]);",
            "an annular Khovanov chain",
        ),
        (
            "khbraid U[1,[]]; show khDMinus2(zero[KhAC(U,0,1,-9223372036854775808)]);",
            "Khovanov annular degree is out of range",
        ),
        (
            "khbraid L[120,[]]; khAnnularCheck[L];",
            "enhanced cube exceeds the work limit",
        ),
    ] {
        let message = diagnostic(interpret(source));
        assert!(message.contains(expected), "{source}\n{message}");
    }
}
