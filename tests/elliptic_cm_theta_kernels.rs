mod common;

use common::{diagnostic, example, interpret, stdout};

const F41: &str = include_str!("../examples/elliptic-cm-theta-kernels.asl");

#[test]
fn nonprincipal_polarization_and_geometric_three_torsion_agree() {
    let output = stdout(example("elliptic-cm-theta-kernels"));
    for expected in [
        "matrix=[[3, 0, 1, 3], [0, 3, -1, 0], [0, -3, 3, 0], [1, 1, 0, 3]], smith=[1, 1, 6, 6], type=[1, 6], dimension=6, kernel_size=36",
        "rvCMThetaLagrangians = count=12, by_exponent={6: 12}",
        "rvCMThetaPair = 1/6",
        "rvCMThetaPair = 5/6",
        "rvCMThetaOrder = 6",
        "exponent=6, isotropic=true, perpendicular_size=6, maximal_isotropic=true, geometric_descent_degree=6",
        "rvCMThetaAction = position=[0, 3], phase=4, root_order=6",
        "rvCMThetaLagrangians = count=15, by_exponent={2: 15}",
        "rvCMThetaLagrangians = count=39, by_exponent={2: 3, 4: 36}",
        "level=3, matrix=[[0, 0, 2, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 2, 0, 0]], kernel_smith=[1, 1, 3, 3]",
    ] {
        assert!(output.contains(expected), "{output}");
    }
}

#[test]
fn original_lattice_pairing_matches_inverse_matrix_and_annihilates_relations() {
    let output = stdout(interpret(&format!(
        "{F41}
        let u=rvCMThetaFromLattice[Theta]{{1,0,0,0}};
        let v=rvCMThetaFromLattice[Theta]{{0,1,0,0}};
        let zeta=rvCMThetaFromLattice[Theta]{{0,0,1,0}};
        let t=rvCMThetaFromLattice[Theta]{{0,0,0,1}};
        let z=zero[RVCMTheta(Theta)];
        assert add(zeta,neg(zeta))==z;
        assert rvCMThetaFromLattice[Theta]{{3,0,0,1}}==z;
        assert rvCMThetaFromLattice[Theta]{{0,3,-3,1}}==z;
        assert rvCMThetaFromLattice[Theta]{{1,-1,3,0}}==z;
        assert rvCMThetaFromLattice[Theta]{{3,0,0,3}}==z;
        rvCMThetaPair(u,v);
        rvCMThetaPair(u,zeta);
        rvCMThetaPair(u,t);
        rvCMThetaPair(v,zeta);
        rvCMThetaPair(v,t);
        rvCMThetaPair(zeta,t);
        rvCMThetaPair(neg(zeta),v);
        rvCMThetaPair(zeta,zeta);
        rvCMThetaPair(u,add(v,zeta));
        rvCMThetaElements[Theta];"
    )));
    let pairings: Vec<_> = output
        .lines()
        .filter_map(|line| line.strip_prefix("rvCMThetaPair = "))
        .skip(2)
        .collect();
    assert_eq!(
        pairings,
        ["1/2", "1/6", "0/1", "1/6", "1/2", "1/2", "1/6", "0/1", "2/3"]
    );
    let line = output
        .lines()
        .find(|line| line.starts_with("rvCMThetaElements"))
        .unwrap();
    assert_eq!(line.matches("[0, 0, ").count(), 36);
}

#[test]
fn subgroup_perpendiculars_nonfree_cases_and_heisenberg_phases_are_exact() {
    let output = stdout(interpret(&format!(
        "{F41}
        let a=rvCMThetaElem[Mixed]{{1,0,0,0}};
        let b=rvCMThetaElem[Mixed]{{0,0,2,0}};
        let c=rvCMThetaElem[Mixed]{{0,0,0,2}};
        let d=rvCMThetaElem[Mixed]{{0,1,1,0}};
        let s=rvCMThetaSubgroup[Mixed]{{a,b,c}};
        assert s==rvCMThetaSubgroup[Mixed]{{c,b,a,c}};
        assert rvCMThetaPerp(s)==s;
        rvCMThetaSubgroupData(s);
        let nonsplit=rvCMThetaSubgroup[Theta]{{e,ftheta}};
        assert rvCMThetaPerp(rvCMThetaPerp(nonsplit))==nonsplit;
        rvCMThetaSubgroupData(nonsplit);
        rvCMThetaSubgroupData(rvCMThetaLagrangian[Mixed,38]);
        let A1=rvCMThetaOp[0](a);
        let B1=rvCMThetaOp[2](b);
        let C1=rvCMThetaOp[5](d);
        assert rvCMThetaMul(rvCMThetaMul(A1,B1),C1)==rvCMThetaMul(A1,rvCMThetaMul(B1,C1));
        let central=rvCMThetaOp[2](zero[RVCMTheta(Mixed)]);
        assert rvCMThetaMul(A1,C1)==rvCMThetaMul(central,rvCMThetaMul(C1,A1));
        rvCMThetaAction[1,3](rvCMThetaMul(C1,A1));"
    )));
    assert!(output.contains(
        "size=8, exponent=2, isotropic=true, perpendicular_size=8, maximal_isotropic=true"
    ));
    assert!(output.contains("size=36, exponent=6, isotropic=false, perpendicular_size=1, maximal_isotropic=false, geometric_descent_degree=unavailable"));
    assert!(output.contains("rvCMThetaAction = position=[0, 0], phase=1, root_order=4"));
}

#[test]
fn geometric_level_check_compensates_for_a_changed_primitive_root() {
    let output = stdout(interpret(&format!(
        "{F41}
         rvjactorsionbasis BFflip[FL,3](rvJacWeil[3](R3,T)){{R3,T}};
         rvCMThetaTorsion[Theta,BE,BFflip];
         rvCMThetaTorsion[Double,BE,BFflip];"
    )));
    assert!(
        output.contains("kernel_smith=[1, 1, 3, 3], right_pairing_scale=2"),
        "{output}"
    );
    assert!(
        output.contains("kernel_smith=[1, 1, 1, 1], right_pairing_scale=2"),
        "{output}"
    );
}

#[test]
fn integral_model_survives_characteristic_dividing_degree() {
    let output = stdout(interpret(&format!(
        "{F41}
         rvcmthetakernel Char(rvCMAbelianDiv[A,1,41](off));
         rvCMThetaData[Char];
         let line=rvCMThetaSubgroup[Char]{{rvCMThetaElem[Char]{{0,0,1,0}}}};
         rvCMThetaSubgroupData(line);
         rvcmthetakernel Principal(rvCMAbelianDiv[A,1,1](off));
         rvCMThetaLagrangians[Principal];
         rvCMThetaData[Principal];"
    )));
    for expected in [
        "type=[1, 41], dimension=41, kernel_size=1681, characteristic=41, geometric_etale=false",
        "size=41, exponent=41, isotropic=true, perpendicular_size=41, maximal_isotropic=true, geometric_descent_degree=unavailable",
        "rvCMThetaLagrangians = count=1, by_exponent={1: 1}",
        "type=[1, 1], dimension=1, kernel_size=1, characteristic=41, geometric_etale=true",
    ] {
        assert!(output.contains(expected), "{output}");
    }
}

#[test]
fn nominal_integrality_positivity_and_enumeration_guards() {
    for (query, expected) in [
        (
            "rvcmthetakernel Bad(rvCMAbelianScale[1,2](D));",
            "requires an integral divisor",
        ),
        (
            "rvcmthetakernel Bad(neg(D));",
            "requires a positive definite divisor",
        ),
        (
            "rvcmthetakernel Bad(rvCMAbelianDiv[A,1,1](rvCMHomCoords[H](g)));",
            "requires a positive definite divisor",
        ),
        (
            "rvcmthetakernel Bad(rvCMAbelianMukai[1,1](D));",
            "expects a marked abelian divisor",
        ),
        (
            "rvCMThetaPair(e,zero[RVCMTheta(Mixed)]);",
            "CM theta pairing argument",
        ),
        (
            "show rvCMThetaMul(op,rvCMThetaOp[0](zero[RVCMTheta(Mixed)]));",
            "CM theta operator argument",
        ),
        (
            "show rvCMThetaSubgroup[Double]{e};",
            "CM theta subgroup generator",
        ),
        ("rvCMThetaAction[1,0](op);", "outside the Schrödinger basis"),
        (
            "show rvCMThetaLagrangian[Theta,12];",
            "Lagrangian index is out of range",
        ),
        (
            "rvcmthetakernel Huge(rvCMAbelianDiv[A,3,61](off)); rvCMThetaElements[Huge];",
            "4096-element enumeration limit",
        ),
        (
            "rvCMThetaTorsion[Theta,BE,BE];",
            "does not give the declared endpoint curves",
        ),
        ("rvcmthetakernel Theta(D);", "already defined"),
    ] {
        let error = diagnostic(interpret(&format!("{F41}\n{query}")));
        assert!(error.contains(expected), "{error}");
    }
}

#[test]
fn odd_fundamental_discriminant_uses_the_half_integral_maximal_order_basis() {
    let output = stdout(interpret(
        r"
        rvunram k[7,1]{0,1};
        rvec E[k](rvUnramElem[k]{0},rvUnramElem[k]{3});
        rvec F[k](rvUnramElem[k]{0},rvUnramElem[k]{3});
        let f=rvECChart[E,F](rvUnramElem[k]{1});
        let g=rvECChart[E,F](rvUnramElem[k]{2});
        let fd=rvECChart[F,E](rvUnramElem[k]{1});
        let gd=rvECChart[F,E](rvUnramElem[k]{4});
        rvcmhomframe H(f,g);
        rvcmhomdual HD[H](fd,gd);
        rvcmhomframe RE(rvECMultMap[E,1],rvECMapCompose(fd,g));
        rvcmhomframe RF(rvECMultMap[F,1],rvECMapCompose(g,fd));
        rvcmabelian A[H,HD,RE,RF];
        rvcmthetakernel K(rvCMAbelianDiv[A,2,2](rvCMHomCoord[H]{0,1}));
        rvCMThetaData[K];
        rvCMThetaLagrangians[K];",
    ));
    assert!(
        output.contains(
            "matrix=[[2, 0, 0, 1], [0, 2, -1, 0], [0, -1, 2, 0], [1, 0, 0, 2]], smith=[1, 1, 3, 3]"
        ),
        "{output}"
    );
    assert!(output.contains("count=4, by_exponent={3: 4}"), "{output}");
}
