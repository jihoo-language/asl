mod common;

use std::collections::HashSet;
use std::fmt::Write as _;

use common::{diagnostic, example, interpret, stdout};

const BASE: &str = "rvunram B[3,7]{0,1}; rvcyclotomiclocal K[B,1];";

fn setup(invariant: usize, precision: usize) -> String {
    format!(
        "{BASE}
        rvcubicdiv D[K,{invariant},{precision}];
        let c = rvCubicBaseChamber[D];
        let t0 = rvCubicIwT[D,[0]];
        let t1 = rvCubicIwT[D,[1]];
        let id = rvCubicIwT[D,[]];
        let t01 = rvCubicIwT[D,[0,1]];
        let t10 = rvCubicIwT[D,[1,0]];
        let om = rvCubicIwOmega[D,1];
        let delta = rvCubicIwDelta(c);"
    )
}

#[test]
fn runnable_example_composes_actual_galleries_and_transforms_spherical_products() {
    let output = stdout(example("gl2-division-iwahori-hecke"));
    assert_eq!(
        output,
        concat!(
            "rvCubicGallery = length=2, left_descent=0, right_descent=1\n",
            "rvCubicIwKernel = 1\n",
            "rvCubicIwVectorCoeff = 1\n",
            "RVCubicSatake(D){terms=[(-2, \"729\"), (0, \"26\"), (2, \"1\")]}\n",
            "rvCubicSatakeCounts = radius=2, heights=[(-2, 729), (0, 26), (2, 1)]\n",
        )
    );
}

#[test]
fn panels_have_exactly_27_neighbors_and_affine_relations_act_on_real_chambers() {
    for t in [1, 2] {
        let mut program = setup(t, 5);
        for panel in 0..=1 {
            writeln!(
                program,
                "let square{panel} = rvCubicHeckeMul(t{panel},t{panel});
                 assert square{panel} == add(rvCubicHeckeScale[26,1](t{panel}),rvCubicHeckeScale[27,1](id));
                 assert rvCubicIwAct(square{panel},delta) == rvCubicIwAct(t{panel},rvCubicIwAct(t{panel},delta));"
            ).unwrap();
            for j in 0..27 {
                writeln!(program, "show rvCubicChamberPanel[{panel},{j}](c);").unwrap();
            }
        }
        program.push_str(
            "
             let u = rvCubicChamberPanel[0,0](c);
             let v = rvCubicChamberPanel[1,0](u);
             let w = rvCubicChamberPanel[0,0](v);
             let t010 = rvCubicIwT[D,[0,1,0]];
             let t101 = rvCubicIwT[D,[1,0,1]];
             let z = rvCubicIwZ[D];
             assert z == add(add(t01,t10),rvCubicHeckeScale[-26,1](add(t0,t1)));
             assert rvCubicHeckeMul(z,t0) == rvCubicHeckeMul(t0,z);
             assert rvCubicHeckeMul(z,t1) == rvCubicHeckeMul(t1,z);
             assert rvCubicHeckeMul(z,om) == rvCubicHeckeMul(om,z);
             assert rvCubicHeckeMul(t01,t10) ==
               add(rvCubicHeckeScale[26,1](t010),
                   add(rvCubicHeckeScale[702,1](t0),rvCubicHeckeScale[729,1](id)));
             assert add(t0,neg(t0)) == zero[RVCubicIw(D)];
             rvCubicIwWordData[[1,0,1]];
             rvCubicGallery(c,w);
             rvCubicIwKernel(t010,c,w);
             rvCubicIwKernel(t101,c,w);
             rvCubicIwKernel(rvCubicHeckeMul(t0,t0),c,c);
             rvCubicIwKernel(rvCubicHeckeMul(t0,t0),c,u);
             rvCubicIwKernel(rvCubicIwCompress(rvCubicSph[D,2]),c,v);
             show rvCubicIwAct(t01,delta);
             rvCubicIwVectorCoeff(rvCubicIwAct(t01,delta),v);",
        );
        let output = stdout(interpret(&program));
        let rows: Vec<_> = output.lines().collect();
        assert_eq!(rows[..54].iter().copied().collect::<HashSet<_>>().len(), 54);
        assert_eq!(
            rows[54..],
            [
                "rvCubicIwWordData = length=3, left_descent=1, right_descent=1",
                "rvCubicGallery = length=3, left_descent=0, right_descent=0",
                "rvCubicIwKernel = 1",
                "rvCubicIwKernel = 0",
                "rvCubicIwKernel = 27",
                "rvCubicIwKernel = 26",
                "rvCubicIwKernel = 1/28",
                "RVCubicChamberVec(D){support=729}",
                "rvCubicIwVectorCoeff = 1",
            ]
        );
    }
}

#[test]
fn geometric_omega_has_a_nontrivial_square_and_a_trivial_sixth_power() {
    for t in [1, 2] {
        let mut program = setup(t, 3);
        program.push_str(
            "
            let u = rvCubicChamberPanel[0,0](c);
            let v = rvCubicChamberPanel[1,0](u);
            let om2 = rvCubicHeckeMul(om,om);
            let om6 = rvCubicHeckeMul(rvCubicHeckeMul(om2,om2),om2);
            let psi = rvCubicIwOmega[D,-1];
            assert rvCubicHeckeMul(om,psi) == id;
            assert rvCubicHeckeMul(psi,om) == id;
            assert rvCubicHeckeMul(om,t0) == rvCubicHeckeMul(t1,om);
            assert rvCubicOmega[1](c) == c;
            assert rvCubicOmega[-1](rvCubicOmega[1](v)) == v;
            assert rvCubicOmega[6](v) == v;
            assert rvCubicIwAct(om,rvCubicIwDelta(u)) == rvCubicIwDelta(rvCubicOmega[1](u));
            assert rvCubicIwAct(om6,rvCubicIwDelta(v)) == rvCubicIwDelta(v);
            assert rvCubicIwAct(rvCubicHeckeMul(t0,om),delta) ==
               rvCubicIwAct(t1,rvCubicIwAct(om,delta));
            rvCubicIwCoeff[6,[]](om6);
            rvCubicIwCoeff[0,[]](om6);
            rvCubicGallery(c,rvCubicOmega[1](u));",
        );
        for j in 0..28 {
            writeln!(
                program,
                "let x{j}=rvCubicTreeAt[D,1,{j}]; show x{j};
                 show rvCubicOmega[2](x{j});
                 assert rvCubicOmega[2](x{j}) == rvCubicTreeApply(rvCubicTreePi[D],x{j});"
            )
            .unwrap();
        }
        let output = stdout(interpret(&program));
        let lines: Vec<_> = output.lines().collect();
        assert_eq!(
            lines[..3],
            [
                "rvCubicIwCoeff = 1",
                "rvCubicIwCoeff = 0",
                "rvCubicGallery = length=1, left_descent=1, right_descent=1",
            ]
        );
        assert_eq!(
            lines[3..]
                .chunks_exact(2)
                .filter(|xs| xs[0] == xs[1])
                .count(),
            4
        );
    }
}

#[test]
fn spherical_satake_and_affine_compression_are_multiplicative() {
    let mut program = setup(1, 4);
    program.push_str("
            let a0 = rvCubicSph[D,0]; let a1 = rvCubicSph[D,1];
            let a2 = rvCubicSph[D,2]; let a3 = rvCubicSph[D,3];
            let a4 = rvCubicSph[D,4];
            let e = rvCubicIwE[D]; let r2 = rvCubicIwCompress(a2);
            assert e == rvCubicHeckeScale[1,28](add(id,t0));
            assert rvCubicHeckeMul(e,e) == e;
            assert rvCubicHeckeMul(e,r2) == r2;
            assert rvCubicHeckeMul(r2,e) == r2;
            assert rvCubicHeckeMul(a1,a1) == add(a2,rvCubicHeckeScale[28,1](a0));
            assert rvCubicHeckeMul(a1,a2) == add(a3,rvCubicHeckeScale[27,1](a1));
            let p = add(a2,rvCubicHeckeScale[2,3](a0));
            let q = add(a4,neg(a2));
            let pq = rvCubicHeckeMul(p,q);
            assert pq == rvCubicHeckeMul(q,p);
            assert rvCubicIwCompress(pq) == rvCubicHeckeMul(rvCubicIwCompress(p),rvCubicIwCompress(q));
            assert rvCubicSatake(pq) == rvCubicHeckeMul(rvCubicSatake(p),rvCubicSatake(q));
            let s1 = rvCubicSatake(a1); let s3 = rvCubicSatake(a3);
            rvCubicSphCoeff[0](rvCubicHeckeMul(a1,a1));
            rvCubicSatakeCoeff[-1](s1); rvCubicSatakeCoeff[1](s1);
            rvCubicSatakeCoeff[-3](s3); rvCubicSatakeCoeff[-1](s3);
            rvCubicSatakeCoeff[1](s3); rvCubicSatakeCoeff[3](s3);
            let ray = rvCubicTreeAt[D,3,19683];
            rvCubicBusemann(rvCubicTreeAncestor[1](ray),ray);
            rvCubicBusemann(rvCubicTreeAt[D,3,0],ray);
            rvCubicSatakeCounts(ray);");
    let output = stdout(interpret(&program));
    assert_eq!(
        output,
        concat!(
            "rvCubicSphCoeff = 28\n",
            "rvCubicSatakeCoeff = 27\n",
            "rvCubicSatakeCoeff = 1\n",
            "rvCubicSatakeCoeff = 19683\n",
            "rvCubicSatakeCoeff = 702\n",
            "rvCubicSatakeCoeff = 26\n",
            "rvCubicSatakeCoeff = 1\n",
            "rvCubicBusemann = 1\n",
            "rvCubicBusemann = -3\n",
            "rvCubicSatakeCounts = radius=3, heights=[(-3, 19683), (-1, 702), (1, 26), (3, 1)]\n",
        )
    );
}

#[test]
fn nominal_types_precision_and_finite_vector_bounds_are_enforced() {
    let fixtures = setup(2, 2);
    for (statement, expected) in [
        ("show rvCubicIwT[D,[0,0]];", "must alternate"),
        ("show rvCubicChamber(rvCubicTreeRoot[D],rvCubicTreeAt[D,2,0]);", "adjacent lattice vertices"),
        ("show rvCubicChamberVertex[2](c);", "panel must be 0 or 1"),
        ("show rvCubicChamberPanel[0,27](c);", "panel index must be in 0..27"),
        ("show rvCubicIwCompress(rvCubicSph[D,3]);", "even spherical support"),
        ("show rvCubicHeckeScale[1,0](id);", "denominator must be positive"),
        ("rvcubicdiv E[K,2,2]; show rvCubicHeckeMul(id,rvCubicIwT[E,[]]);", "RVCubicIw(D)"),
        ("rvcubicdiv E[K,2,2]; show rvCubicIwAct(t0,rvCubicIwDelta(rvCubicBaseChamber[E]));", "RVCubicChamberVec(D)"),
        ("show rvCubicHeckeMul(delta,delta);", "an affine, spherical or Satake series"),
        ("show rvCubicOmega[1](rvCubicTreeAt[D,2,729]);", "one more outward digit"),
        ("show rvCubicChamberPanel[0,0](rvCubicChamber(rvCubicTreeAt[D,2,0],rvCubicTreeAt[D,1,0]));", "full chamber panel needs another outward digit"),
        ("rvCubicBusemann(rvCubicTreeAt[D,2,0],rvCubicTreeAt[D,1,0]);", "at least the target vertex radius"),
    ] {
        let error = diagnostic(interpret(&format!("{fixtures} {statement}")));
        assert!(error.contains(expected), "{statement}\n{error}");
    }
    stdout(interpret(&format!(
        "{fixtures}
        let last = rvCubicTreeAt[D,2,729];
        assert rvCubicOmega[6](last) == last;
        assert rvCubicOmega[2](last) == rvCubicTreeApply(rvCubicTreePi[D],last);"
    )));
    let large = setup(2, 5);
    let error = diagnostic(interpret(&format!(
        "{large} show rvCubicIwAct(rvCubicIwT[D,[0,1,0]],delta);"
    )));
    assert!(
        error.contains("explicit cubic Iwahori vector exceeds 8192"),
        "{error}"
    );
}
