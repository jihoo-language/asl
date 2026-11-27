mod common;

use common::{diagnostic, example, interpret, stdout};

fn fixtures(t: usize, precision: usize) -> String {
    format!(
        "rvunram B[3,7]{{0,1}}; rvcyclotomiclocal K[B,1];
         rvcubicdiv D[K,{t},{precision}];
         let o = rvCubicTreeRoot[D];
         let l = rvCubicTreeAt[D,1,0];
         let r = rvCubicTreeAt[D,1,27];
         let e1 = rvCubicBoundaryCylinder(l);
         let e2 = rvCubicBoundaryCylinder(r);
         let one = rvCubicBoundaryConstant[D,1,1];
         let th = rvCubicDivTheta[D];
         let pi = rvCubicDivPi[D];
         let unit = rvCubicDivElem[D]([1]);
         let g = rvCubicTreeMatrix(unit,th,pi,unit);
         let h = rvCubicTreePi[D];"
    )
}

#[test]
fn sample_computes_actual_poisson_transform_and_the_three_knapp_stein_scalars() {
    let output = stdout(example("gl2-principal-series-knapp-stein"));
    assert_eq!(
        output,
        concat!(
            "RVCubicBoundary(D){level=3, minimal_level=1, nonzero_cylinders=2}\n",
            "rvCubicBoundaryMass = 1/28\n",
            "rvCubicBoundaryValue = 1\n",
            "rvCubicPoisson = 2/21\n",
            "rvCubicPoissonEigen = eigenvalue=12, value=2/21, adjacency=8/7\n",
            "rvCubicKnappSpectrum = level=1, eigenvalue=10/3\n",
            "rvCubicKnappSpectrum = level=2, eigenvalue=10\n",
            "rvCubicKnappSpectrum = level=3, eigenvalue=30\n",
            "rvCubicKnappKernel = 13/4\n",
            "rvCubicPoissonCovariance = z=9, value=2/21, inverse_radius=1\n",
        )
    );
}

#[test]
fn cylinder_refinement_integrals_coarsening_and_martingale_decomposition() {
    let output = stdout(interpret(&format!(
        "{}
         let x = rvCubicTreeChild[1](l);
         let y = rvCubicTreeChild[4](x);
         let a = rvCubicBoundaryCylinder(y);
         let e = rvCubicBoundaryRefine[5](a);
         assert e == a;
         assert rvCubicBoundaryCoarsen[3](e) == a;
         let u0 = rvCubicBoundaryDifference[0](a);
         let u1 = rvCubicBoundaryDifference[1](a);
         let u2 = rvCubicBoundaryDifference[2](a);
         let u3 = rvCubicBoundaryDifference[3](a);
         assert a == add(add(u0,u1),add(u2,u3));
         assert rvCubicBoundaryDifference[4](a) == zero[RVCubicBoundary(D)];
         assert rvCubicBoundaryAverage[0](u3) == zero[RVCubicBoundary(D)];
         assert rvCubicBoundaryAverage[2](u3) == zero[RVCubicBoundary(D)];
         assert rvCubicBoundaryAverage[2](a) ==
            rvCubicBoundaryScale[1,27](rvCubicBoundaryCylinder(x));
         assert rvCubicBoundaryAverage[1](a) == rvCubicBoundaryScale[1,729](e1);
         assert rvCubicBoundaryAverage[0](a) == rvCubicBoundaryConstant[D,1,20412];
         assert rvCubicKnappStein[3,1](a) ==
            add(add(u0,rvCubicBoundaryScale[10,3](u1)),
                add(rvCubicBoundaryScale[10,1](u2),rvCubicBoundaryScale[30,1](u3)));
         assert rvCubicKnappStein[9,1](rvCubicKnappStein[3,1](a)) == a;
         assert rvCubicKnappStein[27,1](a) == u0;
         assert rvCubicKnappStein[-27,1](a) == u0;
         rvCubicBoundaryMass(y);
         rvCubicBoundaryIntegral(a,l);
         rvCubicBoundaryIntegral(a,r);
         rvCubicBoundaryValue(a,rvCubicTreeChild[7](y));
         rvCubicKnappKernel[3,1](l,r);
         rvCubicKnappKernel[9,1](l,r);
         rvCubicKnappKernel[27,1](y,o);
         rvCubicKnappSpectrum[9,1,1];
         rvCubicKnappSpectrum[9,1,2];
         rvCubicKnappSpectrum[9,1,3];",
        fixtures(2, 5)
    )));
    assert_eq!(
        output,
        concat!(
            "rvCubicBoundaryMass = 1/20412\n",
            "rvCubicBoundaryIntegral = 1/20412\n",
            "rvCubicBoundaryIntegral = 0\n",
            "rvCubicBoundaryValue = 1\n",
            "rvCubicKnappKernel = -1/12\n",
            "rvCubicKnappKernel = 1/40\n",
            "rvCubicKnappKernel = 1/20412\n",
            "rvCubicKnappSpectrum = level=1, eigenvalue=3/10\n",
            "rvCubicKnappSpectrum = level=2, eigenvalue=1/10\n",
            "rvCubicKnappSpectrum = level=3, eigenvalue=1/30\n",
        )
    );
}

#[test]
fn actual_root_change_transports_complements_and_lambda_digits() {
    for t in [1, 2] {
        let output = stdout(interpret(&format!(
            "{}
             let x = rvCubicTreeAt[D,3,1];
             let y = rvCubicTreeChild[9](x);
             let f = rvCubicBoundaryCylinder(y);
             let image = rvCubicTreeTranslate(y);
             assert rvCubicPrincipalTranslate[3,1](f) ==
                rvCubicBoundaryScale[3,1](rvCubicBoundaryCylinder(image));
             assert rvCubicPrincipalTranslateBack[3,1](rvCubicPrincipalTranslate[3,1](f)) == f;
             assert rvCubicPrincipal[3,1](g,f) ==
                rvCubicBoundaryCylinder(rvCubicTreeApply(g,y));
             assert rvCubicPrincipal[3,1](rvCubicTreeCompose(g,h),f) ==
                rvCubicPrincipal[3,1](g,rvCubicPrincipal[3,1](h,f));
             let outward = add(one,neg(e2));
             assert rvCubicPrincipalTranslate[3,1](e1) == rvCubicBoundaryScale[1,3](outward);
             assert rvCubicPrincipalTranslate[-3,1](e1) == rvCubicBoundaryScale[-1,3](outward);
             assert rvCubicPrincipalTranslate[3,1](one) ==
                add(rvCubicBoundaryScale[3,1](e2),rvCubicBoundaryScale[1,3](outward));
             let omega2 = rvCubicPrincipalOmega[3,1,2](f);
             assert omega2 == rvCubicPrincipal[3,1](h,f);
             assert rvCubicPrincipalOmega[3,1,6](f) == f;
             rvCubicPoissonCovariance[3,1](g,f,y);
             rvCubicPoissonCovariantTranslate[3,1](e1,o);
             rvCubicPoissonCovariantBack[3,1](e2,o);
             rvCubicPoissonCovariantOmega[-3,1,1](e1,o);",
            fixtures(t, 6)
        )));
        assert!(output.contains("z=9, value=3/28, inverse_radius=1"));
        assert!(output.contains("z=-9, value=-3/28, inverse_radius=1"));
        assert_eq!(output.lines().count(), 4);
    }
}

#[test]
fn normalized_intertwiner_crosses_both_root_moves_and_the_actual_compact_action() {
    for t in [1, 2] {
        let output = stdout(interpret(&format!(
            "{}
             let w = rvCubicTreeChild[1](l);
             let f = add(rvCubicBoundaryScale[2,3](rvCubicBoundaryCylinder(w)),
                         add(rvCubicBoundaryScale[-5,2](e2),rvCubicBoundaryScale[1,7](one)));
             let a = rvCubicKnappStein[5,2](f);
             assert rvCubicKnappStein[54,5](a) == f;
             assert rvCubicKnappStein[5,2](rvCubicPrincipal[5,2](g,f)) ==
                rvCubicPrincipal[54,5](g,a);
             assert rvCubicKnappStein[5,2](rvCubicPrincipalTranslate[5,2](f)) ==
                rvCubicPrincipalTranslate[54,5](a);
             assert rvCubicKnappStein[5,2](rvCubicPrincipalOmega[5,2,-1](f)) ==
                rvCubicPrincipalOmega[54,5,-1](a);
             let b = rvCubicKnappStein[-3,1](f);
             assert rvCubicKnappStein[-9,1](b) == f;
             assert rvCubicKnappStein[-3,1](rvCubicPrincipalTranslate[-3,1](f)) ==
                rvCubicPrincipalTranslate[-9,1](b);
             rvCubicPoissonEigen[5,2](f,w);
             rvCubicPoissonCovariantTranslate[54,5](f,w);
             rvCubicPoissonCovariantOmega[54,5,-1](f,w);",
            fixtures(t, 5)
        )));
        assert!(output.starts_with("rvCubicPoissonEigen = eigenvalue=133/10"));
        assert_eq!(output.matches("z=5/2").count(), 2);
    }
}

#[test]
fn invalid_parameters_unsupported_geometry_and_type_mismatches_are_diagnostics() {
    let base = fixtures(1, 2);
    for (input, message) in [
        (
            "show rvCubicBoundaryRefine[0](e1);",
            "refinement requires a level at least",
        ),
        ("show rvCubicBoundaryRefine[3](e1);", "radius exceeds"),
        (
            "show rvCubicBoundaryCoarsen[0](e1);",
            "function is not constant",
        ),
        (
            "rvCubicBoundaryValue(e1,o);",
            "varies on the requested cylinder",
        ),
        (
            "rvCubicKnappKernel[3,1](l,o);",
            "varies on the requested cylinder",
        ),
        (
            "rvCubicPoisson[0,1](e1,o);",
            "spectral parameter must be nonzero",
        ),
        ("show rvCubicKnappStein[-1,1](e1);", "has a pole"),
        ("rvCubicKnappSpectrum[1,1,3];", "has a pole"),
        (
            "show rvCubicBoundaryConstant[D,1,-1];",
            "denominator must be positive",
        ),
        (
            "rvCubicPoissonEigen[3,1](e1,rvCubicTreeAt[D,2,0]);",
            "needs another outward digit",
        ),
        (
            "show rvCubicPrincipalTranslate[3,1](rvCubicBoundaryCylinder(rvCubicTreeAt[D,2,729]));",
            "one more outward digit",
        ),
        (
            "rvcubicdiv E[K,1,2]; rvCubicBoundaryValue(e1,rvCubicTreeRoot[E]);",
            "RVCubicTreeVertex(D)",
        ),
        (
            "rvcubicdiv E[K,1,2]; show rvCubicPrincipal[3,1](g,rvCubicBoundaryConstant[E,1,1]);",
            "RVCubicBoundary(D)",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{base} {input}")));
        assert!(error.contains(message), "{input}\n{error}");
    }
}
