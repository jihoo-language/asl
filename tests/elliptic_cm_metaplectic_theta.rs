mod common;

use common::{diagnostic, example, interpret, stdout};

const SAMPLE: &str = include_str!("../examples/elliptic-cm-metaplectic-theta.asl");

#[test]
fn genuine_f41_polarizations_give_exact_fourier_chirp_and_mixed_shear() {
    let output = stdout(example("elliptic-cm-metaplectic-theta"));
    for expected in [
        "dimension=6, conductor=12, gauss=[-1, 2, 2, -1]/1, triple_trace=[-1, 2, 2, -1]/36",
        "top_chirp_power_phases=[0, 6, 0, 6, 0, 6], top_chirp_power_clock=[0, 3]",
        "dimension=8, conductor=8, gauss=[0, 2, 0, 2]/1, triple_trace=[0, 1, 0, 1]/32",
        "top_chirp_power_phases=[0, 4, 0, 4, 0, 4, 0, 4], top_chirp_power_clock=[0, 2]",
        "RVCycloQ(8){[8] / 1}",
    ] {
        assert!(output.contains(expected), "{output}");
    }
}

#[test]
fn odd_and_even_cyclic_dimensions_and_double_polarization_match_exact_gauss_sums() {
    let output = stdout(interpret(&format!(
        "{SAMPLE}
        rvcmthetakernel K1(rvCMAbelianDiv[A,1,1](off));
        rvcmthetakernel K2(rvCMAbelianDiv[A,1,2](off));
        rvcmthetakernel K3(rvCMAbelianDiv[A,1,3](off));
        rvcmthetakernel K4(rvCMAbelianDiv[A,1,4](off));
        rvcmthetakernel K5(rvCMAbelianDiv[A,1,5](off));
        rvcmthetakernel K22(rvCMAbelianDiv[A,2,2](off));
        rvCMMetaVerify[K1];
        rvCMMetaVerify[K2];
        rvCMMetaVerify[K3];
        rvCMMetaVerify[K4];
        rvCMMetaVerify[K5];
        rvCMMetaVerify[K22];"
    )));
    for expected in [
        "dimension=1, conductor=2, gauss=[1]/1, triple_trace=[1]/1, top_chirp_power_phases=[0], top_chirp_power_clock=[0, 0]",
        "dimension=2, conductor=4, gauss=[1, 1]/1, triple_trace=[1, 1]/4, top_chirp_power_phases=[0, 2], top_chirp_power_clock=[0, 1]",
        "dimension=3, conductor=6, gauss=[1, -2]/1, triple_trace=[1, -2]/9, top_chirp_power_phases=[0, 0, 0], top_chirp_power_clock=[0, 0]",
        "dimension=4, conductor=8, gauss=[0, 2]/1, triple_trace=[0, 1]/8, top_chirp_power_phases=[0, 4, 0, 4], top_chirp_power_clock=[0, 2]",
        "dimension=5, conductor=10, gauss=[-1, 0, -2, 2]/1, triple_trace=[-1, 0, -2, 2]/25, top_chirp_power_phases=[0, 0, 0, 0, 0]",
        "dimension=4, conductor=4, gauss=[0, 2]/1, triple_trace=[0, 1]/8, top_chirp_power_phases=[0, 2, 2, 0], top_chirp_power_clock=[1, 1]",
    ] {
        assert!(output.contains(expected), "{output}");
    }
}

#[test]
fn shifted_clock_projector_selects_nonzero_pivot_and_retains_scalar_phase() {
    stdout(interpret(&format!(
        "{SAMPLE}
        let neutral=rvCMMetaLift[0](zero[RVCMTheta(K6)]);
        let x=rvCMThetaElem[K6]{{0,0,1,0}};
        let z=rvCMThetaElem[K6]{{0,0,0,1}};
        let shift=rvCMMetaAuto[K6](neutral,neutral,rvCMMetaLift[0](x),rvCMMetaLift[4](z));
        let phase=rvCMMetaAuto[K6](neutral,neutral,rvCMMetaLift[10](x),rvCMMetaLift[0](z));
        let Xsq=rvCMMetaHeis(rvCMMetaLift[0](add(x,x)));
        assert rvCMMetaIntertwine(shift)==Xsq;
        assert rvCMMetaCompose(rvCMMetaIntertwine(phase),Xsq)==rvCMMetaScale(
            rvCMMetaIntertwine(rvCMMetaAutoCompose(phase,shift)),
            rvCycloRat[1](rvCyclo[12,-4]));
        assert rvCMMetaCocycle(phase,shift)==rvCycloRat[1](rvCyclo[12,-4]);
        assert rvCMMetaCocycle(shift,phase)==one12;
        assert rvCMMetaAutoCompose(phase,shift)==rvCMMetaAutoCompose(shift,phase);
        let a=rvCMThetaElem[K6]{{0,0,2,1}};
        let h=rvCMMetaLift[5](a);
        assert rvCMMetaImage(sa,h)==rvCMMetaLift[9](rvCMThetaElem[K6]{{0,0,-1,2}});
        assert rvCMMetaImage(ta,h)==rvCMMetaLift[9](rvCMThetaElem[K6]{{0,0,2,-1}});"
    )));
}

#[test]
fn arbitrary_exact_vectors_and_twisted_lagrangian_eigenlines() {
    stdout(interpret(&format!(
        "{SAMPLE}
        let zero8=zero[RVCycloQ(8)];
        let v=rvCMMetaVector[K8]{{zero8,z8,zero8,zero8,zero8,n8,zero8,zero8}};
        assert rvCMMetaCoeff[5](v)==n8;
        assert rvCMMetaApply(rvCMMetaFInv[K8],rvCMMetaApply(W,v))==v;
        let twistedZ2=rvCMMetaLift[4](rvCMThetaElem[K8]{{0,0,0,2}});
        let Twisted=rvCMMetaProjector[K8]{{X1,rvCMMetaLiftCompose(X2,X2),twistedZ2}};
        assert rvCMMetaCompose(Twisted,Twisted)==Twisted;
        assert rvCMMetaTrace(Twisted)==z8;
        assert rvCMMetaTrace(rvCMMetaCompose(Twisted,Nonfree))==zero8;
        assert rvCMMetaEntry[1,1](Twisted)==rvCycloRat[4](rvCyclo[8,0]);"
    )));
}

#[test]
fn checks_generator_relations_lifted_isotropy_coefficients_and_nominality() {
    for (query, expected) in [
        (
            "show rvCMMetaAuto[K8](rvCMMetaLift[1](rvCMThetaElem[K8]{1,0,0,0}),Z1,X2,Z2);",
            "does not satisfy its order relation",
        ),
        (
            "show rvCMMetaAuto[K8](X1,X1,X2,Z2);",
            "do not preserve the commutator",
        ),
        (
            "show rvCMMetaProjector[K8]{rvCMMetaLift[1](zero[RVCMTheta(K8)])};",
            "nontrivial central phase",
        ),
        (
            "show rvCMMetaProjector[K8]{X1,X2,Z2};",
            "requires a lifted maximal isotropic subgroup",
        ),
        (
            "show rvCMMetaProjector[K8]{X1};",
            "requires a lifted maximal isotropic subgroup",
        ),
        (
            "show rvCMMetaCompose(W,S);",
            "CM metaplectic second argument",
        ),
        (
            "show rvCMMetaScale(W,one12);",
            "CM metaplectic second argument",
        ),
        ("show rvCMMetaEntry[8,0](W);", "basis index is out of range"),
        (
            "show rvCMMetaVector[K8]{one12};",
            "CM metaplectic vector coefficient",
        ),
        (
            "show rvCMMetaVector[K8]{z8};",
            "needs exactly the dimension many coefficients",
        ),
        (
            "rvcmthetakernel Huge(rvCMAbelianDiv[A,1,33](off)); show rvCMMetaF[Huge];",
            "32-dimensional matrix limit",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{SAMPLE}\n{query}")));
        assert!(error.contains(expected), "{error}");
    }
}
