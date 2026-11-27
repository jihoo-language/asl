mod common;

use common::{diagnostic, example, interpret, stdout};

const MARKED: &str = r"
    rvunram F[3,1]{0,1};
    rvunram W[3,1]{1,0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
    rvjacbase L[C,W];
    let zeroW=rvUnramElem[W]{};
    let oneW=rvUnramElem[W]{1};
    let alpha=rvUnramElem[W]{0,1};
    let e1=rvJacPoint[L](zeroW,zeroW);
    let f1=rvJacPoint[L](oneW,zeroW);
    let e2=add(add(e1,f1),rvJacPoint[L](neg(oneW),zeroW));
    let f2=add(add(e1,f1),rvJacPoint[L](alpha,zeroW));
    rvjactorsionbasis B[L,2](neg(oneW)){e1,e2,f1,f2};
";

const PLANE: &str = r"
    let u=rvK3Lattice[B]{[0,2,0,0,2,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]};
    let v=rvK3Lattice[B]{[3,0,0,0,0,-4],[0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1]};
";

#[test]
fn nonsquare_period_uses_the_integral_glued_lattice_and_both_half_spinors() {
    let output = stdout(example("hyperelliptic-k3-kuga-satake-period"));
    assert_eq!(
        output,
        concat!(
            "rvK3PeriodReport[P] = discriminant=40, saturation_index=1, primitive_discriminant=40, ",
            "NS_rank=20, NS_signature=(1, 19), NS_determinant=-40, half_spin_dimension=1024, ",
            "half_spin_CM_rank=512, even_H1_dimension=2097152\n",
            "rvK3PeriodReport[Scaled] = discriminant=1440, saturation_index=6, primitive_discriminant=40, ",
            "NS_rank=20, NS_signature=(1, 19), NS_determinant=-40, half_spin_dimension=1024, ",
            "half_spin_CM_rank=512, even_H1_dimension=2097152\n",
            "rvK3PeriodIs11[P] = true\n",
            "rvK3PeriodIs11[P] = false\n",
            "rvK3PeriodIs11[P] = false\n",
        ),
    );
}

#[test]
fn sheared_nonprimitive_plane_keeps_the_primitive_discriminant() {
    let output = stdout(interpret(&format!(
        "{MARKED}{PLANE}
         rvk3period Original(u,v);
         rvk3period Skew(add(u,v),add(u,add(u,add(u,add(u,v)))));
         assert rvK3PeriodBlade[Skew]==rvCliffordScale[-3,1](rvK3PeriodBlade[Original]);
         rvK3PeriodReport[Skew];"
    )));
    assert!(
        output.starts_with(
            "rvK3PeriodReport[Skew] = discriminant=360, saturation_index=3, primitive_discriminant=40"
        ),
        "{output}",
    );
}

#[test]
fn report_supports_wide_discriminants_and_indices() {
    let large = i128::MAX - 1;
    let u = glued_vector([0, large, 0, 0, large, 0]);
    let v = glued_vector([large, 0, 0, 0, 0, -large]);
    let output = stdout(interpret(&format!(
        "{MARKED} rvk3period Large({u},{v}); rvK3PeriodReport[Large];"
    )));
    // In the pure wedge plane these are multiples of two primitive,
    // orthogonal square-four vectors by (2^127−2)/2.
    let discriminant = "837987995621412318723376562387865382928058357590630106895443550132466638075211699315075594140024596736433396812272342957068459789683076491174045559029776";
    let index = "7237005577332262213973186563042994240659232858142066020734411696778686496769";
    assert!(
        output.starts_with(&format!(
            "rvK3PeriodReport[Large] = discriminant={discriminant}, saturation_index={index}, primitive_discriminant=16"
        )),
        "{output}",
    );
}

#[test]
fn out_of_range_blades_are_diagnosed_without_restricting_period_reports() {
    let max = i128::MAX;
    let prev = max - 1;
    let u = glued_vector([max, max, 1, 0, prev, -prev]);
    let v = glued_vector([max, -max, 0, 0, -prev, -prev]);
    let setup = format!("{MARKED} rvk3period Large({u},{v});");
    let output = stdout(interpret(&format!("{setup} rvK3PeriodReport[Large];")));
    assert!(
        output.starts_with("rvK3PeriodReport[Large] = discriminant="),
        "{output}"
    );
    let error = diagnostic(interpret(&format!("{setup} show rvK3PeriodBlade[Large];")));
    assert!(
        error.contains("Clifford rational coefficient exceeds the supported exact integer range"),
        "{error}"
    );
}

#[test]
fn rejects_degenerate_indefinite_and_nominally_mixed_planes() {
    let setup = format!(
        "{MARKED}{PLANE}
         rvjactorsionbasis Twin[L,2](neg(oneW)){{e1,e2,f1,f2}};
         rvk3period P(u,v);
         let root=rvK3FromKummerDiv(rvKummerNodeDiv[B,0]);"
    );
    for (command, expected) in [
        (
            "rvk3period Invalid(root,v);",
            "first lattice vector to have positive square",
        ),
        (
            "rvk3period Invalid(u,u);",
            "requires a positive Gram determinant",
        ),
        (
            "rvk3period Invalid(u,root);",
            "requires a positive Gram determinant",
        ),
        (
            "rvk3period Invalid(u,rvK3LatticeBasis[Twin,0]);",
            "expected type RVK3Lattice(B)",
        ),
        ("rvk3period P(u,v);", "K3 period 'P' is already defined"),
        ("rvK3PeriodReport[Missing];", "unknown K3 period 'Missing'"),
        (
            "rvK3PeriodIs11[P](rvK3LatticeBasis[Twin,0]);",
            "expected type RVK3Lattice(B)",
        ),
        (
            "show rvK3PeriodCM[P](rvK3SpinorBasis[Twin,1,0]);",
            "expected type RVK3Spinor(B,1)",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}{command}")));
        assert!(error.contains(expected), "{error}");
    }
}

fn glued_vector(w: [i128; 6]) -> String {
    let pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
    let b: [i128; 16] = std::array::from_fn(|node| {
        w.iter()
            .zip(pairs)
            .map(|(entry, (i, j))| (entry & 1) * (((node >> i) & (node >> j) & 1) as i128))
            .sum::<i128>()
            & 1
    });
    format!("rvK3Lattice[B]{{{w:?},{b:?}}}")
}
