mod common;

use common::{diagnostic, example, interpret, stdout};

const F41: &str = include_str!("../examples/elliptic-cm-hermitian-period.asl");

#[test]
fn nonprincipal_f41_period_has_exact_determinant_fourier_cocycle_and_wall() {
    let output = stdout(example("elliptic-cm-hermitian-period"));
    for expected in [
        "determinant=[0, 16, 1], imaginary_determinant=6, factor=[1, 0, 1]",
        "determinant=[0, -1, 16], imaginary_determinant=3/128, factor=[0, 16, 1]",
        "kind=quadratic, center=-5/2, radius_squared=9/4, positive_height_locus=true, coefficients=[1, 5, 4]",
        "charge_quadratics=[[-1, 0, 4], [-1, -2, -1]]",
        "rvCMPeriodWallAt = on_locus=true, orientation=same_ray, first_zero=false, second_zero=false, ratio=5/3, charges=[[0, 15, 2], [0, 9, 2]]",
        "rvCMChargeRay = orientation=same_ray, first_zero=false, second_zero=false, ratio=5/3",
    ] {
        assert!(output.contains(expected), "{output}");
    }
}

#[test]
fn linear_collinear_empty_and_zero_charge_loci_are_distinguished() {
    let output = stdout(interpret(&format!(
        "{F41}
         let single=rvCMAbelianMukai[0,-1](z0);
         let linear=rvCMPeriodWall(I,z0,a,single);
         rvCMPeriodWallData(linear);
         rvCMPeriodWallAt[0,1,1,1](linear);
         let negative=rvCMPeriodWall(I,z0,a,neg(a));
         rvCMPeriodWallData(negative);
         rvCMPeriodWallAt[0,1,1,1](negative);
         let zero_class=rvCMAbelianMukai[1,4](z0);
         let zeros=rvCMPeriodWall(I,z0,zero_class,b);
         rvCMPeriodWallAt[0,1,2,1](zeros);
         let line_charge=rvCMAbelianMukai[0,0](I);
         let none=rvCMPeriodWall(I,z0,line_charge,single);
         rvCMPeriodWallData(none);
         rvCMPeriodWallAt[0,1,1,1](none);
         let no_upper_height=rvCMPeriodWall(I,z0,rvCMAbelianMukai[1,0](z0),line_charge);
         rvCMPeriodWallData(no_upper_height);
         rvCMPeriodWallAt[0,1,1,1](no_upper_height);
         rvCMChargeRay(zero[RVCMCharge(A)],zero[RVCMCharge(A)]);
         let wide=rvCMCharge[A,1267650600228229401496703205376,1267650600228229401496703205376,1];
         rvCMChargeRay(wide,wide);"
    )));
    for expected in [
        "kind=linear, real_coordinate=0, coefficients=[0, 1, 0]",
        "on_locus=true, orientation=same_ray, first_zero=false, second_zero=false, ratio=5",
        "kind=identically_collinear, coefficients=[0, 0, 0]",
        "on_locus=true, orientation=opposite_ray, first_zero=false, second_zero=false, ratio=-1",
        "on_locus=true, orientation=zero, first_zero=true, second_zero=false, ratio=undefined",
        "kind=empty, coefficients=[0, 0, 1]",
        "kind=quadratic, center=0, radius_squared=0, positive_height_locus=false, coefficients=[1, 0, 0]",
        "on_locus=false, orientation=different_rays",
        "rvCMChargeRay = orientation=zero, first_zero=true, second_zero=true",
        "rvCMChargeRay = orientation=same_ray, first_zero=false, second_zero=false, ratio=1",
    ] {
        assert!(output.contains(expected), "{output}");
    }
}

#[test]
fn generic_hermitian_slice_with_nonzero_offset_matches_actual_period_charge() {
    let output = stdout(interpret(&format!(
        "{F41}
         let translated=rvCMPeriodWall(D,J,v,u);
         rvCMPeriodWallData(translated);
         rvCMPeriodWallAt[0,1,1,1](translated);
         assert rvCMPeriodCharge(P,u)==rvCMCharge[A,7,-32,1];"
    )));
    assert!(output
        .contains("coefficients=[6, 111, 95], charge_quadratics=[[-12, -20, 5], [-18, -32, -11]]"));
    assert!(output.contains("on_locus=false, orientation=different_rays"));
    assert!(output.contains("charges=[[17, -20, 1], [7, -32, 1]]"));
}

#[test]
fn central_imaginary_unit_remains_independent_in_the_gaussian_cm_field() {
    let output = stdout(interpret(
        r"
        rvunram k[5,1]{0,1};
        rvec E[k](rvUnramElem[k]{2},rvUnramElem[k]{0});
        rvec F[k](rvUnramElem[k]{2},rvUnramElem[k]{0});
        let f=rvECChart[E,F](rvUnramElem[k]{1});
        let g=rvECChart[E,F](rvUnramElem[k]{2});
        let fd=rvECChart[F,E](rvUnramElem[k]{1});
        let gd=rvECChart[F,E](rvUnramElem[k]{3});
        rvcmhomframe H(f,g);
        rvcmhomdual HD[H](fd,gd);
        rvcmhomframe RE(rvECMultMap[E,1],rvECMapCompose(fd,g));
        rvcmhomframe RF(rvECMultMap[F,1],rvECMapCompose(g,fd));
        rvcmabelian A[H,HD,RE,RF];
        let beta=rvCMAbelianDiv[A,1,0](rvCMHomCoord[H]{0,1});
        let l=rvCMAbelianDiv[A,2,2](rvCMHomCoord[H]{1,0});
        let p=rvCMPeriod(beta,l);
        let s=rvCMPeriodFourier(p);
        let v=rvCMAbelianMukai[1,1](beta);
        assert rvCMPeriodDet(p)==rvCMCharge[A,-4,2,1];
        assert rvCMPeriodCharge(p,v)==rvCMCharge[A,1,0,1];
        assert rvCMChargeMul(rvCMPeriodFactor(s),rvCMPeriodCharge(s,rvCMAbelianFourier(v)))==rvCMCharge[A,1,0,1];
        assert rvCMPeriodFourier(s)==p;
        rvCMPeriodData(s);",
    ));
    assert!(
        output.contains("imaginary_determinant=3/20, factor=[-4, 2, 1]"),
        "{output}"
    );
}

#[test]
fn positivity_and_nominal_types_are_checked_in_every_entry_point() {
    for (query, expected) in [
        ("show rvCMPeriod(I,neg(I));", "must be positive definite"),
        ("show rvCMPeriod(I,rvCMAbelianDiv[A,0,1](rvCMHomCoord[H]{1,0}));", "must be positive definite"),
        ("show rvCMPeriodWall(z0,z0,a,b);", "must be positive definite"),
        ("rvCMPeriodWallAt[0,1,0,1](wall);", "height must be positive"),
        ("rvCMPeriodWallAt[0,1,1,0](wall);", "denominator must be positive"),
        ("show rvCMCharge[A,0,1,-1];", "denominator must be positive"),
        ("let m=rvCMAbelianMap[A](id,rvCMHomMap[HD,0,0],rvCMHomMap[H,0,0],rvECMultMap[F,0]); show rvCMPeriodPull(P,m);", "requires a nonsingular product isogeny"),
        ("rvcmabelian Other[H,HD,R,RF]; show rvCMPeriod(I,zero[RVCMAbelianDiv(Other)]);", "CM period imaginary divisor"),
        ("rvcmabelian Other[H,HD,R,RF]; rvCMChargeRay(rvCMCharge[Other,0,1,1],c);", "CM charge ray arguments"),
        ("rvCMPeriodWallData(c);", "requires a wall"),
    ] {
        let error = diagnostic(interpret(&format!("{F41}\n{query}")));
        assert!(error.contains(expected), "{error}");
    }
}
