mod common;

use common::{diagnostic, example, interpret, stdout};
use std::fmt::Write;

const BASIS: &str = r"
    rvunram F[3,1]{0,1};
    rvunram W[3,1]{1,0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
    rvjacbase L[C,W];
    let zw=rvUnramElem[W]{};
    let ow=rvUnramElem[W]{1};
    let i=rvUnramElem[W]{0,1};
    let e1=rvJacPoint[L](zw,zw);
    let f1=rvJacPoint[L](ow,zw);
    let e2=add(add(e1,f1),rvJacPoint[L](neg(ow),zw));
    let f2=add(add(e1,f1),rvJacPoint[L](i,zw));
    rvjactorsionbasis B[L,2](neg(ow)){e1,e2,f1,f2};
";

#[test]
fn all_integral_basis_vectors_round_trip_and_preserve_pairings() {
    let mut source = format!(
        r"{BASIS}
        rvk3autom T[B]{{translate(3)}};
        rvk3autom S[B]{{switch(5)}};
        rvk3autom P[B]{{project(1)}};
        rvk3autom Forward[B]{{translate(3),switch(5),project(1)}};
        rvk3autom Back[B]{{project(1),switch(5),translate(3)}};
        let p=rvK3Lattice[B]{{[2,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}};
        let q=rvK3Lattice[B]{{[0,2,0,0,-2,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}};
        assert rvK3AutoLattice[T](p)==p;
        assert rvK3AutoLattice[S](p)==neg(p);
        assert rvK3AutoLattice[P](q)==neg(q);
        assert rvK3AutoLattice[Forward](add(p,q))==add(rvK3AutoLattice[Forward](p),rvK3AutoLattice[Forward](q));
        let d=add(rvKummerHyperplane[B],add(rvKummerNodeDiv[B,9],rvKummerTropeDiv[B,10]));
        assert rvK3AutoLattice[Forward](rvK3FromKummerDiv(d))==rvK3FromKummerDiv(rvK3AutoDiv[Forward](d));
        assert rvK3ToKummerDiv(rvK3AutoLattice[Back](rvK3FromKummerDiv(d)))==rvK3AutoDiv[Back](d);
        "
    );
    for index in 0..22 {
        let other = (7 * index + 3) % 22;
        writeln!(
            source,
            "let x{index}=rvK3LatticeBasis[B,{index}];
             let y{index}=rvK3LatticeBasis[B,{other}];
             let fx{index}=rvK3AutoLattice[Forward](x{index});
             let fy{index}=rvK3AutoLattice[Forward](y{index});
             assert fx{index}==rvK3AutoLattice[P](rvK3AutoLattice[S](rvK3AutoLattice[T](x{index})));
             assert rvK3AutoLattice[Back](fx{index})==x{index};
             rvK3LatticePair(x{index},y{index});
             rvK3LatticePair(fx{index},fy{index});"
        )
        .unwrap();
        for generator in ["T", "S", "P"] {
            writeln!(
                source,
                "assert rvK3AutoLattice[{generator}](rvK3AutoLattice[{generator}](x{index}))==x{index};"
            )
            .unwrap();
        }
    }
    let output = stdout(interpret(&source));
    let pairings: Vec<_> = output.lines().collect();
    assert_eq!(pairings.len(), 44);
    assert!(pairings.chunks_exact(2).all(|pair| pair[0] == pair[1]));
}

#[test]
fn every_two_torsion_translation_has_the_expected_invariant_lattice() {
    let mut source = BASIS.to_string();
    for mask in 0..16 {
        writeln!(
            source,
            "rvk3autom T{mask}[B]{{translate({mask})}}; rvK3AutoNikulin[T{mask}];"
        )
        .unwrap();
    }
    let output = stdout(interpret(&source));
    let mut reports = output.lines();
    assert_eq!(
        reports.next().unwrap(),
        "rvK3AutoNikulin[T0] = identity, r=22, a=0, delta=0, discriminant=(Z/2)^0, fixed locus=whole K3 surface"
    );
    for mask in 1..16 {
        assert_eq!(
            reports.next().unwrap(),
            format!("rvK3AutoNikulin[T{mask}] = symplectic, r=14, a=8, delta=0, discriminant=(Z/2)^8, fixed locus=8 isolated points")
        );
    }
    assert!(reports.next().is_none());
}

#[test]
fn checked_wide_action_preserves_cancellation_through_an_entire_word() {
    let high = i128::MAX - 1;
    let zeros = [0; 16];
    let mut orthogonal_nodes = zeros;
    orthogonal_nodes[0] = -high;
    assert_eq!(
        stdout(interpret(&format!(
            r"{BASIS}
            rvk3autom P[B]{{project(0)}};
            rvk3autom Twice[B]{{project(0),project(0)}};
            let large=rvK3Lattice[B]{{[0,{high},0,0,{high},0],{zeros:?}}};
            let orthogonal=rvK3Lattice[B]{{[0,{high},0,0,{high},0],{orthogonal_nodes:?}}};
            assert rvK3AutoLattice[Twice](large)==large;
            assert rvK3AutoLattice[P](orthogonal)==orthogonal;
            rvK3AutoNikulin[Twice];
            "
        ))),
        "rvK3AutoNikulin[Twice] = identity, r=22, a=0, delta=0, discriminant=(Z/2)^0, fixed locus=whole K3 surface\n"
    );
}

#[test]
fn rejects_noninvolutions_cross_marking_and_coordinate_overflow() {
    let high = i128::MAX - 1;
    let long = ["project(0),project(1),project(2)"; 103].join(",");
    for (program, message) in [
        ("rvk3autom Bad[B]{project(0),project(1)}; rvK3AutoNikulin[Bad];".into(), "square is the identity"),
        ("rvk3autom Bad[B]{project(0),project(1),project(2)}; rvK3AutoNikulin[Bad];".into(), "square is the identity"),
        ("rvK3AutoNikulin[Absent];".into(), "unknown K3 automorphism word"),
        ("rvjactorsionbasis Other[L,2](neg(ow)){e1,e2,f1,f2}; show rvK3AutoLattice[P](rvK3LatticeBasis[Other,0]);".into(), "expected type RVK3Lattice(B)"),
        (format!("show rvK3AutoLattice[P](rvK3Lattice[B]{{[0,{high},0,0,{high},0],{:?}}});", [0; 16]), "exceeds signed 128-bit range"),
        (format!("rvk3autom Long[B]{{{long}}}; show rvK3AutoLattice[Long](rvK3FromKummerDiv(rvKummerHyperplane[B]));"), "exceeds the exact intermediate range"),
    ] {
        let problem = diagnostic(interpret(&format!("{BASIS} rvk3autom P[B]{{project(0)}}; {program}")));
        assert!(problem.contains(message), "{problem}");
    }
}

#[test]
fn executable_example_reports_four_fixed_locus_types_and_conjugacy() {
    assert_eq!(
        stdout(example("hyperelliptic-k3-involution-nikulin")),
        concat!(
            "rvK3LatticePair = -2\n",
            "rvK3LatticePair = -2\n",
            "rvK3AutoNikulin[Id] = identity, r=22, a=0, delta=0, discriminant=(Z/2)^0, fixed locus=whole K3 surface\n",
            "rvK3AutoNikulin[Translate] = symplectic, r=14, a=8, delta=0, discriminant=(Z/2)^8, fixed locus=8 isolated points\n",
            "rvK3AutoNikulin[Enriques] = antisymplectic, r=10, a=10, delta=0, discriminant=(Z/2)^10, fixed locus=empty\n",
            "rvK3AutoNikulin[Gauss] = antisymplectic, r=6, a=6, delta=1, discriminant=(Z/2)^6, fixed locus=1 curve of genus 5\n",
            "rvK3AutoNikulin[Project] = antisymplectic, r=16, a=6, delta=1, discriminant=(Z/2)^6, fixed locus=6 rational curves\n",
            "rvK3AutoNikulin[Conjugate] = antisymplectic, r=16, a=6, delta=1, discriminant=(Z/2)^6, fixed locus=6 rational curves\n",
        )
    );
}

#[test]
fn theta_switches_split_into_six_genus_five_and_ten_free_involutions() {
    let mut source = BASIS.to_string();
    for mask in 0..16 {
        writeln!(
            source,
            "rvk3autom S{mask}[B]{{switch({mask})}}; rvK3AutoNikulin[S{mask}];"
        )
        .unwrap();
    }
    let output = stdout(interpret(&source));
    let reports: Vec<_> = output.lines().collect();
    assert_eq!(reports.len(), 16);
    // The branch classes (including infinity) in this marked Jacobian have
    // masks 0, e1, f1, e1+e2+f1, e1+f1+f2 and their sum.
    for (mask, report) in reports.into_iter().enumerate() {
        let expected = if [0, 1, 4, 7, 13, 15].contains(&mask) {
            "antisymplectic, r=6, a=6, delta=1, discriminant=(Z/2)^6, fixed locus=1 curve of genus 5"
        } else {
            "antisymplectic, r=10, a=10, delta=0, discriminant=(Z/2)^10, fixed locus=empty"
        };
        assert_eq!(report, format!("rvK3AutoNikulin[S{mask}] = {expected}"));
    }
}
