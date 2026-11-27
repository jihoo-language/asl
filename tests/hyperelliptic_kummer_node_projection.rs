mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const F3: &str = r"
    rvunram F[3,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let two=rvUnramElem[F]{2};
";

const SPLIT: &str = r"
    rvhyperelliptic C[F]{z,two,z,z,z,o};
    rvunram W[3,1]{1,0,1};
    rvjacbase D[C,W];
    let zw=rvUnramElem[W]{};
    let ow=rvUnramElem[W]{1};
    let tw=neg(ow);
    let i=rvUnramElem[W]{0,1};
    let t0=rvJacPoint[D](zw,zw);
    let t1=rvJacPoint[D](ow,zw);
    let tm=rvJacPoint[D](tw,zw);
    let ti=rvJacPoint[D](i,zw);
    rvjactorsionbasis B[D,2](tw){t0,add(add(t0,t1),tm),t1,add(add(t0,t1),ti)};
    let identity=zero[RVKummer(D)];
";

#[test]
fn exact_line_and_tangent_images_in_characteristic_three() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-node-projection")),
        concat!(
            "RVKummerK3(C){base=[[1], [0], [1], [2]]}\n",
            "RVKummerK3(C){base=[[1], [1], [1], [0]]}\n",
            "RVKummerK3(C){base=[[0], [1], [0], [0]], direction=[[1], [0], [1], [2]]}\n",
        )
    );
}

#[test]
fn every_f3_and_f9_center_acts_involutively_on_its_full_resolution() {
    let mut source = format!("{F3}{SPLIT}");
    source.push_str(
        r"
        rvKummerK3Points[C,1];
        rvKummerK3Points[D,1];
        let c0=rvJacPoint[C](z,z);
        let c1=rvJacPoint[C](o,z);
        let cm=rvJacPoint[C](two,z);
        let cqm=rvjac[C]{[1],[],[1]}{};
        ",
    );
    for center in [
        "zero[RVJac(C)]",
        "c0",
        "c1",
        "cm",
        "add(c0,c1)",
        "add(c0,cm)",
        "add(c1,cm)",
        "cqm",
    ] {
        writeln!(
            source,
            r"assert (\x:RVKummerK3(C)=>rvKummerNodeInvolution(rvJacKummerPoint({center}),
                rvKummerNodeInvolution(rvJacKummerPoint({center}),x)))==(\x:RVKummerK3(C)=>x);"
        )
        .unwrap();
    }
    for mask in 0..16 {
        writeln!(source, "let n{mask}=rvKummerTranslate[B,{mask}](identity);").unwrap();
        writeln!(
            source,
            r"assert (\x:RVKummerK3(D)=>rvKummerNodeInvolution(n{mask},rvKummerNodeInvolution(n{mask},x)))
                ==(\x:RVKummerK3(D)=>x);"
        )
        .unwrap();
    }
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "rvKummerK3Points[C, 1] = singular=8, nodes=8, resolved=32\n",
            "rvKummerK3Points[D, 1] = singular=104, nodes=16, resolved=248\n",
        )
    );
}

#[test]
fn projection_commutes_with_field_transport_and_conjugates_by_node_translation() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F3}{SPLIT}
            let baseCenter=rvJacKummerPoint(rvJacPoint[C](o,z));
            let extended=rvKummerExtend[D](baseCenter);
            assert (\x:RVKummerK3(C)=>rvKummerK3Extend[D](rvKummerNodeInvolution(baseCenter,x)))
                ==(\x:RVKummerK3(C)=>rvKummerNodeInvolution(extended,rvKummerK3Extend[D](x)));

            let node=rvKummerTranslate[B,13](identity);
            let conjugate=rvKummerTranslate[B,15](identity);
            assert rvKummerFrob[C,1](node)==conjugate;
            assert (\x:RVKummerK3(D)=>rvKummerK3Frob[C,1](rvKummerNodeInvolution(node,x)))
                ==(\x:RVKummerK3(D)=>rvKummerNodeInvolution(conjugate,rvKummerK3Frob[C,1](x)));

            let shifted=rvKummerTranslate[B,9](identity);
            assert rvKummerTranslate[B,4](node)==shifted;
            assert (\x:RVKummerK3(D)=>rvKummerK3Translate[B,4](
                rvKummerNodeInvolution(node,rvKummerK3Translate[B,4](x))))
                ==(\x:RVKummerK3(D)=>rvKummerNodeInvolution(shifted,x));
            "
        ))),
        ""
    );
}

#[test]
fn unmarked_projections_also_square_to_identity_over_f7() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram F[7,1]{0,1};
            let z=rvUnramElem[F]{};
            let o=rvUnramElem[F]{1};
            let three=rvUnramElem[F]{3};
            let four=rvUnramElem[F]{4};
            rvhyperelliptic C[F]{z,three,neg(o),z,four,o};
            rvKummerK3Points[C,1];
            let n0=zero[RVKummer(C)];
            let p=rvJacPoint[C](z,z);
            let n1=rvJacKummerPoint(p);
            let n2=rvJacKummerPoint(add(p,rvJacPoint[C](o,z)));
            assert (\x:RVKummerK3(C)=>rvKummerNodeInvolution(n0,rvKummerNodeInvolution(n0,x)))
                ==(\x:RVKummerK3(C)=>x);
            assert (\x:RVKummerK3(C)=>rvKummerNodeInvolution(n1,rvKummerNodeInvolution(n1,x)))
                ==(\x:RVKummerK3(C)=>x);
            assert (\x:RVKummerK3(C)=>rvKummerNodeInvolution(n2,rvKummerNodeInvolution(n2,x)))
                ==(\x:RVKummerK3(C)=>x);
            "
        )),
        "rvKummerK3Points[C, 1] = singular=48, nodes=16, resolved=160\n"
    );
}

#[test]
fn projection_requires_a_node_and_a_resolved_point_on_the_same_named_curve() {
    let setup = r"
        rvhyperelliptic C[F]{z,o,z,z,z,o};
        let center=zero[RVKummer(C)];
        let smooth=rvKummer[C](o,z,o,o);
        let x=rvKummerSmooth(smooth);
    ";
    for (input, expected) in [
        ("show rvKummerNodeInvolution(smooth,x);", "center must be a Kummer node"),
        ("show rvKummerNodeInvolution(center,smooth);", "expected type RVKummerK3(C)"),
        ("show rvKummerNodeInvolution(x,x);", "an RVKummer point"),
        (
            "rvhyperelliptic Twin[F]{z,o,z,z,z,o}; show rvKummerNodeInvolution(zero[RVKummer(Twin)],x);",
            "expected type RVKummerK3(Twin), found RVKummerK3(C)",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{F3}{setup}{input}")));
        assert!(error.contains(expected), "{error}");
    }
}
