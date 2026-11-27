mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn setup() -> &'static str {
    include_str!("../examples/hyperelliptic-ks-torsion-linking.asl")
        .split_once("// Linking in the oriented")
        .unwrap()
        .0
}

#[test]
fn example_certifies_low_and_complementary_degree_orders() {
    assert_eq!(
        stdout(example("hyperelliptic-ks-torsion-linking")),
        concat!(
            "rvQZOrder = 3\n",
            "rvKSLinkCertificate = value=1/3, pairing_order=3, left_exact=true, right_exact=true\n",
            "rvKSLinkCertificate = value=5/6, pairing_order=6, left_exact=true, right_exact=true\n",
        )
    );
}

#[test]
fn representative_changes_preimage_changes_and_bilinearity_are_exact() {
    let source = format!(
        r"{}
        let original=rvKSLink[A,1,3,3](u,w,v,pre);
        let u2=rvKSExteriorBasis[P,128];
        let w2=rvKSExteriorBasis[P,64];
        let v2=rvKSExteriorBasis[P,18446744073709551551];
        let pre2=rvKSExteriorBasis[P,18446744073709551487];
        assert rvKSExteriorDifference[A,1,1](w2)==rvKSExteriorScale[3](u2);
        assert rvKSExteriorDifference[A,1,1](pre2)==rvKSExteriorScale[3](v2);
        let second=rvKSLink[A,1,3,3](u2,w2,v2,pre2);
        assert second==rvQZ[1,3];
        assert rvKSLink[A,1,3,3](u,w,v2,pre2)==zero[RVQZ];
        assert rvKSLink[A,1,3,3](u2,w2,v,pre)==zero[RVQZ];
        let total=rvKSLink[A,1,3,3](add(u,u2),add(w,w2),add(v,v2),add(pre,pre2));
        assert total==add(original,second);
        assert total==rvQZScale[2](original);
        assert rvKSLink[A,1,3,3](neg(u),neg(w),v,pre)==neg(original);
        let movedu=add(u,rvKSExteriorDifference[A,1,1](w));
        let movedw=add(rvKSExteriorScale[4](w),u);
        let movedv=add(v,rvKSExteriorDifference[A,1,1](pre));
        let movedpre=add(rvKSExteriorScale[4](pre),v);
        assert rvKSLink[A,1,3,3](movedu,movedw,movedv,movedpre)==original;
        assert rvKSLink[A,1,6,9](u,rvKSExteriorScale[2](w),v,rvKSExteriorScale[3](pre))==original;
        let coboundary=rvKSExteriorDifference[A,1,1](w);
        assert rvKSLink[A,1,1,3](coboundary,w,v,pre)==zero[RVQZ];
        let coboundaryRight=rvKSExteriorDifference[A,1,1](pre);
        assert rvKSLink[A,1,3,1](u,w,coboundaryRight,pre)==zero[RVQZ];
        rvKSLinkCertificate[A,1,6,3](u,rvKSExteriorScale[2](w),v,pre);
        ",
        setup()
    );
    assert_eq!(
        stdout(interpret(&source)),
        "rvKSLinkCertificate = value=1/3, pairing_order=3, left_exact=false, right_exact=true\n"
    );
}

#[test]
fn exterior_orientation_produces_the_odd_exchange_sign() {
    let source = format!(
        r"{}
        let source=rvKSExteriorBasis[P,5];
        let target=rvKSExteriorBasis[P,6];
        let other=rvKSExteriorBasis[P,18446744073709551610];
        let otherSource=rvKSExteriorBasis[P,18446744073709551609];
        assert rvKSExteriorDifference[A,1,1](source)==rvKSExteriorScale[3](target);
        assert rvKSExteriorDifference[A,1,1](otherSource)==rvKSExteriorScale[3](other);
        let direct=rvKSLink[A,1,3,3](target,source,other,otherSource);
        let reverse=rvKSLink[A,1,3,3](other,otherSource,target,source);
        assert direct==rvQZ[-1,3];
        assert reverse==neg(direct);
        ",
        setup()
    );
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn middle_linking_is_alternating_on_genuine_order_two_classes() {
    let active = [0, 6, 10, 12, 18, 20, 24, 30, 34, 36, 40, 46, 48, 54, 58, 60];
    let occupied: u64 = active.iter().map(|i| 3u64 << i).sum();
    let inert: u64 = (0..64)
        .map(|i| 1u64 << i)
        .filter(|bit| bit & occupied == 0)
        .take(31)
        .sum();
    let pre1 = inert | 1;
    let image1 = inert | 2;
    let pre2 = !image1;
    let image2 = !pre1;
    let source = format!(
        r"{}
        let pre1=rvKSExteriorBasis[P,{pre1}];
        let pre2=rvKSExteriorBasis[P,{pre2}];
        let image1=rvKSExteriorBasis[P,{image1}];
        let image2=rvKSExteriorBasis[P,{image2}];
        let class1=rvKSExteriorScale[3](image1);
        let class2=rvKSExteriorScale[3](image2);
        assert rvKSExteriorDifference[A,2,1](pre1)==rvKSExteriorScale[2](class1);
        assert rvKSExteriorDifference[A,2,1](pre2)==rvKSExteriorScale[2](class2);
        let forward=rvKSLink[A,2,2,2](class1,pre1,class2,pre2);
        let backward=rvKSLink[A,2,2,2](class2,pre2,class1,pre1);
        assert forward==rvQZ[1,2];
        assert backward==neg(forward);
        assert rvKSLink[A,2,2,2](class1,pre1,class1,pre1)==zero[RVQZ];
        assert rvKSLink[A,2,2,2](class2,pre2,class2,pre2)==zero[RVQZ];
        let sum=add(class1,class2);
        let sumPre=add(pre1,pre2);
        assert rvKSLink[A,2,2,2](sum,sumPre,sum,sumPre)==zero[RVQZ];
        rvKSLinkCertificate[A,2,2,2](class1,pre1,class2,pre2);
        ",
        setup()
    );
    assert_eq!(
        stdout(interpret(&source)),
        "rvKSLinkCertificate = value=1/2, pairing_order=2, left_exact=true, right_exact=true\n"
    );

    // Two singly occupied pairs create a lift whose self numerator is
    // nonzero before taking Q/Z. Its integral self-evaluation uses parity.
    let double_inert = inert & (inert - 1);
    let other_inert = !(double_inert | 195);
    let source = format!(
        r"{}
        let lift=add(rvKSExteriorBasis[P,{}],rvKSExteriorBasis[P,{}]);
        let single1=add(rvKSExteriorBasis[P,{}],rvKSExteriorBasis[P,{}]);
        let single2=add(rvKSExteriorBasis[P,{}],rvKSExteriorBasis[P,{}]);
        let double=add(rvKSExteriorBasis[P,{}],rvKSExteriorBasis[P,{}]);
        let cls=add(rvKSExteriorScale[3](add(single1,single2)),rvKSExteriorScale[18](double));
        assert rvKSExteriorDifference[A,2,1](lift)==rvKSExteriorScale[2](cls);
        let partner=rvKSExteriorScale[3](rvKSExteriorBasis[P,{}]);
        let partnerPre=rvKSExteriorBasis[P,{}];
        rvKSLinkCertificate[A,2,2,2](cls,lift,partner,partnerPre);
        assert rvKSLink[A,2,2,2](cls,lift,cls,lift)==zero[RVQZ];
        rvKSExteriorCoeff[18446744073709551615](rvKSExteriorCup(cls,cls));
        rvKSExteriorCoeff[18446744073709551615](rvKSExteriorCup(lift,cls));
        ",
        setup(),
        double_inert | 65,
        other_inert | 65,
        double_inert | 66,
        other_inert | 66,
        double_inert | 129,
        other_inert | 129,
        double_inert | 130,
        other_inert | 130,
        other_inert | 130,
        other_inert | 66,
    );
    let actual = stdout(interpret(&source));
    let mut lines = actual.lines();
    assert_eq!(
        lines.next().unwrap(),
        "rvKSLinkCertificate = value=1/2, pairing_order=2, left_exact=true, right_exact=true"
    );
    let mut coeffs = lines.map(|line| line.split_once(" = ").unwrap().1.parse::<i32>().unwrap());
    let square = coeffs.next().unwrap();
    assert_eq!(square.abs(), 36);
    assert_eq!(coeffs.next(), Some(-square));
}

#[test]
fn wide_witness_products_and_reusable_qz_scalars_are_exact() {
    let cover = -(1i128 << 126);
    let bound = cover.unsigned_abs() * 3;
    let big = i128::MAX;
    let source = format!(
        r"{}
        let ordinary=rvKSLink[A,{cover},{bound},{bound}](u,neg(w),v,neg(pre));
        assert ordinary==rvQZ[-1,{bound}];
        let hugeU=rvKSExteriorScale[{big}](rvKSExteriorScale[{big}](u));
        let hugeW=rvKSExteriorScale[{big}](rvKSExteriorScale[{big}](neg(w)));
        let wide=rvKSLink[A,{cover},{bound},{bound}](hugeU,hugeW,v,neg(pre));
        assert wide==rvQZScale[{big}](rvQZScale[{big}](ordinary));
        rvQZOrder(ordinary);
        let a=rvQZ[1,340282366920938463463374607431768211455];
        let b=rvQZ[1,340282366920938463463374607431768211454];
        let sum=add(a,b);
        assert add(sum,neg(b))==a;
        assert add(sum,neg(sum))==zero[RVQZ];
        assert rvQZ[-7,6]==rvQZ[5,6];
        assert rvQZScale[-170141183460469231731687303715884105728](rvQZ[1,170141183460469231731687303715884105728])==zero[RVQZ];
        assert add(rvQZ[1,2],rvQZ[1,3])==rvQZ[5,6];
        rvQZOrder(zero[RVQZ]);
        ",
        setup()
    );
    assert_eq!(
        stdout(interpret(&source)),
        format!("rvQZOrder = {bound}\nrvQZOrder = 1\n")
    );
}

#[test]
fn zero_and_witness_diagnostics_respect_cover_degree_and_marking() {
    let source = format!(
        r"{}
        let z1=zero[RVKSExterior(P,1)];
        let z63=zero[RVKSExterior(P,63)];
        assert rvKSLink[A,0,3,7](z1,w,z63,pre)==zero[RVQZ];
        assert rvKSLink[A,1,3,3](z1,u,z63,v)==zero[RVQZ];
        let unit=rvKSExteriorBasis[P,0];
        let top=rvKSExteriorBasis[P,18446744073709551615];
        assert rvKSLink[A,1,1,1](zero[RVKSExterior(P,0)],unit,zero[RVKSExterior(P,64)],top)==zero[RVQZ];
        ",
        setup()
    );
    assert_eq!(stdout(interpret(&source)), "");
    let invalid = [
        (
            "show rvKSLink[A,1,0,3](u,w,v,pre);",
            "witness bounds must be positive",
        ),
        (
            "show rvKSLink[A,2,3,6](u,w,v,pre);",
            "left Kuga–Satake linking witness",
        ),
        (
            "show rvKSLink[A,1,3,3](u,w,v,v);",
            "right Kuga–Satake linking witness",
        ),
        (
            "show rvKSLink[A,0,3,3](u,w,v,pre);",
            "identity monodromy at m=0",
        ),
        (
            "show rvKSLink[A,1,3,3](u,w,u,w);",
            "expected type RVKSExterior(P,63)",
        ),
        (
            "show rvKSLink[A,1,3,3](zero[RVKSExterior(P,65)],w,v,pre);",
            "fiber degree must be between 0 and 64",
        ),
        (
            "rvog6cusp Q(h,e,f);show rvKSLink[A,1,3,3](rvKSExteriorBasis[Q,2],w,v,pre);",
            "expected type RVKSExterior(P,1)",
        ),
        (
            "show rvKSLink[Missing,1,3,3](u,w,v,pre);",
            "unknown OG6 nilpotent",
        ),
        ("show rvQZ[1,0];", "denominator must be positive"),
        ("rvQZOrder(u);", "expected type RVQZ"),
    ];
    for (source, error) in invalid {
        let mut program = setup().to_owned();
        writeln!(program, "{source}").unwrap();
        let actual = diagnostic(interpret(&program));
        assert!(actual.contains(error), "{actual}");
    }
}
