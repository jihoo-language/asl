mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const MARKED: &str = r"
    rvunram W[3,1]{1,0,1};
    let z=rvUnramElem[W]{};
    let o=rvUnramElem[W]{1};
    let i=rvUnramElem[W]{0,1};
    rvhyperelliptic C[W]{z,neg(o),z,z,z,o};
    let e1=rvJacPoint[C](z,z);
    let f1=rvJacPoint[C](o,z);
    let e2=add(add(e1,f1),rvJacPoint[C](neg(o),z));
    let f2=add(add(e1,f1),rvJacPoint[C](i,z));
    rvjactorsionbasis B[C,2](neg(o)){e1,e2,f1,f2};
";

#[test]
fn chern_numbers_todd_genus_and_exceptional_fourth_power() {
    assert_eq!(
        stdout(example("hyperelliptic-hilbert-square-ring")),
        concat!(
            "rvHilbSquareIntegral = 12\n",
            "rvHilbSquareIntegral = 828\n",
            "rvHilbSquareIntegral = -60\n",
            "rvHilbSquareIntegral = 324\n",
            "rvHilbSquareIntegral = 3\n",
        )
    );
}

#[test]
fn actual_integral_marking_and_every_symmetric_monomial_use_the_advertised_coordinates() {
    let mut source = MARKED.to_owned();
    for i in 0..23 {
        writeln!(source, "let x{i}=rvHilbSquareBasis[B,2,{i}];").unwrap();
    }
    for i in 0..22 {
        writeln!(
            source,
            "assert rvHilbSquareDiv[0](rvK3LatticeBasis[B,{i}])==x{i};"
        )
        .unwrap();
    }
    source.push_str("assert rvHilbSquareDiv[1](zero[RVK3Lattice(B)])==x22;\n");
    let mut index = 0;
    for i in 0..23 {
        for j in i..23 {
            writeln!(
                source,
                "assert rvHilbSquareCup(x{j},x{i})==rvHilbSquareBasis[B,4,{index}];"
            )
            .unwrap();
            index += 1;
        }
    }
    assert_eq!(index, 276);
    source.push_str(
        r"
        let theta=rvK3Lattice[B]{[0,2,0,0,2,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]};
        let theta2=rvHilbSquareCup(rvHilbSquareDiv[0](theta),rvHilbSquareDiv[0](theta));
        rvK3LatticePair(theta,theta);
        rvHilbSquareIntegral(rvHilbSquareCup(theta2,theta2));
        rvHilbSquareIntegral(rvHilbSquareCup(theta2,rvHilbSquareCup(x22,x22)));
        rvHilbSquareIntegral(rvHilbSquareCup(rvHilbSquareC2[B],theta2));
        assert rvHilbSquareDiv[-2](add(rvK3LatticeBasis[B,1],rvK3LatticeBasis[B,19]))
            ==add(add(x1,x19),rvHilbSquareScale[-2,1](x22));
    ",
    );
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "rvK3LatticePair = 4\n",
            "rvHilbSquareIntegral = 48\n",
            "rvHilbSquareIntegral = -8\n",
            "rvHilbSquareIntegral = 120\n"
        )
    );
}

#[test]
fn associative_products_and_c2_contraction_hold_for_rational_linear_combinations() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
            let x=rvHilbSquareScale[-2,3](add(rvHilbSquareBasis[B,2,0],rvHilbSquareBasis[B,2,22]));
            let y=rvHilbSquareScale[7,5](add(rvHilbSquareBasis[B,2,5],rvHilbSquareBasis[B,2,22]));
            let w=rvHilbSquareScale[4,-7](add(rvHilbSquareBasis[B,2,7],rvHilbSquareBasis[B,2,11]));
            let yDual=rvHilbSquareScale[7,5](add(rvHilbSquareBasis[B,6,5],rvHilbSquareBasis[B,6,22]));
            let c=rvHilbSquareC2[B];
            let xy=rvHilbSquareCup(x,y);
            let zz=rvHilbSquareCup(w,w);
            assert rvHilbSquareCup(xy,w)==rvHilbSquareCup(x,rvHilbSquareCup(y,w));
            assert rvHilbSquareCup(rvHilbSquareCup(xy,w),w)==rvHilbSquareCup(xy,zz);
            assert rvHilbSquareCup(xy,c)==rvHilbSquareCup(x,rvHilbSquareCup(y,c));
            assert rvHilbSquareCup(xy,c)==rvHilbSquareScale[30,1](rvHilbSquareCup(x,yDual));
            assert rvHilbSquareCup(add(x,y),w)==add(rvHilbSquareCup(x,w),rvHilbSquareCup(y,w));
            assert add(x,neg(x))==zero[RVHilbSquare(B,2)];
            assert rvHilbSquareCup(rvHilbSquareBasis[B,0,0],xy)==xy;
            let point=rvHilbSquareBasis[B,8,0];
            assert rvHilbSquareCup(rvHilbSquareScale[2,3](rvHilbSquareBasis[B,0,0]),point)==rvHilbSquareScale[2,3](point);
            assert rvHilbSquareCup(x,point)==zero[RVHilbSquare(B,10)];
            assert rvHilbSquareCup(zero[RVHilbSquare(B,1)],zero[RVHilbSquare(B,1)])==zero[RVHilbSquare(B,2)];
            rvHilbSquareIntegral(rvHilbSquareScale[-6,9](point));
            "
        ))),
        "rvHilbSquareIntegral = -2/3\n"
    );
}

#[test]
fn weighted_cup_product_cancels_before_multiplying_and_annihilates_zero() {
    let maximum = i128::MAX;
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
            let d=rvHilbSquareBasis[B,2,22];
            let d2=rvHilbSquareCup(d,d);
            let large=rvHilbSquareScale[{maximum},1](d2);
            let small=rvHilbSquareScale[1,12](d2);
            assert rvHilbSquareCup(zero[RVHilbSquare(B,2)],large)==zero[RVHilbSquare(B,6)];
            assert rvHilbSquareCup(d,rvHilbSquareBasis[B,6,22])
                ==rvHilbSquareScale[-2,1](rvHilbSquareBasis[B,8,0]);
            rvHilbSquareIntegral(rvHilbSquareCup(large,small));
            "
        ))),
        format!("rvHilbSquareIntegral = {maximum}\n")
    );
}

#[test]
fn dimensions_and_nominal_errors_are_checked() {
    for (degree, dimension) in [(0, 1), (2, 23), (4, 276), (6, 23), (8, 1), (1, 0), (10, 0)] {
        if dimension > 0 {
            assert_eq!(
                stdout(interpret(&format!(
                    "{MARKED} let x=rvHilbSquareBasis[B,{degree},{}];",
                    dimension - 1
                ))),
                ""
            );
        }
        let error = diagnostic(interpret(&format!(
            "{MARKED} show rvHilbSquareBasis[B,{degree},{dimension}];"
        )));
        assert!(
            error.contains(&format!("degree {degree} of dimension {dimension}")),
            "{error}"
        );
    }
    for (source, expected) in [
        ("show rvHilbSquareScale[1,0](rvHilbSquareBasis[B,0,0]);", "denominator cannot vanish"),
        ("rvHilbSquareIntegral(rvHilbSquareBasis[B,6,0]);", "requires cohomological degree eight"),
        ("show rvHilbSquareDiv[0](zero[RVJac(C)]);", "an RVK3Lattice vector"),
        ("show rvHilbSquareCup(zero[RVJac(C)],rvHilbSquareBasis[B,2,0]);", "an RVHilbSquare cohomology class"),
        ("rvjactorsionbasis D[C,2](neg(o)){e1,e2,f1,f2}; show rvHilbSquareCup(zero[RVHilbSquare(B,2)],zero[RVHilbSquare(D,4)]);", "requires the same named marking"),
        ("show add(zero[RVHilbSquare(B,2)],zero[RVHilbSquare(B,4)]);", "expected type RVHilbSquare(B, 2)"),
        ("show rvHilbSquareC2[Missing];", "unknown Jacobian torsion basis"),
        ("show rvHilbSquareScale[2,1](rvHilbSquareScale[170141183460469231731687303715884105727,1](rvHilbSquareBasis[B,8,0]));", "exceeds the signed 128-bit coefficient range"),
    ] {
        let error = diagnostic(interpret(&format!("{MARKED}{source}")));
        assert!(error.contains(expected), "{error}");
    }
}
