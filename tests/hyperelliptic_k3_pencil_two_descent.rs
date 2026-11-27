mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn f13() -> &'static str {
    include_str!("../examples/hyperelliptic-k3-pencil-two-descent.asl")
        .split_once("rvK3PencilTorsion2[P];")
        .unwrap()
        .0
}

const POINT: &str = "rvK3PencilSectionPoint[P]{[0],[1],[0],[6]}{[1]}{[0],[0],[7],[6],[3]}{[1]}";

#[test]
fn torsion_is_deterministic_and_its_kummer_images_obey_the_group_law() {
    let source = format!(
        r"{}
        rvK3PencilTorsion2[P];
        let zeroPoint=zero[RVK3PencilSection(P)];
        let zeroClass=zero[RVK3PencilKummer2(P)];
        let t0=rvK3PencilTorsion2At[P,0];
        let t1=rvK3PencilTorsion2At[P,1];
        let t2=rvK3PencilTorsion2At[P,2];
        assert t0==rvK3PencilSectionPoint[P]{{[0],[1],[9],[6]}}{{[1]}}{{}}{{[1]}};
        assert t1==rvK3PencilSectionPoint[P]{{[0],[6],[8],[10]}}{{[1]}}{{}}{{[1]}};
        assert t2==rvK3PencilSectionPoint[P]{{[0],[8],[6],[9]}}{{[1]}}{{}}{{[1]}};
        assert add(t0,t0)==zeroPoint;
        assert add(t1,t1)==zeroPoint;
        assert add(t2,t2)==zeroPoint;
        assert add(t0,t1)==t2;
        assert add(t0,t2)==t1;
        assert add(t1,t2)==t0;
        assert add(zeroPoint,t0)==t0;
        assert neg(t0)==t0;
        assert rvK3PencilDescent2(zeroPoint)==zeroClass;
        assert rvK3PencilDivisible2(zeroPoint)==S.hi;
        let d0=rvK3PencilDescent2(t0);
        let d1=rvK3PencilDescent2(t1);
        let d2=rvK3PencilDescent2(t2);
        assert add(d0,d1)==d2;
        assert add(d0,d0)==zeroClass;
        assert neg(d2)==d2;
        assert rvK3PencilDivisible2(t0)==S.lo;
        assert rvK3PencilDivisible2(t1)==S.lo;
        assert rvK3PencilDivisible2(t2)==S.lo;
        show d0;
        ",
        f13(),
    );
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "rvK3PencilTorsion2[P] = order=4, affine x-roots=[[[0], [1], [9], [6]], [[0], [6], [8], [10]], [[0], [8], [6], [9]]]
",
            "RVK3PencilKummer2(P)[(square, [[[1], [1]], [[7], [1]], [[9], [1]], [[11], [1]]]), (square, [[[0], [1]], [[7], [1]], [[9], [1]]]), (square, [[[0], [1]], [[1], [1]], [[11], [1]]])]
",
        )
    );
}

#[test]
fn normalized_rational_sections_match_independent_doubling_and_addition() {
    let source = format!(
        r"{}
        let g={POINT};
        let scaled=rvK3PencilSectionPoint[P]{{[0],[2],[2],[12],[12],[0]}}{{[2],[2]}}{{[0],[0],[8],[0],[1],[9]}}{{[3],[3],[0]}};
        assert scaled==g;
        let origin=zero[RVK3PencilSection(P)];
        assert add(g,neg(g))==origin;
        let twice=add(g,g);
        let explicit2=rvK3PencilSectionPoint[P]{{[3],[1],[6],[6],[4]}}{{[1]}}{{[12],[2],[11],[11],[1],[7],[5]}}{{[1]}};
        assert twice==explicit2;
        let triple=add(twice,g);
        let explicit3=rvK3PencilSectionPoint[P]
            {{[0],[3],[5],[2],[11],[6],[6],[10],[6],[3],[6],[6]}}
            {{[3],[0],[12],[0],[7],[0],[3],[0],[1]}}
            {{[0],[0],[5],[3],[8],[8],[6],[7],[11],[4],[1],[5],[9],[7],[3],[10],[9]}}
            {{[12],[0],[7],[0],[10],[0],[2],[0],[9],[0],[11],[0],[1]}};
        assert triple==explicit3;
        assert add(triple,neg(twice))==g;
        assert add(twice,g)==add(g,twice);
        let t0=rvK3PencilTorsion2At[P,0];
        let t1=rvK3PencilTorsion2At[P,1];
        assert add(add(g,t0),t1)==add(g,add(t0,t1));
        assert rvK3PencilDescent2(triple)==rvK3PencilDescent2(g);
        assert rvK3PencilDescent2(neg(g))==rvK3PencilDescent2(g);
        assert rvK3PencilDescent2(add(g,t0))==add(rvK3PencilDescent2(g),rvK3PencilDescent2(t0));
        assert rvK3PencilDivisible2(twice)==S.hi;
        assert rvK3PencilDivisible2(g)==S.lo;
        assert rvK3PencilDivisible2(triple)==S.lo;
        show origin;
        ",
        f13(),
    );
    assert_eq!(
        stdout(interpret(&source)),
        "RVK3PencilSection(P){O}
"
    );
}

#[test]
fn a_section_has_different_prime_to_characteristic_specialization_orders() {
    let mut source = format!("{} let g={POINT}; show g;", f13());
    let parameters = [1, 5, 7, 9, 10, 11];
    for parameter in parameters {
        writeln!(
            source,
            "rvK3FiberJacobian[P](rvK3PencilParam[P](o,rvUnramElem[F]{{{parameter}}}));"
        )
        .unwrap();
    }
    let output = stdout(interpret(&source));
    let mut lines = output.lines();
    assert_eq!(
        lines.next().unwrap(),
        "RVK3PencilSection(P){x=[[0], [1], [0], [6]]/[[1]], y=[[0], [0], [7], [6], [3]]/[[1]]}"
    );
    let orders: Vec<_> = parameters
        .into_iter()
        .zip(lines)
        .map(|(t, fiber)| {
            assert!(fiber.ends_with("status=smooth"));
            let coefficients = ["a2=[", "a4=["].map(|label| {
                fiber
                    .split_once(label)
                    .unwrap()
                    .1
                    .split_once(']')
                    .unwrap()
                    .0
                    .parse()
                    .unwrap()
            });
            let x = residue(t + 6 * t * t * t);
            let y = residue(7 * t * t + 6 * t * t * t + 3 * t * t * t * t);
            let a6: i64 = fiber
                .split_once("a6=[")
                .unwrap()
                .1
                .split_once(']')
                .unwrap()
                .0
                .parse()
                .unwrap();
            assert_eq!(
                residue(y * y),
                residue(x * x * x + coefficients[0] * x * x + coefficients[1] * x + a6)
            );
            let mut acc = None;
            for order in 1..=20 {
                acc = add_affine(coefficients, acc, (x, y));
                if acc.is_none() {
                    return order;
                }
            }
            unreachable!("finite elliptic order exceeds its Hasse bound")
        })
        .collect();
    // Reduction is injective on prime-to-13 torsion. These orders certify
    // infinite order of the supplied section, irrespective of 13-primary torsion.
    assert_eq!(orders, [3, 4, 6, 6, 4, 3]);
}

fn residue(n: i64) -> i64 {
    n.rem_euclid(13)
}

fn add_affine(
    [a2, a4]: [i64; 2],
    left: Option<(i64, i64)>,
    (u, v): (i64, i64),
) -> Option<(i64, i64)> {
    let Some((x, y)) = left else {
        return Some((u, v));
    };
    if x == u && residue(y + v) == 0 {
        return None;
    }
    let (n, d) = if x == u {
        (3 * x * x + 2 * a2 * x + a4, 2 * y)
    } else {
        (v - y, u - x)
    };
    let inverse = (1..13).find(|i| residue(i * d) == 1).unwrap();
    let slope = residue(n * inverse);
    let result_x = residue(slope * slope - a2 - x - u);
    Some((result_x, residue(slope * (x - result_x) - y)))
}

#[test]
fn nonsquare_constants_and_characteristic_three_are_supported() {
    let prime = f13()
        .replace("F[13,1]", "F[5,1]")
        .replace("z,eleven,two,nine,three,o", "z,neg(o),z,z,z,o");
    let output = stdout(interpret(&format!(
        r"{prime}
        rvK3PencilTorsion2[P];
        show rvK3PencilDescent2(rvK3PencilTorsion2At[P,1]);
        "
    )));
    assert!(output.contains("order=4"));
    assert!(output.contains("(nonsquare,"));

    let extension = r"
        rvunram F[3,1]{1,0,1};
        let z=rvUnramElem[F]{};
        let o=rvUnramElem[F]{1};
        let i=rvUnramElem[F]{0,1};
        rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
        let r0=rvJacPoint[C](z,z);
        let r1=rvJacPoint[C](o,z);
        let r2=rvJacPoint[C](neg(o),z);
        let r3=rvJacPoint[C](i,z);
        rvjactorsionbasis B[C,2](neg(o)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),r3)};
        rvk3pencil P[B,0,1];
        rvK3PencilTorsion2[P];
        let t0=rvK3PencilTorsion2At[P,0];
        let t1=rvK3PencilTorsion2At[P,1];
        let t2=rvK3PencilTorsion2At[P,2];
        assert add(t0,t1)==t2;
        assert add(t0,t0)==zero[RVK3PencilSection(P)];
        assert add(rvK3PencilDescent2(t0),rvK3PencilDescent2(t1))==rvK3PencilDescent2(t2);
    ";
    assert!(stdout(interpret(extension)).contains("order=4"));
}

#[test]
fn rejects_invalid_coordinates_and_cross_pencil_arithmetic() {
    for (source, error) in [
        (
            "show rvK3PencilSectionPoint[P]{[0]}{}{[0]}{[1]};",
            "zero denominator",
        ),
        (
            "show rvK3PencilSectionPoint[P]{[1]}{[1]}{[0]}{[1]};",
            "does not lie on the generic K3 pencil Jacobian",
        ),
        (
            "show rvK3PencilSectionPoint[P]{[0,1]}{[1]}{[0]}{[1]};",
            "at most 1 coordinates",
        ),
        (
            "show rvK3PencilTorsion2At[P,3];",
            "affine two-torsion index is out of range",
        ),
        (
            "show rvK3PencilDivisible2(o);",
            "an RVK3PencilSection point",
        ),
        (
            "rvk3pencil Q[B,1,0]; show add(zero[RVK3PencilSection(P)],zero[RVK3PencilSection(Q)]);",
            "expected type RVK3PencilSection(P)",
        ),
        (
            "rvk3pencil Q[B,1,0]; show add(zero[RVK3PencilKummer2(P)],zero[RVK3PencilKummer2(Q)]);",
            "expected type RVK3PencilKummer2(P)",
        ),
    ] {
        let message = diagnostic(interpret(&format!("{} {source}", f13())));
        assert!(message.contains(error), "{message}");
    }
}

#[test]
fn executable_example_shows_nontrivial_and_divisible_sections() {
    let output = stdout(example("hyperelliptic-k3-pencil-two-descent"));
    assert!(output.contains("order=4"));
    assert!(output.contains("RVK3PencilKummer2(P)"));
    assert!(output.ends_with(
        "S.lo
"
    ));
}
