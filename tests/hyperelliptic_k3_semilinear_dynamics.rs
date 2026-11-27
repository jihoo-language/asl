mod common;

use common::{diagnostic, example, interpret, stdout};

const SHEAR: &str = r"
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
fn identity_matches_k3_counts_and_zeta() {
    let output = stdout(interpret(&format!(
        r"{SHEAR}
        rvk3autom Idle[B]{{}};
        rvKummerK3Points[C,1];
        rvK3FrobFixed[Idle,C,1];
        rvKummerK3Points[C,2];
        rvK3FrobFixed[Idle,C,2];
        rvKummerK3Points[C,40];
        rvK3FrobFixed[Idle,C,40];
        rvK3FrobZeta[Idle,C];
        rvK3FrobZeta[Idle,L];
        "
    )));
    let mut lines = output.lines();
    for degree in [1, 2, 40] {
        let expected = lines.next().unwrap().split_once("resolved=").unwrap().1;
        assert_eq!(
            lines.next().unwrap(),
            format!("rvK3FrobFixed[Idle, C, {degree}] = {expected}")
        );
    }
    let over_base = lines.next().unwrap();
    assert!(over_base.contains("(1-9*T)D(3*T)(1-3*T)h(T)"));
    let (marked, cross) = factors(over_base);
    let mut node_factors: Vec<&[i128]> = vec![&[1, -1]; 9];
    node_factors.extend([&[1, 0, -1][..]; 4]);
    assert_eq!(marked, polynomial_product(&node_factors));
    assert_eq!(cross, [1, 8, 30, 72, 81]);

    let over_splitting = lines.next().unwrap();
    assert!(over_splitting.contains("(1-81*T)D(9*T)(1-9*T)h(T)"));
    assert_eq!(
        factors(over_splitting).0,
        polynomial_product(&[&[1, -1][..]; 17])
    );
    assert!(lines.next().is_none());
}

#[test]
fn enriques_frobenius_matches_the_whole_geometric_action_over_f9() {
    let output = stdout(interpret(&format!(
        r"{SHEAR}
        rvk3autom Free[B]{{switch(2)}};
        rvK3FrobFixed[Free,C,1];
        rvK3FrobFixed[Free,C,2];
        rvK3FrobZeta[Free,C];
        let psi=\x:RVKummerK3(L)=>rvK3AutoPoint[Free](rvKummerK3Frob[C,1](x));
        assert (\x:RVKummerK3(L)=>psi(psi(x)))==(\x:RVKummerK3(L)=>x);
        show psi;
        "
    )));
    let mut lines = output.lines();
    assert_eq!(lines.next().unwrap(), "rvK3FrobFixed[Free, C, 1] = 12");
    assert_eq!(lines.next().unwrap(), "rvK3FrobFixed[Free, C, 2] = 248");
    let zeta = lines.next().unwrap();
    assert!(zeta.contains("(1+3*T)h(-T)"));
    assert_eq!(initial_counts(zeta, -1, 2), [12, 248]);
    let map = lines.next().unwrap();
    let points: Vec<_> = map
        .match_indices("RVKummerK3(L){")
        .map(|(start, _)| {
            let end = start + map[start..].find('}').unwrap() + 1;
            &map[start..end]
        })
        .collect();
    assert_eq!(points.len(), 2 * 248);
    assert_eq!(
        points
            .chunks_exact(2)
            .filter(|pair| pair[0] == pair[1])
            .count(),
        12
    );
    assert!(lines.next().is_none());
}

#[test]
fn noncommuting_translation_uses_the_product_before_taking_powers() {
    let output = stdout(interpret(&format!(
        r"{SHEAR}
        rvk3autom Shift[B]{{translate(8)}};
        let psi=\x:RVKummerK3(L)=>rvK3AutoPoint[Shift](rvKummerK3Frob[C,1](x));
        let square=\x:RVKummerK3(L)=>psi(psi(x));
        assert square==(\x:RVKummerK3(L)=>rvKummerK3Translate[B,2](x));
        let x=rvKummerExceptional(rvJacKummerPoint(e1),ow,zw,zw,zw);
        assert rvKummerBlowdown(square(x))==rvJacKummerPoint(add(e1,e2));
        assert square(square(x))==x;
        rvK3FrobFixed[Shift,C,1];
        rvK3FrobFixed[Shift,C,2];
        rvK3FrobFixed[Shift,C,3];
        rvK3FrobFixed[Shift,C,4];
        rvK3FrobZeta[Shift,C];
        "
    )));
    let mut lines = output.lines();
    let counts = [8, 104, 776, 8216];
    for (index, count) in counts.into_iter().enumerate() {
        assert_eq!(
            lines.next().unwrap(),
            format!("rvK3FrobFixed[Shift, C, {}] = {count}", index + 1)
        );
    }
    let zeta = lines.next().unwrap();
    let node_orbit = &[1, 0, 0, 0, -1][..];
    assert_eq!(
        factors(zeta).0,
        polynomial_product(&[&[1, -1], node_orbit, node_orbit, node_orbit, node_orbit])
    );
    assert_eq!(initial_counts(zeta, 1, 4), counts);
    assert!(lines.next().is_none());
}

#[test]
fn runnable_example_includes_a_positive_entropy_word() {
    let output = stdout(example("hyperelliptic-k3-semilinear-dynamics"));
    assert!(output.contains("rvK3FrobFixed[Enriques, C, 1] = 12\n"));
    assert!(output.contains("rvK3FrobFixed[Noncommuting, C, 2] = 104\n"));
    assert!(output.contains("rvK3FrobFixed[Hyperbolic, C, 1] = 84\n"));
    assert!(output.contains("rvK3FrobFixed[Hyperbolic, C, 2] = 3128\n"));
    let zeta = output.lines().last().unwrap();
    let mut marked_factors: Vec<&[i128]> = vec![&[1, -1]; 10];
    marked_factors.extend([&[1, 1][..]; 5]);
    marked_factors.push(&[1, -18, 1]);
    assert_eq!(factors(zeta).0, polynomial_product(&marked_factors));
    assert_eq!(initial_counts(zeta, -1, 2), [84, 3128]);
}

#[test]
fn a_large_conjugating_word_preserves_small_exact_traces() {
    let conjugator = ["project(0),project(1),project(2)"; 12].join(",");
    let inverse = ["project(2),project(1),project(0)"; 12].join(",");
    let output = stdout(interpret(&format!(
        r"{SHEAR}
        rvk3autom Project[B]{{project(3)}};
        rvk3autom Conjugate[B]{{{conjugator},project(3),{inverse}}};
        rvK3FrobFixed[Project,C,40];
        rvK3FrobFixed[Conjugate,C,40];
        rvK3FrobZeta[Project,C];
        rvK3FrobZeta[Conjugate,C];
        "
    )));
    let values: Vec<_> = output
        .lines()
        .map(|line| line.split_once(" = ").unwrap().1)
        .collect();
    assert_eq!(values.len(), 4);
    assert_eq!(values[0], values[1]);
    assert_eq!(values[2], values[3]);
}

#[test]
fn rejects_wrong_descent_and_iteration_ranges() {
    for (program, message) in [
        (
            "rvK3FrobFixed[Idle,C,0];",
            "iteration count must be positive",
        ),
        (
            "rvK3FrobFixed[Idle,C,41];",
            "exceeds unsigned 128-bit range",
        ),
        (
            "rvK3FrobFixed[Idle,C,340282366920938463463374607431768211455];",
            "exceeds unsigned 128-bit range",
        ),
        (
            "rvhyperelliptic Other[F]{o,neg(o),z,z,z,o}; rvK3FrobZeta[Idle,Other];",
            "not the declared base extension",
        ),
    ] {
        let error = diagnostic(interpret(&format!(
            "{SHEAR} rvk3autom Idle[B]{{}}; {program}"
        )));
        assert!(error.contains(message), "{error}");
    }
}

fn factors(report: &str) -> (Vec<i128>, Vec<i128>) {
    fn coefficients(list: &str) -> Vec<i128> {
        list.strip_prefix('[')
            .unwrap()
            .strip_suffix(']')
            .unwrap()
            .split(", ")
            .map(|entry| entry.parse().unwrap())
            .collect()
    }

    let (marked, cross) = report
        .split_once("D=")
        .unwrap()
        .1
        .split_once(", h=")
        .unwrap();
    (coefficients(marked), coefficients(cross))
}

fn polynomial_product(factors: &[&[i128]]) -> Vec<i128> {
    factors.iter().fold(vec![1], |product, factor| {
        let mut result = vec![0; product.len() + factor.len() - 1];
        for (i, &a) in product.iter().enumerate() {
            for (j, &b) in factor.iter().enumerate() {
                result[i + j] += a * b;
            }
        }
        result
    })
}

fn initial_counts(zeta: &str, primitive_sign: i128, order: usize) -> Vec<i128> {
    let (mut marked, mut cross) = factors(zeta);
    for (i, coefficient) in marked.iter_mut().enumerate() {
        *coefficient *= 3i128.pow(i as u32);
    }
    for (i, coefficient) in cross.iter_mut().enumerate() {
        *coefficient *= primitive_sign.pow(i as u32);
    }
    let denominator = polynomial_product(&[
        &[1, -1],
        &[1, -9],
        &marked,
        &[1, -3 * primitive_sign],
        &cross,
    ]);
    // The coefficients of -T D'(T)/D(T) are the fixed-point counts.
    let mut counts = vec![0];
    for k in 1..=order {
        let convolution: i128 = (1..k).map(|i| denominator[i] * counts[k - i]).sum();
        counts.push(-(k as i128) * denominator[k] - convolution);
    }
    counts[1..].to_vec()
}
