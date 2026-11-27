mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const F3: &str = r"
    rvunram F[3,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let two=rvUnramElem[F]{2};
";

const F7: &str = r"
    rvunram F[7,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let two=rvUnramElem[F]{2};
    let three=rvUnramElem[F]{3};
    let four=rvUnramElem[F]{4};
    let five=rvUnramElem[F]{5};
    let six=neg(o);
    rvhyperelliptic C[F]{z,three,six,z,four,o};
    let t0=rvJacPoint[C](z,z);
    let t1=rvJacPoint[C](o,z);
    let t2=rvJacPoint[C](two,z);
    let t3=rvJacPoint[C](three,z);
    rvjactorsionbasis B[C,2](six){t0,add(add(t0,t1),t2),t1,add(add(t0,t1),t3)};
    let identity=zero[RVKummer(C)];
";

#[test]
fn the_example_has_exact_images_on_the_smooth_and_exceptional_charts() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-switch")),
        concat!(
            "RVKummerK3(C){base=[[1], [2], [2], [0]]}\n",
            "RVKummerK3(D){base=[[0, 0], [0, 0], [0, 0], [1, 0]], direction=[[1, 0], [2, 2], [0, 2], [0, 0]]}\n",
        )
    );
}

#[test]
fn the_switch_is_an_involution_and_commutes_with_extension_and_frobenius() {
    let output = stdout(interpret(&format!(
        r"{F3}
        rvhyperelliptic C[F]{{z,two,z,z,z,o}};
        rvunram W[3,1]{{1,0,1}};
        rvjacbase D[C,W];
        rvKummerK3Points[C,1];
        rvKummerK3Points[D,1];
        assert (\x:RVKummerK3(C)=>rvKummerSwitch(rvKummerSwitch(x)))==(\x:RVKummerK3(C)=>x);
        assert (\x:RVKummerK3(D)=>rvKummerSwitch(rvKummerSwitch(x)))==(\x:RVKummerK3(D)=>x);
        assert (\x:RVKummerK3(C)=>rvKummerSwitch(rvKummerK3Extend[D](x)))
            ==(\x:RVKummerK3(C)=>rvKummerK3Extend[D](rvKummerSwitch(x)));
        assert (\x:RVKummerK3(D)=>rvKummerSwitch(rvKummerK3Frob[C,1](x)))
            ==(\x:RVKummerK3(D)=>rvKummerK3Frob[C,1](rvKummerSwitch(x)));
        show (\x:RVKummerK3(D)=>(rvKummerIsExceptional(x),rvKummerIsExceptional(rvKummerSwitch(x))));
        "
    )));
    let mut lines = output.lines();
    assert_eq!(
        lines.next(),
        Some("rvKummerK3Points[C, 1] = singular=8, nodes=8, resolved=32")
    );
    assert_eq!(
        lines.next(),
        Some("rvKummerK3Points[D, 1] = singular=104, nodes=16, resolved=248")
    );
    let transitions = lines.next().unwrap();
    for (label, count) in [
        ("(S.lo, S.lo)", 24),
        ("(S.lo, S.hi)", 64),
        ("(S.hi, S.lo)", 64),
        ("(S.hi, S.hi)", 96),
    ] {
        assert_eq!(transitions.matches(label).count(), count, "{label}");
    }
    assert!(lines.next().is_none());
}

fn squarefree_over_f3(coefficients: &[i32; 6]) -> bool {
    // A repeated factor of a quintic has degree at most two, so it
    // shares a root with the derivative over F9 = F3[i]/(i²+1).
    for a in 0..3 {
        for b in 0..3 {
            let horner = |[r, s]: [i32; 2], [u, v]: [i32; 2]| {
                [
                    (a * r - b * s + u).rem_euclid(3),
                    (a * s + b * r + v).rem_euclid(3),
                ]
            };
            let mut value = [0, 0];
            let mut derivative = [0, 0];
            for &coefficient in coefficients.iter().rev() {
                derivative = horner(derivative, value);
                value = horner(value, [coefficient, 0]);
            }
            if value == [0, 0] && derivative == [0, 0] {
                return false;
            }
        }
    }
    true
}

#[test]
fn the_unmarked_switch_is_total_on_every_smooth_monic_quintic_over_f3() {
    let mut source = F3.to_owned();
    let mut count = 0;
    for (index, mut digits) in (0..243).enumerate() {
        let coefficients = std::array::from_fn(|i| {
            if i == 5 {
                1
            } else {
                let coefficient = digits % 3;
                digits /= 3;
                coefficient
            }
        });
        if !squarefree_over_f3(&coefficients) {
            continue;
        }
        count += 1;
        let coefficients = coefficients.map(|coefficient| ["z", "o", "two"][coefficient as usize]);
        writeln!(
            source,
            "rvhyperelliptic C{index}[F]{{{}}};",
            coefficients.join(",")
        )
        .unwrap();
        writeln!(
            source,
            r"assert (\x:RVKummerK3(C{index})=>rvKummerSwitch(rvKummerSwitch(x)))==(\x:RVKummerK3(C{index})=>x);"
        )
        .unwrap();
    }
    assert_eq!(count, 162);
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn all_split_exceptional_curves_and_strict_tropes_are_exchanged() {
    let mut source = F7.to_owned();
    source.push_str(
        r"
        rvKummerK3Points[C,1];
        assert (\x:RVKummerK3(C)=>rvKummerSwitch(rvKummerSwitch(x)))==(\x:RVKummerK3(C)=>x);
        assert (\x:RVKummerK3(C)=>rvKummerSwitch(rvKummerK3Translate[B,6](x)))
            ==(\x:RVKummerK3(C)=>rvKummerK3Translate[B,6](rvKummerSwitch(x)));
        ",
    );
    let branches = [
        ("z", "o", 0),
        ("o", "z", 1),
        ("o", "o", 4),
        ("o", "two", 7),
        ("o", "three", 13),
        ("o", "four", 15),
    ];
    for a in 0..16 {
        // The six branch parameters are the intersections E_(a xor b) ∩ T_a.
        for (u, v, b) in branches {
            let other = a ^ b;
            let x = format!("x_{a}_{b}");
            let y = format!("y_{a}_{b}");
            writeln!(source, "let {x}=rvKummerTropeParam[B,{a}]({u},{v});").unwrap();
            writeln!(source, "assert rvKummerOnTrope[B,{a}]({x})==S.hi;").unwrap();
            writeln!(
                source,
                "assert rvKummerBlowdown({x})==rvKummerTranslate[B,{other}](identity);"
            )
            .unwrap();
            writeln!(source, "let {y}=rvKummerSwitch({x});").unwrap();
            writeln!(
                source,
                "assert {y}==rvKummerTropeParam[B,{other}]({u},{v});"
            )
            .unwrap();
            writeln!(source, "assert rvKummerIsExceptional({y})==S.hi;").unwrap();
            writeln!(
                source,
                "assert rvKummerBlowdown({y})==rvKummerTranslate[B,{a}](identity);"
            )
            .unwrap();
            writeln!(source, "assert rvKummerOnTrope[B,{other}]({y})==S.hi;").unwrap();
        }
        for v in ["five", "six"] {
            let x = format!("x_{a}_{v}");
            let y = format!("y_{a}_{v}");
            writeln!(source, "let {x}=rvKummerTropeParam[B,{a}](o,{v});").unwrap();
            writeln!(source, "assert rvKummerIsExceptional({x})==S.lo;").unwrap();
            writeln!(source, "let {y}=rvKummerSwitch({x});").unwrap();
            writeln!(source, "assert rvKummerIsExceptional({y})==S.hi;").unwrap();
            writeln!(
                source,
                "assert rvKummerBlowdown({y})==rvKummerTranslate[B,{a}](identity);"
            )
            .unwrap();
            writeln!(source, "assert rvKummerSwitch({y})=={x};").unwrap();
        }
    }
    let output = stdout(interpret(&source));
    assert_eq!(
        output,
        "rvKummerK3Points[C, 1] = singular=48, nodes=16, resolved=160\n"
    );
}

#[test]
fn the_switch_requires_a_resolved_point_and_preserves_its_nominal_curve() {
    let cases = [
        ("show rvKummerSwitch(zero[RVKummer(C)]);", "an RVKummerK3 point"),
        (
            "rvhyperelliptic D[F]{z,o,z,o,o,o}; assert rvKummerSwitch(rvKummerSmooth(rvKummer[C](o,z,two,o)))==rvKummerSmooth(rvKummer[D](o,two,two,z));",
            "expected type RVKummerK3(C), found RVKummerK3(D)",
        ),
    ];
    for (input, expected) in cases {
        let error = diagnostic(interpret(&format!(
            "{F3}rvhyperelliptic C[F]{{z,o,z,o,o,o}};{input}"
        )));
        assert!(error.contains(expected), "{error}");
    }
}
