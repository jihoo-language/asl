#[allow(dead_code)]
mod common;

use std::fmt::Write;
use std::process::Command;

use common::{diagnostic, interpret, stdout};

fn generated(left: usize, right: usize, prime: u32, explicit: bool) -> String {
    stdout(
        Command::new("awk")
            .current_dir(env!("CARGO_MANIFEST_DIR"))
            .args([
                "-v",
                &format!("left={left}"),
                "-v",
                &format!("right={right}"),
                "-v",
                &format!("prime={prime}"),
                "-v",
                if explicit { "explicit=1" } else { "explicit=0" },
                "-f",
                "examples/sheaf-higher-yoneda.awk",
            ])
            .output()
            .unwrap(),
    )
}

#[test]
fn two_unequal_arms_detect_the_same_relation_with_the_correct_parity() {
    for (left, right, p) in [(3, 4, 2), (3, 4, 3), (3, 5, 5), (4, 6, 7)] {
        let mut source = generated(left, right, p, true);
        let arity = (1..=left)
            .map(|i| format!("a{i}"))
            .collect::<Vec<_>>()
            .join(",");
        writeln!(source, "assert shInMassey({arity},resulta) == S.hi;").unwrap();
        writeln!(
            source,
            "assert shInMassey({arity},zero[ShExt(Fbottom,Ftop,2)]) == S.lo;"
        )
        .unwrap();
        writeln!(source, "let w = shMasseyWitness[0,2]({arity});").unwrap();
        source.push_str("assert shD(w) == shScale[-1](shExtCompose(shRep(a2),shRep(a1)));\n");
        let output = stdout(interpret(&source));
        let lines: Vec<_> = output.lines().collect();
        assert_eq!(lines[0], "shExtDim[Fbottom, Ftop, 2] = 1");
        assert_eq!(lines.len(), 3);
        assert!(!lines[1].ends_with("{0}"), "{output}");
        assert!(!lines[2].ends_with("{0}"), "{output}");
        assert_eq!(
            lines[1] == lines[2],
            p == 2 || (left + right) % 2 == 1,
            "{output}"
        );
    }
}

#[test]
fn explicit_curvature_and_exhaustive_membership_agree_on_nontrivial_indeterminacy() {
    let mut source = generated(2, 3, 5, true);
    source.push_str(
        r"
        let z0 = zero[ShExt(Fbottom,Fbottom,1)];
        let z1 = zero[ShExt(Fbottom,Fa1,1)];
        let q = shScale[3](resulta);
        let aw = shMasseyFind[0,2](z0,z1,a2,q);
        let bw = shMasseyFind[1,3](z0,z1,a2,q);
        assert aw == shScale[3](shRep(a1));
        assert bw == zero[ShHomC(Fbottom,Ftop,1)];
        assert shMassey(z0,z1,a2){0<2=>aw,1<3=>bw} == q;
        let raw = shCurvature(z0,z1,shRep(a2)){0<2=>aw};
        assert shD(raw) == zero[ShHomC(Fbottom,Ftop,3)];
        assert shClass(raw) == q;
        assert shMasseySolve(z0,z1,a2) == zero[ShExt(Fbottom,Ftop,2)];
        let contains = \v : ShExt(Fbottom,Ftop,2) => shInMassey(z0,z1,a2,v);
        card contains;
        let z2 = zero[ShExt(Ftop,Ftop,1)];
        assert shInMassey(z0,z1,a2,z2,q) == S.lo;
        assert shInMassey(z0,z1,a2,z2,zero[ShExt(Fbottom,Ftop,2)]) == S.hi;
        assert shInMassey(a1,a2,z2,zero[ShExt(Fbottom,Ftop,2)]) == S.lo;
    ",
    );
    let output = stdout(interpret(&source));
    assert!(output.ends_with("card = 5 / 5\n"), "{output}");
}

#[test]
fn defining_systems_reject_wrong_equations_intervals_and_types() {
    let source = generated(3, 4, 3, false);
    for (statement, reason) in [
        ("show shMassey(a1,a2,a3){};", "equation fails on interval [0,2]"),
        ("show shCurvature(a1,a2,a3){0<1=>shRep(a1)};", "invalid proper interval [0,1]"),
        ("show shCurvature(a1,a2,a3){0<3=>shRep(a1)};", "invalid proper interval [0,3]"),
        ("show shMassey(a1,a2,a3){0<2=>shMasseyWitness[0,2](a1,a2,a3),0<2=>shMasseyWitness[0,2](a1,a2,a3)};", "duplicate witness"),
        ("show shMassey(a1,a2,a3){0<2=>shRep(a1)};", "witness for interval [0,2]"),
        ("show shMasseySolve(a2,a1);", "requires composable inputs"),
        ("show shMasseySolve(resulta,zero[ShExt(Ftop,Ftop,1)]);", "must have degree one"),
        ("show shMasseySolve(shHomTerm[Fbottom,Ftop,[bottom,top],0,0],zero[ShExt(Ftop,Ftop,1)]);", "adjacent sheaf defining-system input must be closed"),
        ("show shMasseyFind[0,2](a1,a2,a3,zero[ShExt(Fbottom,Ftop,2)]);", "no sheaf defining system realizes the requested class"),
        ("show shInMassey(a1,a2,a3,a1);", "shInMassey: expected type ShExt(Fbottom, Ftop, 2)"),
    ] {
        let error = diagnostic(interpret(&format!("{source}{statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
    let source = generated(2, 4, 3, false);
    let error = diagnostic(interpret(&format!(
        "{source}show shMasseySolve(a1,a2,zero[ShExt(Ftop,Ftop,1)]);"
    )));
    assert!(error.contains("no sheaf defining system exists"), "{error}");
}

#[test]
fn exact_changes_of_the_defining_cochain_preserve_the_top_class() {
    let mut source = generated(2, 3, 5, false);
    source.push_str(
        r"
        sheaf E[Arms,5] {bottom:2,a1:1,b1:1,b2:1,top:1}
            {bottom<a1:[0,1],a1<top:[1],bottom<b1:[0,1],b1<b2:[1],b2<top:[1]};
        let z0 = zero[ShExt(E,E,1)];
        let z1 = zero[ShExt(E,Fa1,1)];
        let x = shHomTerm[E,Fa1,[bottom,a1],0,0];
        let beta = shHomTerm[E,Fa1,[a1],0,0];
        let ordinary = shCurvature(z0,z1,a2){0<2=>x};
        let shifted = shCurvature(z0,z1,a2){0<2=>add(x,shD(beta))};
        let difference = shExtCompose(shRep(a2),shD(beta));
        assert difference == shHomTerm[E,Ftop,[bottom,a1,top],1,0,-1];
        assert difference == shScale[-1](shD(shExtCompose(shRep(a2),beta)));
        assert shifted == add(ordinary,difference);
        assert shClass(shifted) == shClass(ordinary);
        assert shInMassey(z0,z1,a2,shClass(shifted)) == S.hi;
    ",
    );
    assert_eq!(stdout(interpret(&source)).lines().count(), 3);
}

#[test]
fn actual_arity_sixty_four_fits_the_default_budgets() {
    let program = generated(64, 35, 5, true);
    assert_eq!(
        stdout(interpret(&program)),
        concat!(
            "shExtDim[Fbottom, Ftop, 2] = 1\n",
            "ShExt(Fbottom, Ftop, 2){4*g0}\n",
            "ShExt(Fbottom, Ftop, 2){4*g0}\n",
        )
    );
}
