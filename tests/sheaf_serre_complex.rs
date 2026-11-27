mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn append_stalks(code: &mut String, complex: &str, tag: &str, degrees: usize, points: &[&str]) {
    for degree in 0..degrees {
        writeln!(code, "shcohomology {tag}H{degree}[{complex},{degree}];").unwrap();
        for point in points {
            let name = format!("{tag}H{degree}{point}");
            writeln!(
                code,
                "shpull {name}[at{point},{tag}H{degree}];\nshdim[{name},0];"
            )
            .unwrap();
        }
    }
}

fn dims(output: &str, tag: &str) -> Vec<usize> {
    output
        .lines()
        .filter(|line| line.starts_with(&format!("shdim[{tag}H")))
        .map(|line| line.rsplit_once(" = ").unwrap().1.parse().unwrap())
        .collect()
}

#[test]
fn serre_distinguishes_the_actual_diamond_k_invariant() {
    assert_eq!(
        stdout(example("sheaf-serre-complex")),
        concat!(
            "shserreshift[TN] = 2\n",
            "shdim[N1d, 0] = 0\n",
            "shdim[N2d, 0] = 0\n",
            "shdim[S1d, 0] = 1\n",
            "shdim[S2d, 0] = 1\n",
        )
    );
    let fixture = include_str!("../examples/sheaf-serre-complex.asl");
    for prime in [2, 3, 5] {
        let mut code = fixture.replace("Diamond,3]", &format!("Diamond,{prime}]"));
        for (point, name) in [("a", "ata"), ("b", "atb"), ("c", "atc"), ("d", "atd")] {
            writeln!(
                code,
                "let {name} = tab One -> Diamond {{pt=>Diamond.{point}}};"
            )
            .unwrap();
        }
        let points = ["a", "b", "c", "d"];
        append_stalks(&mut code, "Nonsplit", "IN", 2, &points);
        append_stalks(&mut code, "Split", "IS", 2, &points);
        append_stalks(&mut code, "TN", "ON", 4, &points);
        append_stalks(&mut code, "TS", "OS", 4, &points);
        code.push_str(
            r"
            assert shD(shTerm[ONH2,[a],0]) == add(shTerm[ONH2,[a,b],0,-1],shTerm[ONH2,[a,c],0,-1]);
            assert shD(shTerm[OSH2,[a],0]) == add(shTerm[OSH2,[a,d],0,-1],
                add(shTerm[OSH2,[a,b],0,-1],shTerm[OSH2,[a,c],0,-1]));
        ",
        );
        let output = stdout(interpret(&code));
        assert_eq!(dims(&output, "IN"), [0, 0, 0, 1, 1, 0, 0, 0]);
        assert_eq!(dims(&output, "IN"), dims(&output, "IS"));
        assert_eq!(
            dims(&output, "ON"),
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0]
        );
        assert_eq!(
            dims(&output, "OS"),
            [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0]
        );
    }
}

#[test]
fn a_signed_contractible_summand_disappears_after_the_serre_construction() {
    let fixture = include_str!("../examples/sheaf-serre-complex.asl");
    let mut code = fixture.replace("Diamond,3]", "Diamond,5]");
    code.push_str(
        r"
        shmap minus[Pa,Pa] {a:[-1],b:[-1],c:[-1],d:[-1]};
        shcomplex Contractible {Pa,minus,Pa};
        shserrecomplex TContractible[Contractible];
        sheaf Sum0[Diamond,5] {a:1,b:2,c:2,d:3}
            {a<b:[0,1],a<c:[0,1],b<d:[1,0,0,0,0,1],c<d:[0,0,1,0,0,1]};
        sheaf Sum1[Diamond,5] {a:2,b:2,c:2,d:2}
            {a<b:[1,0,0,1],a<c:[1,0,0,1],b<d:[1,0,0,1],c<d:[1,0,0,1]};
        shmap blockSum[Sum0,Sum1] {a:[0,-1],b:[1,0,0,-1],c:[1,0,0,-1],d:[1,1,0,0,0,-1]};
        shcomplex Stabilized {Sum0,blockSum,Sum1};
        shserrecomplex TStabilized[Stabilized];
        let ata = tab One -> Diamond {pt=>Diamond.a};
        let atb = tab One -> Diamond {pt=>Diamond.b};
        let atc = tab One -> Diamond {pt=>Diamond.c};
        let atd = tab One -> Diamond {pt=>Diamond.d};
    ",
    );
    let points = ["a", "b", "c", "d"];
    append_stalks(&mut code, "TContractible", "Acyclic", 4, &points);
    append_stalks(&mut code, "TStabilized", "Stable", 4, &points);
    code.push_str(
        r"
        let to = shExtGen[N2,StableH2,0,0];
        let from = shExtGen[StableH2,N2,0,0];
        show shExtCompose(from,to);
        show shExtCompose(to,from);
    ",
    );
    let output = stdout(interpret(&code));
    assert_eq!(dims(&output, "Acyclic"), [0; 16]);
    assert_eq!(
        dims(&output, "Stable"),
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0]
    );
    assert!(output
        .lines()
        .rev()
        .take(2)
        .all(|line| line.starts_with("ShExt(") && !line.ends_with("{0}")));
}

fn star_case(prime: usize, support: [usize; 4], arrows: [&str; 3]) -> String {
    let mut code = "space StarD4 {c,a,b,d,c<a,c<b,c<d};\n".to_string();
    let points = ["c", "a", "b", "d"];
    for point in points {
        writeln!(
            code,
            "let at{point} = tab One -> StarD4 {{pt=>StarD4.{point}}};"
        )
        .unwrap();
    }
    let [center, a, b, d] = support;
    let [ma, mb, md] = arrows;
    writeln!(code, "sheaf F[StarD4,{prime}] {{c:{center},a:{a},b:{b},d:{d}}} {{c<a:[{ma}],c<b:[{mb}],c<d:[{md}]}};").unwrap();
    code.push_str(
        r"
        shcomplex K {F};
        shserrecomplex T1[K];
        shserrecomplex T2[T1];
        shserrecomplex T3[T2];
        shserreshift[T1];
        shserreshift[T2];
        shserreshift[T3];
    ",
    );
    append_stalks(&mut code, "T3", "Return", 4, &points);
    code.push_str(
        r"
        shExtDim[F,F,0];
        shExtDim[ReturnH1,ReturnH1,0];
        shExtDim[F,ReturnH1,0];
        shExtDim[ReturnH1,F,0];
        let to = shExtGen[F,ReturnH1,0,0];
        let from = shExtGen[ReturnH1,F,0,0];
        show shExtCompose(from,to);
        show shExtCompose(to,from);
    ",
    );
    code
}

#[test]
fn d4_half_period_returns_simples_projectives_and_the_highest_root() {
    for prime in [2, 3, 5] {
        for (support, arrows) in [
            ([1, 0, 0, 0], ["", "", ""]),
            ([0, 1, 0, 0], ["", "", ""]),
            ([0, 0, 1, 0], ["", "", ""]),
            ([0, 0, 0, 1], ["", "", ""]),
            ([1, 1, 1, 1], ["1", "1", "1"]),
            ([2, 1, 1, 1], ["1,0", "0,1", "1,1"]),
        ] {
            let output = stdout(interpret(&star_case(prime, support, arrows)));
            let expected = [[0; 4], support, [0; 4], [0; 4]].concat();
            assert_eq!(dims(&output, "Return"), expected);
            assert!(output.lines().take(3).all(|line| line.ends_with(" = 1")));
            let end: Vec<_> = output.lines().rev().take(6).collect();
            assert!(
                end[..2]
                    .iter()
                    .all(|line| line.starts_with("ShExt(") && !line.ends_with("{0}")),
                "{output}"
            );
            assert!(
                end[2..]
                    .iter()
                    .all(|line| line.starts_with("shExtDim[") && line.ends_with(" = 1")),
                "{output}"
            );
        }
    }
}

#[test]
fn serre_complex_requires_an_ordinary_declared_complex() {
    let source = r"
        sheaf F[One,5] {pt:1} {};
        shcomplex K {F};
    ";
    for (statement, expected) in [
        ("shserrecomplex T[F];", "unknown sheaf complex 'F'"),
        ("shserrecomplex K[K];", "already defined"),
        ("shserre T[K];", "unknown sheaf 'K'"),
    ] {
        assert!(diagnostic(interpret(&format!("{source}{statement}"))).contains(expected));
    }
}

#[test]
fn the_internal_arrow_has_the_koszul_sign_in_odd_bar_length() {
    let program = r"
        space X {a<b};
        sheaf F[X,5] {a:1,b:1} {a<b:[1]};
        shmap minus[F,F] {a:[-1],b:[-1]};
        shcomplex K {F,minus,F};
        shserrecomplex T[K];
        let u = shTotEmbed[T,0](shTerm[T__0,[a],0]);
        let v = shTotEmbed[T,0](shTerm[T__0,[b],0]);
        assert shTotPart[1](shD(u)) == add(shTerm[T__1,[a],0,-1],
            add(shTerm[T__1,[a],1],shTerm[T__1,[a],2]));
        assert shTotPart[1](shD(v)) == add(shTerm[T__1,[b],0],shTerm[T__1,[b],1]);
        assert shD(shD(u)) == zero[ShTot(T,2)];
        assert shD(shD(v)) == zero[ShTot(T,2)];
    ";
    assert_eq!(stdout(interpret(program)), "");
}
