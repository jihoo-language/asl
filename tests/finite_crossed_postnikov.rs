mod common;

use std::collections::BTreeMap;
use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn list(values: impl IntoIterator<Item = usize>) -> String {
    let entries = values
        .into_iter()
        .map(|v| v.to_string())
        .collect::<Vec<_>>()
        .join(",");
    format!("{{{entries}}}")
}

fn table(rows: impl IntoIterator<Item = impl IntoIterator<Item = usize>>) -> String {
    let entries = rows.into_iter().map(list).collect::<Vec<_>>().join(",");
    format!("{{{entries}}}")
}

fn cyclic_crossed(p: usize, active: usize) -> String {
    let n = p * p;
    let action = table((0..n).map(|g| (0..n).map(move |h| ((1 + p * active * g) * h) % n)));
    let boundary = list((0..n).map(|h| p * h % n));
    format!(
        "rvgroup H[{n}];rvgroup E[{n}];rvcrossed C[H,E]{boundary}{action};rvpostnikov P[C,{p}];\n"
    )
}

#[test]
fn example_shows_nonsplit_and_genuinely_twisted_skeleta() {
    assert_eq!(stdout(example("finite-crossed-postnikov")), concat!(
        "rvCrossedData[C] = kernel=[0, 2], quotient_representatives=[0, 1]\n",
        "rvPostnikovData[P] = quotient=2, prime=2, action=[1, 1], rank_d2=0, rank_d3=0, dim_h3=1, nonzero=true\n",
        "rvPostnikovClass = nonzero=true, reduced=[(0, 1)]\n",
        "rvPostnikovCoeff[1, 1, 1] = 1\n",
        "rvPostnikovData[Q] = quotient=2, prime=3, action=[1, 2], rank_d2=1, rank_d3=0, dim_h3=0, nonzero=false\n",
        "rvPostnikovCoeff[1, 1, 1] = 1\n",
        "rvPostnikovClass = nonzero=false, reduced=[]\n",
    ));
}

fn integral_cycle(p: usize) -> Vec<[usize; 3]> {
    let mut boundary = BTreeMap::<[usize; 2], i32>::new();
    for i in 1..p {
        for (pair, sign) in [
            ([i, 1], 1),
            ([(i + 1) % p, 1], -1),
            ([1, (i + 1) % p], 1),
            ([1, i], -1),
        ] {
            if !pair.contains(&0) {
                *boundary.entry(pair).or_default() += sign;
            }
        }
    }
    assert!(boundary.values().all(|&x| x == 0));
    (0..p).map(|i| [1, i, 1]).collect()
}

fn triples(p: usize) -> impl Iterator<Item = [usize; 3]> {
    (1..p).flat_map(move |x| (1..p).flat_map(move |y| (1..p).map(move |z| [x, y, z])))
}

#[test]
fn cyclic_primes_match_the_exact_class_cycle_and_changed_section() {
    for (p, rank_d2, rank_d3) in [(2, 0, 0), (3, 2, 5), (5, 12, 51)] {
        let n = p * p;
        for active in 0..=1 {
            let mut source = cyclic_crossed(p, active);
            let expected = list(triples(p).map(|[x, y, z]| (active * x * ((y + z) / p)) % p));
            writeln!(
                source,
                "let k=rvPostnikovK[P]; assert k==rvPostnikovCochain[P,3]{expected};"
            )
            .unwrap();
            let section = |x| (x + p * x * x) % n;
            let lift = |x, y| {
                (section(x) + section(y) + n - section((x + y) % p)) % n / p
                    + p * ((x * y * (x + 1)) % p)
            };
            let lifts = table((0..p).map(|x| (0..p).map(move |y| lift(x, y))));
            let section_values = list((0..p).map(section));
            let alternate = list(triples(p).map(|[x, y, z]| {
                let left = (1 + p * active * section(x)) * lift(y, z) + lift(x, (y + z) % p);
                let right = lift((x + y) % p, z) + lift(x, y);
                let image = (left + 2 * n - right) % n;
                assert_eq!(image % p, 0);
                image / p
            }));
            writeln!(source, "let alt=rvPostnikovLifted[P]{section_values}{lifts}; assert alt==rvPostnikovCochain[P,3]{alternate};").unwrap();
            source.push_str("let witness=rvPostnikovGaugeWitness(k,alt); assert rvPostnikovGauge(k,witness)==alt; assert rvPostnikovCoboundary(witness)==add(alt,neg(k)); rvPostnikovData[P]; rvPostnikovClass(k); rvPostnikovClass(alt);\n");
            let b = |x: usize, y: usize| if x == 0 || y == 0 { 0 } else { (x * x + y) % p };
            let gauge = list((1..p).flat_map(|x| (1..p).map(move |y| b(x, y))));
            let coboundary = list(triples(p).map(|[x, y, z]| {
                (b(y, z) + b(x, (y + z) % p) + 2 * p - b((x + y) % p, z) - b(x, y)) % p
            }));
            writeln!(source, "let b=rvPostnikovCochain[P,2]{gauge}; let delta=rvPostnikovCoboundary(b); assert delta==rvPostnikovCochain[P,3]{coboundary}; let changed=rvPostnikovGauge(k,b);").unwrap();
            for [x, y, z] in integral_cycle(p) {
                writeln!(source,"rvPostnikovCoeff[{x},{y},{z}](k); rvPostnikovCoeff[{x},{y},{z}](alt); rvPostnikovCoeff[{x},{y},{z}](delta);").unwrap();
            }
            let actual = stdout(interpret(&source));
            let mut lines = actual.lines();
            let data = lines.next().unwrap();
            assert!(
                data.ends_with(&format!(
                    "rank_d2={rank_d2}, rank_d3={rank_d3}, dim_h3=1, nonzero={}",
                    active != 0
                )),
                "{data}"
            );
            assert_eq!(lines.next(), lines.next());
            let mut evaluations = [0; 3];
            for _ in 0..p {
                for entry in &mut evaluations {
                    *entry += lines
                        .next()
                        .unwrap()
                        .split_once(" = ")
                        .unwrap()
                        .1
                        .parse::<usize>()
                        .unwrap();
                }
            }
            assert_eq!(evaluations.map(|x| x % p), [active, active, 0]);
            assert!(lines.next().is_none());
        }
    }
}

fn pentagon(
    source: &mut String,
    n: usize,
    product: impl Fn(usize, usize) -> usize,
    name: &str,
    cochain: &str,
) {
    for g in 0..n {
        for h in 0..n {
            for k in 0..n {
                for l in 0..n {
                    let (gh, hk, kl) = (product(g, h), product(h, k), product(k, l));
                    writeln!(source,r"
                        assert rvPostnikovVertical(rvPostnikovAssoc[{g},{h},{kl}]({cochain}),rvPostnikovAssoc[{gh},{k},{l}]({cochain}))
                        ==rvPostnikovVertical(rvPostnikovHorizontal(rvPostnikovCell[{name},{g},0],rvPostnikovAssoc[{h},{k},{l}]({cochain})),
                            rvPostnikovVertical(rvPostnikovAssoc[{g},{hk},{l}]({cochain}),rvPostnikovHorizontal(rvPostnikovAssoc[{g},{h},{k}]({cochain}),rvPostnikovCell[{name},{l},0])));
                        ").unwrap();
                }
            }
        }
    }
}

#[test]
fn all_typed_pentagons_hold_on_nonsplit_and_twisted_cochains() {
    let mut cyclic = cyclic_crossed(3, 1);
    cyclic.push_str("let k=rvPostnikovK[P]; let b=rvPostnikovCochain[P,2]{1,0,0,2}; let changed=rvPostnikovGauge(k,b);\n");
    pentagon(&mut cyclic, 3, |x, y| (x + y) % 3, "P", "k");
    pentagon(&mut cyclic, 3, |x, y| (x + y) % 3, "P", "changed");
    assert_eq!(stdout(interpret(&cyclic)), "");

    let mut twist = r"
        rvgroup H[3];rvgroup E[2];
        rvcrossed C[H,E]{0,0,0}{{0,1,2},{0,2,1}};rvpostnikov P[C,1];
        let k=rvPostnikovGauge(rvPostnikovK[P],rvPostnikovCochain[P,2]{1});
        assert rvCrossedKernelAct[C,1](rvGroupElem[H,1])==rvGroupElem[H,2];
        assert rvPostnikovHorizontal(rvPostnikovCell[P,1,0],rvPostnikovCell[P,0,1])==rvPostnikovCell[P,1,2];
        assert rvPostnikovVertical(rvPostnikovCell[P,1,2],rvPostnikovAssoc[1,1,1](k))==rvPostnikovCell[P,1,0];
        let a=rvPostnikovCell[P,1,1]; let b=rvPostnikovCell[P,0,1];
        let c=rvPostnikovCell[P,1,2]; let d=rvPostnikovCell[P,0,1];
        assert rvPostnikovVertical(rvPostnikovHorizontal(a,b),rvPostnikovHorizontal(c,d))
            ==rvPostnikovHorizontal(rvPostnikovVertical(a,c),rvPostnikovVertical(b,d));
        assert (\x : RVPostnikovCochain(P,2) => add(x,neg(x)))
            ==(\x : RVPostnikovCochain(P,2) => zero[RVPostnikovCochain(P,2)]);
        ".to_owned();
    pentagon(&mut twist, 2, |x, y| (x + y) % 2, "P", "k");
    assert_eq!(stdout(interpret(&twist)), "");
}

fn s3_mul(x: usize, y: usize) -> usize {
    let (a, b, c, d) = (x % 3, x / 3, y % 3, y / 3);
    (a + if b == 0 { c } else { 3 - c }) % 3 + 3 * ((b + d) % 2)
}

fn s3() -> String {
    format!(
        "rvgroup S3{};\n",
        table((0..6).map(|x| (0..6).map(move |y| s3_mul(x, y))))
    )
}

#[test]
fn nonabelian_tables_peiffer_action_and_bar_quotient_are_genuine() {
    let mut source = s3();
    let conjugations = table((0..6).map(|g| {
        let inv = (0..6).find(|&h| s3_mul(g, h) == 0).unwrap();
        (0..6).map(move |x| s3_mul(s3_mul(g, x), inv))
    }));
    writeln!(
        source,
        "rvcrossed Inner[S3,S3]{{0,1,2,3,4,5}}{conjugations};"
    )
    .unwrap();
    source.push_str(
        r"
        assert rvGroupMul(rvGroupElem[S3,3],rvGroupElem[S3,1])==rvGroupElem[S3,5];
        assert rvGroupMul(rvGroupElem[S3,1],rvGroupElem[S3,3])==rvGroupElem[S3,4];
        assert rvGroupInv(rvGroupElem[S3,1])==rvGroupElem[S3,2];
        assert (\x : RVGroup(S3) => rvGroupMul(x,rvGroupElem[S3,0]))==(\x : RVGroup(S3) => x);
        assert rvCrossedAct[Inner](rvGroupElem[S3,3],rvGroupElem[S3,1])==rvGroupElem[S3,2];
        assert rvCrossedAssoc[Inner,0,0,0]==rvGroupElem[S3,0];
        rvCrossedData[Inner];
        rvgroup Z[2];
        rvcrossed Product[Z,S3]{0,0}{{0,1},{0,1},{0,1},{0,1},{0,1},{0,1}};
        rvpostnikov P[Product,1];
        rvCrossedQuotientMul[Product,3,1];
        rvPostnikovData[P];
        ",
    );
    let cocycle = list(triples(6).map(|[x, y, z]| (x / 3) * (y / 3) * (z / 3)));
    writeln!(source,"let pullback=rvPostnikovCochain[P,3]{cocycle}; rvPostnikovClass(pullback); rvPostnikovCoeff[3,3,3](pullback);").unwrap();
    let actual = stdout(interpret(&source));
    let mut lines = actual.lines();
    assert_eq!(
        lines.next().unwrap(),
        "rvCrossedData[Inner] = kernel=[0], quotient_representatives=[0]"
    );
    assert_eq!(
        lines.next().unwrap(),
        "rvCrossedQuotientMul[Product,3,1] = 5"
    );
    assert!(lines.next().unwrap().ends_with("dim_h3=1, nonzero=false"));
    assert!(lines
        .next()
        .unwrap()
        .starts_with("rvPostnikovClass = nonzero=true"));
    assert_eq!(lines.next().unwrap(), "rvPostnikovCoeff[3, 3, 3] = 1");
    assert!(lines.next().is_none());
    let wrong = diagnostic(interpret(&format!(
        "{source}show rvPostnikovGaugeWitness(rvPostnikovK[P],pullback);"
    )));
    assert!(wrong.contains("not gauge equivalent"));
}

#[test]
fn malformed_tables_crossed_rules_choices_and_types_are_rejected() {
    let base = cyclic_crossed(3, 1);
    let s3 = s3();
    let trivial_s3 = table((0..6).map(|_| 0..6));
    let oversized_bar = format!(
        "rvgroup H[2];rvgroup E[25];rvcrossed C[H,E]{{0,0}}{};rvpostnikov P[C,1];",
        table((0..25).map(|_| [0, 1]))
    );
    for (source, expected) in [
        ("rvgroup Bad[129];".into(), "order must be between 1 and 128"),
        (oversized_bar, "degree-four size limit"),
        ("rvgroup H[4];rvgroup E[1];rvcrossed C[H,E]{0,0,0,0}{{0,1,2,3}};rvpostnikov P[C,1];".into(), "cyclic kernel of prime order"),
        ("rvgroup Bad{{0,1},{1}};".into(), "square table"),
        ("rvgroup Bad{{1,0},{0,1}};".into(), "identity 0"),
        ("rvgroup Bad{{0,1,2},{1,0,2},{2,1,0}};".into(), "is not associative"),
        ("rvgroup H[3];rvgroup E[2];rvcrossed C[H,E]{0,1,0}{{0,1,2},{0,1,2}};".into(), "boundary must be a group homomorphism"),
        ("rvgroup H[3];rvgroup E[2];rvcrossed C[H,E]{0,0,0}{{0,1,2},{0,2,2}};".into(), "action rows must be automorphisms"),
        (format!("{s3}rvcrossed C[S3,S3]{{0,1,2,3,4,5}}{trivial_s3};"), "first Peiffer identity"),
        (format!("{s3}rvgroup Z[1];rvcrossed C[S3,Z]{{0,0,0,0,0,0}}{{{{0,1,2,3,4,5}}}};"), "second Peiffer identity"),
        (format!("{base}rvpostnikov Q[C,1];"), "nonidentity kernel element"),
        (format!("{base}show rvPostnikovLifted[P]{{0,2,1}}{{{{0,0,0}},{{0,0,0}},{{0,0,0}}}};"), "normalized representative"),
        (format!("{base}show rvPostnikovLifted[P]{{0,1,2}}{{{{0,0,0}},{{0,0,0}}}};"), "lifts must form a square quotient table"),
        (format!("{base}show rvPostnikovLifted[P]{{0,1,2}}{{{{0,0,0}},{{0,0,0}},{{0,0,0}}}};"), "Postnikov lift at"),
        (format!("{base}show rvGroupMul(rvGroupElem[H,0],rvGroupElem[E,0]);"), "expected type RVGroup(H)"),
        (format!("{base}show rvPostnikovVertical(rvPostnikovCell[P,1,0],rvPostnikovCell[P,2,0]);"), "requires the same object"),
        (format!("{base}rvpostnikov Q[C,3]; show rvPostnikovGauge(rvPostnikovK[P],zero[RVPostnikovCochain(Q,2)]);"), "expected type RVPostnikovCochain(P,2)"),
        (format!("{base}rvPostnikovClass(rvPostnikovCochain[P,3]{{1,0,0,0,0,0,0,0}});"), "fails the twisted cocycle"),
    ] {
        let actual = diagnostic(interpret(&source));
        assert!(actual.contains(expected), "{source}\n{actual}");
    }
}
