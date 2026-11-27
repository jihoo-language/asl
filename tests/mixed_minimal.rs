mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn the_odd_prime_staircase_transfers_to_a_single_higher_operator() {
    assert_eq!(
        stdout(example("mixed-minimal")),
        concat!(
            "mixmindim[Min, 0] = 1\n",
            "mixmindim[Min, 2] = 0\n",
            "mixmindim[Min, 5] = 1\n",
            "MixedMinU(Min, -1){u^3*g1}\n",
            "MixedU(Stair, 0){u^0*x0 + u^1*x1 + u^2*x2}\n",
        )
    );
    let split = include_str!("../examples/mixed-minimal.asl");
    let split = split
        .split_once("assert mixDelta[3]")
        .unwrap()
        .0
        .replace(", x2 => [(y3,1)]", "");
    let source = format!(
        "{split}
        assert mixDelta[3](x) == zero[MixedMin(Min,5)];
        assert mixD(X) == zero[MixedMinU(Min,-1)];
        assert mixMinI(X) == lift;
        assert mixD(lift) == zero[MixedU(Stair,-1)];
        assert mixMinP[Min](lift) == X;
        "
    );
    stdout(interpret(&source));
}

#[test]
fn the_contraction_and_perturbation_hold_on_every_scrambled_generator() {
    let fixture = include_str!("../examples/mixed-bars.asl");
    let generators = fixture
        .split_once('{')
        .unwrap()
        .1
        .split_once('}')
        .unwrap()
        .0;
    let mut source = format!("{fixture} mixminimal N[M];");
    for (i, part) in generators.split(',').enumerate() {
        let (name, degree) = part.split_once(':').unwrap();
        let name = name.trim();
        let b: i64 = degree.trim().parse().unwrap();
        source.push_str(&format!(
            "let raw{i} = mixchain[M,{b}][({name},1)];
             let term{i} = mixU[2](raw{i});
             let h{i} = mixMinH[N](raw{i});
             let H{i} = mixMinH[N](term{i});
             let p{i} = mixMinP[N](raw{i});
             let P{i} = mixMinP[N](term{i});
             assert add(mixb(h{i}),mixMinH[N](mixb(raw{i}))) == add(raw{i},pscale[2](mixMinI(p{i})));
             assert mixMinP[N](mixMinI(p{i})) == p{i};
             assert mixMinP[N](h{i}) == zero[MixedMin(N,{up})];
             assert mixMinH[N](h{i}) == zero[MixedChain(M,{up2})];

             assert add(mixD(H{i}),mixMinH[N](mixD(term{i}))) == add(term{i},pscale[2](mixMinI(P{i})));
             assert mixMinP[N](mixMinI(P{i})) == P{i};
             assert mixMinP[N](H{i}) == zero[MixedMinU(N,{total_up})];
             assert mixMinH[N](H{i}) == zero[MixedU(M,{total_up2})];
             assert mixMinH[N](mixMinI(P{i})) == zero[MixedU(M,{total_up})];
             assert mixD(P{i}) == mixMinP[N](mixD(term{i}));
             assert mixMinI(mixD(P{i})) == mixD(mixMinI(P{i}));
             assert mixD(mixD(P{i})) == zero[MixedMinU(N,{total_down2})];
             assert mixMinP[N](mixU[3](term{i})) == mixU[3](P{i});
            ",
            up=b+1,up2=b+2,total_up=b-3,total_up2=b-2,total_down2=b-6,
        ));
    }
    stdout(interpret(&source));
}

#[test]
fn minimal_values_support_finite_fields_and_validate_source_and_grades() {
    let fixture = include_str!("../examples/mixed-minimal.asl");
    let extra = r"
        mixedcomplex T[3] {a:-2,b:-2} {} {};
        mixminimal H[T];
        card \x: MixedMin(H,-2) => S.hi;
        card \x: MixedMinU(H,-4) => S.hi;
    ";
    assert!(
        stdout(interpret(&format!("{fixture}{extra}"))).ends_with("card = 9 / 9\ncard = 9 / 9\n")
    );
    for (extra, expected) in [
        ("show mixMinGen[Min,2,0];", "generator index"),
        ("show mixMinI(z);", "a MixedMin or MixedMinU value"),
        ("show mixDelta[1](X);", "a MixedMin class"),
        ("show mixD(x);", "a MixedU, MixedE, or MixedMinU value"),
        (
            "show mixDelta[9223372036854775807](x);",
            "degree is out of range",
        ),
        (
            "mixedcomplex C[3] {a:0} {} {}; show mixMinP[Min](mixchain[C,0][(a,1)]);",
            "different source complex",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{fixture} {extra}")));
        assert!(error.contains(expected), "{extra}: {error}");
    }
}
