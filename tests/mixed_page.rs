mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn actual_odd_prime_page_lifts_and_differentials() {
    assert_eq!(
        stdout(example("mixed-pages")),
        concat!(
            "mixdim[Stair, 4, 0, 0] = 0\n",
            "MixedU(Stair, 0){u^0*x0 + u^1*x1 + u^2*x2}\n",
            "MixedE(Stair, 3, -1, 3){g0}\n",
        )
    );
    let src = include_str!("../examples/mixed-pages.asl");
    for (expression, diagnostic_text) in [
        ("mixNext(e3)", "requires a cycle for d_3"),
        ("mixPageClass[3,0](a)", "fails the page cycle condition"),
        ("mixPageClass[1,1](a)", "below the requested filtration"),
        ("mixPrev(e0)", "requires a positive page number"),
        ("mixU[0](e0)", "requires a MixedChain"),
    ] {
        let error = diagnostic(interpret(&format!("{src} show {expression};")));
        assert!(error.contains(diagnostic_text), "{expression}: {error}");
    }
}

#[test]
fn declaration_checks_signed_anticommuting_composites_and_prime_field_scalars() {
    let source = r"
        mixedcomplex Square[3] { a:0, b:1, c:-1, d:0 }
          { b => [(a,1)], d => [(c,1)] }
          { c => [(a,1)], d => [(b,-1)] };
        let d = mixchain[Square,0][(d,-2),(a,1),(a,2)];
        let a = mixchain[Square,0][(a,1)];
        assert mixb(mixB(d)) == pscale[2](a);
        assert mixB(mixb(d)) == a;
        assert add(mixb(mixB(d)),mixB(mixb(d))) == zero[MixedChain(Square,0)];
        assert mixD(mixD(mixU[3](d))) == zero[MixedU(Square,-8)];
        assert pscale[3](d) == zero[MixedChain(Square,0)];
        card \x : MixedChain(Square,0) => S.hi;
        card \x : MixedE(Square,0,0,0) => S.hi;
        show mixb(mixB(d));
        ";
    assert_eq!(
        stdout(interpret(source)),
        "card = 9 / 9\ncard = 9 / 9\nMixedChain(Square, 0){2*a}\n"
    );
    let wrong_sign = source.replace("d => [(b,-1)]", "d => [(b,1)]");
    assert!(diagnostic(interpret(&wrong_sign)).contains("bB + Bb is nonzero"));
}

#[test]
fn invalid_presentations_and_extreme_degrees_have_clear_semantics() {
    for (src, text) in [
        (
            "mixedcomplex P[3] {a:0,b:1,c:2}{b=>[(a,1)],c=>[(b,1)]}{};",
            "b squared",
        ),
        (
            "mixedcomplex P[3] {a:0,b:1,c:2}{}{a=>[(b,1)],b=>[(c,1)]};",
            "B squared",
        ),
        (
            "mixedcomplex P[3] {a:0,b:1}{a=>[(b,1)]}{};",
            "wrong homological degree",
        ),
        ("mixedcomplex P[9] {a:0}{}{};", "not prime"),
        ("mixedcomplex P[4] {a:0}{}{};", "requires a prime field"),
        (
            "mixedcomplex P[2] {a:0,a:0}{}{};",
            "duplicate mixed complex generator",
        ),
        (
            "mixedcomplex P[2] {a:0,b:-1}{a=>[(b,1)],a=>[]}{};",
            "duplicate mixed complex differential column",
        ),
        (
            "mixedcomplex P[5] {a:0}{z=>[]}{};",
            "unknown mixed complex generator",
        ),
    ] {
        let error = diagnostic(interpret(src));
        assert!(error.contains(text), "{src}: {error}");
    }
    let program = r"
        mixedcomplex E[2] { x:0, y:9223372036854775807 } {} {};
        let raw = mixchain[E,0][(x,1)];
        let z = mixU[4611686018427387904](raw);
        assert mixPageClass[31,4611686018427387904](z) == mixgen[E,31,-9223372036854775808,4611686018427387904,0];
        mixdim[E,100,9223372036854775807,0];
        mixdim[E,100,-9223372036854775808,9223372036854775807];
        ";
    assert_eq!(stdout(interpret(program)),
        "mixdim[E, 100, 9223372036854775807, 0] = 1\nmixdim[E, 100, -9223372036854775808, 9223372036854775807] = 0\n");
}

#[test]
fn the_declared_columns_produce_d20_over_multiple_prime_fields_with_a_control() {
    for (p, control) in [(2, false), (5, false), (7, true)] {
        let r = 20;
        let mut source = format!("mixedcomplex Stair[{p}] {{");
        let generators = (0..r)
            .map(|i| format!("x{i}:{},y{}:{}", 2 * i, i + 1, 2 * i + 1))
            .collect::<Vec<_>>();
        source.push_str(&generators.join(","));
        source.push_str("}{");
        let b = (1..r)
            .map(|i| format!("x{i} => [(y{i},-1)]"))
            .collect::<Vec<_>>();
        source.push_str(&b.join(","));
        source.push_str("}{");
        let big_b = (0..r - usize::from(control))
            .map(|i| format!("x{i} => [(y{},1)]", i + 1))
            .collect::<Vec<_>>();
        source.push_str(&big_b.join(","));
        source.push_str("};");
        for i in 0..r {
            source.push_str(&format!(
                "assert mixD(mixgen[Stair,{i},0,0,0]) == zero[MixedE(Stair,{i},-1,{i})];"
            ));
        }
        source.push_str(&format!("let s = mixgen[Stair,{r},0,0,0];"));
        if control {
            source.push_str(&format!(
                "assert mixNext(s) == mixgen[Stair,{},0,0,0];",
                r + 1
            ));
        } else {
            source.push_str(&format!("assert mixD(s) == mixgen[Stair,{r},-1,{r},0];"));
        }
        stdout(interpret(&source));
    }
}

#[test]
fn next_page_removes_the_incoming_and_outgoing_classes_with_field_coefficients() {
    let source = r"
        mixedcomplex C[5] { a:0, z:0, x:2, y1:1, y2:3 }
            { x=>[(y1,-1)] }
            { a=>[(y1,1)], x=>[(y2,3)] };
        let a = mixgen[C,2,0,0,0];
        let z = mixgen[C,2,0,0,1];
        let y = mixgen[C,2,-1,2,0];
        assert mixD(a) == pscale[3](y);
        assert mixD(z) == zero[MixedE(C,2,-1,2)];
        assert mixNext(pscale[4](z)) == pscale[4](mixgen[C,3,0,0,0]);
        assert mixPrev(mixgen[C,3,0,0,0]) == z;
        assert mixNext(y) == zero[MixedE(C,3,-1,2)];
        let pair = add(pscale[2](a),pscale[4](z));
        assert mixPageClass[2,0](mixPageRep(pair)) == pair;
        card \x: MixedE(C,2,0,0) => S.hi;

        mixedcomplex Large[4294967291] { x:0, y:1 } {} { x=>[(y,-1)] };
        assert mixD(mixgen[Large,1,0,0,0]) == pscale[4294967290](mixgen[Large,1,-1,1,0]);
        ";
    assert_eq!(stdout(interpret(source)), "card = 25 / 25\n");
}
