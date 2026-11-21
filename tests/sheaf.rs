mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn circle_sheaves_compute_constants_monodromy_and_zero_stalks() {
    assert_eq!(
        stdout(example("sheaf-cohomology")),
        concat!(
            "shdim[Constant, 0] = 1\n",
            "shdim[Constant, 1] = 1\n",
            "shdim[Unipotent, 0] = 1\n",
            "shdim[Unipotent, 1] = 1\n",
            "shdim[Twisted, 0] = 0\n",
            "shdim[Twisted, 1] = 0\n",
            "shdim[Edges, 0] = 0\n",
            "shdim[Edges, 1] = 3\n",
            "ShC(Edges, 1){4*[a<x]#0}\n",
        )
    );
}

#[test]
fn roos_differential_gets_rectangular_transitive_maps_and_signs_right() {
    let program = r"
        space Interval {c,b,a,a<b<c};
        sheaf F[Interval,5] {a:1,b:2,c:3} {
            a<b:[1,-2],
            b<c:[2,0,1,1,0,2]
        };
        let from_a = shTerm[F,[a],0];
        let to_ab = add(shTerm[F,[a,b],0,-1],shTerm[F,[a,b],1,2]);
        let to_ac = add(shTerm[F,[a,c],0,3],add(shTerm[F,[a,c],1],shTerm[F,[a,c],2,4]));
        assert shD(from_a) == add(to_ab,to_ac);
        assert shD(shTerm[F,[a,b],1,-1]) ==
            add(shTerm[F,[a,b,c],1,4],shTerm[F,[a,b,c],2,3]);
        assert shD(shTerm[F,[a,c],0]) == shTerm[F,[a,b,c],0,-1];
        assert shD(shTerm[F,[b,c],2]) == shTerm[F,[a,b,c],2];
        assert shD(shD(from_a)) == zero[ShC(F,2)];
        assert shClass(shD(from_a)) == zero[ShH(F,1)];
        shdim[F,0];
        shdim[F,1];
        shdim[F,2];
        shdim[F,3];
    ";
    assert_eq!(
        stdout(interpret(program)),
        "shdim[F, 0] = 1\nshdim[F, 1] = 0\nshdim[F, 2] = 0\nshdim[F, 3] = 0\n"
    );
}

#[test]
fn cochains_and_classes_are_first_class_prime_field_values() {
    let program = r"
        space Isolated {a,b};
        sheaf F[Isolated,3] {a:1,b:1} {};
        let from = \v : ShC(F,0) => shClass(v);
        let lift = \h : ShH(F,0) => shRep(h);
        let a = shTerm[F,[a],0,-1];
        let b = shTerm[F,[b],0,-9223372036854775808];
        let c = add(a,b);
        assert lift(from(c)) == c;
        assert from(lift(shgen[F,0,1])) == shgen[F,0,1];
        assert shScale[-1](a) == shTerm[F,[a],0];
        assert shScale[-4](from(c)) == from(shScale[-1](c));
        assert add(shScale[-1](from(c)),from(c)) == zero[ShH(F,0)];
        assert shScale[0](from(c)) == zero[ShH(F,0)];
        assert shScale[3](c) == zero[ShC(F,0)];
        card (\v : ShH(F,0) => S.hi);
        card (\v : ShC(F,0) => S.hi);
        show c;
        show from(c);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "card = 9 / 9\ncard = 9 / 9\nShC(F, 0){2*[a]#0 + [b]#0}\nShH(F, 0){2*g0 + g1}\n"
    );
}

#[test]
fn empty_posets_large_stalks_and_high_degrees_do_not_force_raw_enumeration() {
    let program = r"
        space Empty {};
        sheaf E[Empty,2] {} {};
        sheaf Z[One,7] {pt:0} {};
        sheaf L[One,5] {pt:1000000000} {};
        show shTerm[L,[pt],999999999,-1];
        assert shD(shTerm[L,[pt],999999999]) == zero[ShC(L,1)];
        assert shClass(zero[ShC(L,0)]) == zero[ShH(L,0)];
        shdim[E,0];
        shdim[Z,0];
        shdim[L,1000000000];
        card (\v : ShH(E,0) => S.hi);
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "ShC(L, 0){4*[pt]#999999999}\n",
            "shdim[E, 0] = 0\n",
            "shdim[Z, 0] = 0\n",
            "shdim[L, 1000000000] = 0\n",
            "card = 1 / 1\n",
        )
    );
    let too_many = diagnostic(interpret(&format!("{program}\nshdim[L,0];")));
    assert!(too_many.contains("cochain rank exceeds"), "{too_many}");
}

#[test]
fn validate_hasse_maps_their_shapes_and_path_independence() {
    let diamond = "space D { d,a,b,c,a<b,a<c,b<d,c<d };";
    let stalks = "{a:1,b:1,c:1,d:1}";
    for (statement, reason) in [
        (
            format!("sheaf Bad[D,3] {stalks} {{a<b:[1],a<c:[1],b<d:[1],c<d:[2]}};"),
            "restriction paths from 'a' to 'd' do not commute",
        ),
        (
            format!("sheaf Bad[D,3] {stalks} {{a<d:[1]}};"),
            "is not a Hasse cover",
        ),
        (
            format!("sheaf Bad[D,3] {stalks} {{a<b:[1,2]}};"),
            "requires 1 rows and 1 columns",
        ),
        (
            format!("sheaf Bad[D,3] {stalks} {{a<b:[1],a<b:[1]}};"),
            "duplicate restriction",
        ),
        (
            format!("sheaf Bad[D,3] {stalks} {{a<b:[1]}};"),
            "missing restriction",
        ),
        (format!("sheaf Bad[D,9] {stalks} {{}};"), "is not prime"),
        ("sheaf Bad[D,3] {a:1,b:1,d:1} {};".into(), "missing stalk"),
        ("sheaf Bad[D,3] {a:1,a:0} {};".into(), "duplicate stalk"),
        ("sheaf Bad[E,3] {} {};".into(), "unknown named sheaf base"),
        ("sheaf Bad[D,3] {e:1} {};".into(), "unknown point 'e'"),
    ] {
        let error = diagnostic(interpret(&format!("{diamond}{statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }

    let good = format!("{diamond}sheaf F[D,3] {stalks} {{a<b:[1],a<c:[1],b<d:[-1],c<d:[2]}};");
    for (statement, reason) in [
        ("show shTerm[F,[],0];", "nonempty strict chain"),
        ("show shTerm[F,[b,c],0];", "strict increasing chain"),
        ("show shTerm[F,[a,c],1];", "terminal stalk coordinate"),
        (
            "show shClass(shTerm[F,[a],0]);",
            "requires a closed sheaf cochain",
        ),
        (
            "show shgen[F,1,0];",
            "cohomology generator 0 is out of range",
        ),
        ("show shRep(shTerm[F,[a],0]);", "requires a ShH class"),
    ] {
        let error = diagnostic(interpret(&format!("{good}{statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
}

#[test]
fn commuting_paths_are_verified_without_listing_the_exponentially_many_chains() {
    let mut vertices = Vec::new();
    let mut covers = Vec::new();
    for stage in 0..60 {
        for branch in 0..2 {
            let vertex = format!("p{stage}_{branch}");
            vertices.push(vertex.clone());
            if stage != 0 {
                for previous_branch in 0..2 {
                    covers.push(format!("p{}_{previous_branch}<{vertex}", stage - 1));
                }
            }
        }
    }
    let stalks: Vec<_> = vertices.iter().map(|point| format!("{point}:1")).collect();
    let maps: Vec<_> = covers.iter().map(|edge| format!("{edge}:[1]")).collect();
    let program = format!(
        "space X {{{}}};\nsheaf F[X,3] {{{}}} {{{}}};\nshdim[F,600];",
        covers.join(","),
        stalks.join(","),
        maps.join(","),
    );
    assert_eq!(stdout(interpret(&program)), "shdim[F, 600] = 0\n");
}
