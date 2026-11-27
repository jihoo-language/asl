mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn fork_sheafifies_to_six_sections_and_extends_to_an_arbitrary_named_sheaf() {
    assert_eq!(
        stdout(example("set-sheafification")),
        concat!(
            "setjcheck[F] = no\n",
            "setjcheck[J] = yes\n",
            "setsections[F] = 0\n",
            "setsections[J] = 6\n",
            "setjcheck[G] = yes\n",
            "{root: [0, 2, 1, 3, 5, 4], u: [1, 0], v: [2, 0, 1]}\n",
            "setsieve[O, root, 3] = {u, v}\n",
            "setdense[O, root, 3] = yes\n",
            "setdense[O, root, 1] = no\n",
        ),
    );
}

#[test]
fn dense_sieves_match_double_negation_and_reflection_is_functorial() {
    let output = stdout(interpret(
        r"
        space Fork {x<y,x<z};
        setomega O[Fork];
        setfamily T[Fork] {x:1,y:1,z:1} {x<y:[0],x<z:[0]};
        setnat dense[T,O] {x:[3],y:[1],z:[1]};
        setnat left[T,O] {x:[1],y:[1],z:[0]};
        setnat right[T,O] {x:[2],y:[0],z:[1]};
        assert setNot(setNot(dense)) == setTop[T,O];
        assert setNot(left) == right;
        assert setNot(setNot(left)) == left;
        setdense[O,x,0];
        setdense[O,x,1];
        setdense[O,x,2];
        setdense[O,x,3];
        setdense[O,x,4];

        setfamily F[Fork] {x:2,y:2,z:2} {x<y:[0,1],x<z:[0,1]};
        setj J[F];
        setsections[F];
        setsections[J];
        setnat flip[F,F] {x:[1,0],y:[1,0],z:[1,0]};
        let jf = setJExtend[J](setComp(setJUnit[J],flip));
        assert setJExtend[J](setJUnit[J]) == setId[J];
        assert setComp(jf,setJUnit[J]) == setComp(setJUnit[J],flip);
        assert setComp(jf,jf) == setId[J];
        show setJUnit[J];
        show jf;
    ",
    ));
    assert_eq!(
        output,
        concat!(
            "setdense[O, x, 0] = no\n",
            "setdense[O, x, 1] = no\n",
            "setdense[O, x, 2] = no\n",
            "setdense[O, x, 3] = yes\n",
            "setdense[O, x, 4] = yes\n",
            "setsections[F] = 2\n",
            "setsections[J] = 4\n",
            "{x: [0, 3], y: [0, 1], z: [0, 1]}\n",
            "{x: [3, 2, 1, 0], y: [1, 0], z: [1, 0]}\n",
        ),
    );
}

#[test]
fn empty_sites_and_empty_maximal_fibers() {
    let output = stdout(interpret(
        r"
        space Empty {};
        setfamily F[Empty] {} {};
        setj J[F];
        setjcheck[J];
        assert setJInverse[J] == setJExtend[J](setId[F]);
        setsections[J];

        space Fork {x<y,x<z};
        setfamily A[Fork] {x:0,y:3,z:0} {x<y:[],x<z:[]};
        setj B[A];
        setjcheck[B];
        setsections[B];
        assert setJInverse[B] == setJExtend[B](setId[A]);
    ",
    ));
    assert_eq!(
        output,
        "setjcheck[J] = yes\nsetsections[J] = 1\nsetjcheck[B] = yes\nsetsections[B] = 0\n",
    );
}

#[test]
fn checks_nonsheaf_targets_and_saturates_huge_diagram_products() {
    let program = r"
        setfamily F[S] {lo:2,hi:2} {lo<hi:[0,0]};
        setj J[F];
    ";
    assert_eq!(
        stdout(interpret(&format!("{program}setjcheck[F]; setjcheck[J];"))),
        "setjcheck[F] = no\nsetjcheck[J] = yes\n",
    );
    for (expression, expected) in [
        (
            "show setJInverse[J];",
            "original family to be a double-negation sheaf",
        ),
        (
            "show setJExtend[J](setId[F]);",
            "target must be a double-negation sheaf",
        ),
        (
            "show setJExtend[J](setId[J]);",
            "requires a map from the original family",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{program}{expression}")));
        assert!(error.contains(expected), "{error}");
    }
    let edges = (0..75)
        .map(|i| format!("x<y{i}"))
        .collect::<Vec<_>>()
        .join(",");
    let sizes = (0..75)
        .map(|i| format!("y{i}:2"))
        .collect::<Vec<_>>()
        .join(",");
    let arrows = (0..75)
        .map(|i| format!("x<y{i}:[]"))
        .collect::<Vec<_>>()
        .join(",");
    let program = format!("space X {{{edges}}}; setfamily F[X] {{x:0,{sizes}}} {{{arrows}}};");
    assert_eq!(
        stdout(interpret(&format!("{program}setjcheck[F];"))),
        "setjcheck[F] = no\n",
    );
    let error = diagnostic(interpret(&format!("{program}setj J[F];")));
    assert!(error.contains("limit of 4096 occurrences"), "{error}");
}
