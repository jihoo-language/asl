mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn iterated_deligne_on_the_double_torus_suspension_has_complementary_dimensions() {
    let expected = concat!(
        "shhyperdim[Lower, 0] = 1\nshhyperdim[Upper, 4] = 1\n",
        "shhyperdim[Lower, 1] = 0\nshhyperdim[Upper, 3] = 0\n",
        "shhyperdim[Lower, 2] = 0\nshhyperdim[Upper, 2] = 0\n",
        "shhyperdim[Lower, 3] = 2\nshhyperdim[Upper, 1] = 2\n",
        "shhyperdim[Lower, 4] = 1\nshhyperdim[Upper, 0] = 1\n",
    );
    assert_eq!(stdout(example("sheaf-iterated-deligne")), expected);
    let fixture = include_str!("../examples/sheaf-iterated-deligne.asl");
    for prime in [2, 5] {
        let source = fixture.replace("[Link,3]", &format!("[Link,{prime}]"));
        assert_eq!(stdout(interpret(&source)), expected);
    }
}

#[test]
fn untruncated_second_direct_image_compares_every_actual_source_class() {
    let fixture = include_str!("../examples/sheaf-iterated-deligne.asl");
    let (prefix, _) = fixture.split_once("shdelignecomplex Lower").unwrap();
    for (source, ranks) in [
        ("FirstLower", [1, 0, 2, 1, 0]),
        ("FirstUpper", [1, 2, 0, 1, 0]),
    ] {
        let mut program = format!("{prefix}shleraycomplex Full[k,{source}];\n");
        let mut result = String::new();
        for (degree, rank) in ranks.into_iter().enumerate() {
            program.push_str(&format!(
                "shhyperdim[{source},{degree}]; shhyperdim[Full,{degree}];\n"
            ));
            result.push_str(&format!(
                "shhyperdim[{source}, {degree}] = {rank}\nshhyperdim[Full, {degree}] = {rank}\n"
            ));
            for i in 0..rank {
                program.push_str(&format!(
                    "let u{degree}_{i} = shLerayIncl[Full](shRep(shhypergen[{source},{degree},{i}]));\n\
                     assert shD(u{degree}_{i}) == zero[ShTot(Full,{next})];\n\
                     assert shClass(u{degree}_{i}) == shhypergen[Full,{degree},{i}];\n",
                    next = degree + 1,
                ));
            }
        }
        assert_eq!(stdout(interpret(&program)), result);
    }
}

#[test]
fn a_general_map_to_a_point_preserves_the_total_cochain_differential() {
    let program = r"
        space X {a<b};
        sheaf A[X,5] {a:1,b:1} {a<b:[1]};
        sheaf B[X,5] {a:2,b:2} {a<b:[1,0,0,1]};
        shmap d[A,B] {a:[2,2],b:[2,2]};
        shcomplex K {A,d,B};
        shcomplex Other {A,d,B};
        let to_point = \x : X => One.pt;
        let id = \x : X => x;
        shleraycomplex P[to_point,K];
        let raw = shTotEmbed[K,0](shTerm[A,[a,b],0]);
        let image = shLerayIncl[P](raw);
        assert shD(image) == shLerayIncl[P](shD(raw));
        assert shD(shD(image)) == zero[ShTot(P,3)];
        assert shClass(shLerayIncl[P](shRep(shhypergen[K,1,0]))) == shhypergen[P,1,0];
        shhyperdim[P,0];
        shhyperdim[P,1];
        shhyperdim[P,2];
        shdelignecomplex T[id,K,1];
        shhyperdim[T,0];
        shhyperdim[T,1];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "shhyperdim[P, 0] = 0\nshhyperdim[P, 1] = 1\nshhyperdim[P, 2] = 0\n",
            "shhyperdim[T, 0] = 0\nshhyperdim[T, 1] = 1\n",
        )
    );
    for (statement, reason) in [
        (
            "shleraycomplex P[to_point,K];",
            "sheaf complex 'P' is already defined",
        ),
        ("shleraycomplex Z[to_point,A];", "unknown sheaf complex 'A'"),
        (
            "shdelignecomplex Z[to_point,K,1];",
            "shdelignecomplex requires an injective open embedding",
        ),
        (
            "shleraycomplex Z[id,P];",
            "shleraycomplex requires a complex on the map's named source",
        ),
        (
            "show shLerayIncl[P](zero[ShC(A,0)]);",
            "shLerayIncl requires a ShTot cochain",
        ),
        (
            "show shLerayIncl[P](zero[ShTot(Other,0)]);",
            "requires source complex 'K'",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{program}\n{statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
}
