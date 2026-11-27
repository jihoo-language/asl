mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn signed_cech_operations_and_genus_ten_final_types() {
    assert_eq!(stdout(example("rank-as-odd-derham")), concat!(
        "ASOddDeRham(Curve){alpha(0,1) + 2*eta(0,1)}\n",
        "ASOddDeRham(Curve){alpha(0,0) + alpha(0,1)}\n",
        "ASOddDeRham(Curve){eta(0,0) + 2*eta(0,1)}\n",
        "PrimeScalar(3){2}\n",
        "PrimeScalar(3){1}\n",
        "rvASOddGram[Curve] = [[0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1], [2, 0, 0, 0, 0, 2], [0, 2, 0, 0, 0, 0], [0, 0, 2, 1, 0, 0]]\n",
        "rvASOddEO[Curve] = [0, 0, 1]\n",
        "rvASOddEO[Pure] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]\n",
        "rvASOddaNumber[Cubic] = 6\n",
        "rvASOddaNumber[Quartic] = 6\n",
        "rvASOddEO[Cubic] = [0, 0, 0, 1, 2, 2, 2, 2, 3, 4]\n",
        "rvASOddEO[Quartic] = [0, 1, 1, 2, 2, 2, 2, 2, 3, 4]\n",
    ));
}

#[test]
fn frobenius_verschiebung_and_deck_action_preserve_the_polarized_relations() {
    let mut source = String::from("rvasodd C[3] {(4,1),(2,1)};\n");
    for i in 0..6 {
        writeln!(source, "let x{i}=rvASOddDeRhamGen[C,{i}];").unwrap();
        writeln!(source, "let t{i}=rvASOddDeRhamTranslate[1](x{i});").unwrap();
        writeln!(
            source,
            "assert rvASOddF(rvASOddV(x{i})) == zero[ASOddDeRham(C)];"
        )
        .unwrap();
        writeln!(
            source,
            "assert rvASOddV(rvASOddF(x{i})) == zero[ASOddDeRham(C)];"
        )
        .unwrap();
        writeln!(
            source,
            "assert rvASOddF(t{i}) == rvASOddDeRhamTranslate[1](rvASOddF(x{i}));"
        )
        .unwrap();
        writeln!(
            source,
            "assert rvASOddV(t{i}) == rvASOddDeRhamTranslate[1](rvASOddV(x{i}));"
        )
        .unwrap();
        writeln!(
            source,
            "assert rvASOddDeRhamTranslate[1](rvASOddDeRhamTranslate[1](t{i})) == x{i};"
        )
        .unwrap();
    }
    for i in 0..6 {
        writeln!(
            source,
            "assert rvASOddDeRhamPair(x{i},x{i}) == zero[PrimeScalar(3)];"
        )
        .unwrap();
        for j in 0..6 {
            writeln!(source, "assert rvASOddDeRhamPair(rvASOddF(x{i}),x{j}) == rvASOddDeRhamPair(x{i},rvASOddV(x{j}));").unwrap();
            writeln!(
                source,
                "assert rvASOddDeRhamPair(t{i},t{j}) == rvASOddDeRhamPair(x{i},x{j});"
            )
            .unwrap();
            writeln!(source, "assert add(rvASOddDeRhamPair(x{i},x{j}),rvASOddDeRhamPair(x{j},x{i})) == zero[PrimeScalar(3)];").unwrap();
        }
    }
    source.push_str(
        r"
        assert rvASOddF(pscale[2](x3)) == pscale[2](rvASOddF(x3));
        assert rvASOddV(add(x3,x5)) == add(rvASOddV(x3),rvASOddV(x5));
        assert rvASOddDeRhamPair(x3,x5) == rvPrimeElt[3,-1];
        assert rvASOddDeRhamTranslate[-1](t3) == x3;
        assert rvASOddDeRhamIncl(rvASOddOmega[C,1,0]) == x1;
    ",
    );
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn zero_and_sparse_large_genus_and_scalar_types() {
    let input = r"
        rvasodd Line[3] {(1,1)};
        rvASOddEO[Line];
        rvASOddGram[Line];
        assert rvASOddV(zero[ASOddDeRham(Line)]) == zero[ASOddDeRham(Line)];
        card (\x : ASOddDeRham(Line) => S.hi);
        card (\x : PrimeScalar(5) => S.hi);
    ";
    assert_eq!(
        stdout(interpret(input)),
        "rvASOddEO[Line] = []\nrvASOddGram[Line] = []\ncard = 1 / 1\ncard = 5 / 5\n"
    );
    let degree = (1usize << (usize::BITS / 2)) + 1;
    let setup = format!("rvasodd Huge[3] {{({degree},1)}};");
    assert_eq!(
        stdout(interpret(&format!(
            "{setup} show rvASOddV(rvASOddDeRhamGen[Huge,2]);"
        ))),
        "ASOddDeRham(Huge){alpha(0,0)}\n",
    );
    for query in ["rvASOddEO[Huge];", "rvASOddGram[Huge];"] {
        let error = diagnostic(interpret(&format!("{setup} {query}")));
        assert!(error.contains("storage limit"), "{error}");
    }
    let degree = usize::MAX - 1;
    let error = diagnostic(interpret(&format!(
        "rvasodd Oversized[3] {{({degree},1)}}; rvASOddEO[Oversized];"
    )));
    assert!(
        error.contains("de Rham dimension is out of range"),
        "{error}"
    );
}

#[test]
fn wrong_cover_basis_and_operation_are_rejected() {
    let setup = "rvasodd A[3] {(4,1)}; rvasodd B[3] {(4,1)};";
    for (statement, reason) in [
        (
            "show rvASOddDeRhamGen[A,6];",
            "out of range for dimension 6",
        ),
        (
            "show rvASOddF(rvASOddOmega[A,0,0]);",
            "an ASOddDeRham class",
        ),
        (
            "show rvASOddDeRhamIncl(rvASOddDeRhamGen[A,0]);",
            "an ASOddDiff differential",
        ),
        (
            "show rvASOddDeRhamPair(zero[ASOddDeRham(A)],zero[ASOddDeRham(B)]);",
            "expected type ASOddDeRham(A)",
        ),
        ("show zero[PrimeScalar(9)];", "9 is not prime"),
        (
            "show rvASOddDeRhamGen[Missing,0];",
            "unknown odd Artin–Schreier cover",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup} {statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
}
