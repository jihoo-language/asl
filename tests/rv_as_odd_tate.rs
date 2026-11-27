mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn signed_tate_coordinates_and_pairing() {
    assert_eq!(
        stdout(example("rank-as-odd-tate")),
        concat!(
            "rvASOddTateJordan[Curve] = [2, 2, 2, 2, 2, 2]\n",
            "rvASOddTatedim[Curve, 0] = 6\n",
            "rvASOddTatedim[Curve, -1] = 6\n",
            "ASOddTate(Curve, 0){2*g3}\n",
            "ASOddTate(Curve, 0){2*g1}\n",
            "PrimeScalar(3){1}\n",
            "PrimeScalar(3){2}\n",
        )
    );
}

#[test]
fn differential_from_the_deck_group_and_induced_maps() {
    let mut source = String::from("rvasodd C[3] {(7,1)};\n");
    for i in 0..12 {
        writeln!(source, "let z{i}=rvASOddDeRhamGen[C,{i}];").unwrap();
        writeln!(
            source,
            "assert rvASOddTateD[0](z{i}) == add(rvASOddDeRhamTranslate[1](z{i}),pscale[2](z{i}));"
        )
        .unwrap();
        writeln!(source, "assert rvASOddTateD[-1](z{i}) == add(z{i},add(rvASOddDeRhamTranslate[1](z{i}),rvASOddDeRhamTranslate[2](z{i})));").unwrap();
        for degree in [0, 1] {
            writeln!(
                source,
                "assert rvASOddTateD[{}](rvASOddTateD[{degree}](z{i})) == zero[ASOddDeRham(C)];",
                degree + 1
            )
            .unwrap();
            writeln!(source, "assert rvASOddTateClass[{}](rvASOddTateD[{degree}](z{i})) == zero[ASOddTate(C,{})];", degree + 1, degree + 1).unwrap();
        }
    }
    for (symbol, degree) in [('x', 0), ('y', -1)] {
        for i in 0..6 {
            writeln!(source, "let {symbol}{i}=rvASOddTategen[C,{degree},{i}];").unwrap();
            writeln!(
                source,
                "assert rvASOddTateClass[{degree}](rvASOddTateRep({symbol}{i})) == {symbol}{i};"
            )
            .unwrap();
            writeln!(
                source,
                "assert rvASOddTateRep(rvASOddTategen[C,{},{}]) == rvASOddTateRep({symbol}{i});",
                degree + 6,
                i
            )
            .unwrap();
            for (operation, raw) in [("F", "F"), ("V", "V")] {
                writeln!(source, "assert rvASOddTate{operation}({symbol}{i}) == rvASOddTateClass[{degree}](rvASOdd{raw}(rvASOddTateRep({symbol}{i})));").unwrap();
            }
        }
    }
    for i in 0..6 {
        for j in 0..6 {
            writeln!(source, "assert rvASOddTatePair(rvASOddTateF(x{i}),y{j}) == rvASOddTatePair(x{i},rvASOddTateV(y{j}));").unwrap();
            writeln!(source, "assert rvASOddTatePair(rvASOddTateF(y{i}),x{j}) == rvASOddTatePair(y{i},rvASOddTateV(x{j}));").unwrap();
            writeln!(source, "assert add(rvASOddTatePair(x{i},y{j}),rvASOddTatePair(y{j},x{i})) == zero[PrimeScalar(3)];").unwrap();
        }
    }
    source.push_str(
        r"
        assert rvASOddTateD[0](z4) == z0;
        assert rvASOddTateClass[-1](z0) == zero[ASOddTate(C,-1)];
        assert rvASOddTateF(y2) == y4;
        assert rvASOddTateF(y3) == y0;
        assert rvASOddTateV(y5) == y0;
        assert rvASOddTatePair(x0,y2) == rvPrimeElt[3,1];
        assert rvASOddTateF(pscale[2](x4)) == pscale[2](rvASOddTateF(x4));
        assert rvASOddTatePair(add(x0,x1),y2) == add(rvASOddTatePair(x0,y2),rvASOddTatePair(x1,y2));
        assert rvASOddTateD[-9223372036854775808](z4) == rvASOddTateD[0](z4);
        assert rvASOddTateD[9223372036854775807](z4) == rvASOddTateD[-1](z4);
        ",
    );
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn other_prime_empty_module_and_finite_enumeration() {
    let mut source = String::from(
        r"
        rvasodd C[5] {(4,1)};
        rvASOddTateJordan[C];
        rvASOddTatedim[C,0];
        rvASOddTatedim[C,1];
        rvasodd Half[5] {(2,1)};
        rvASOddTateJordan[Half];
        card (\x : ASOddTate(Half,-1) => S.hi);
        rvasodd Line[3] {(1,1)};
        rvASOddTateJordan[Line];
        rvASOddTatedim[Line,0];
        card (\x : ASOddTate(Line,8) => S.hi);
        ",
    );
    for degree in 0..2 {
        for i in 0..3 {
            for operation in ["F", "V"] {
                writeln!(source, "assert rvASOddTate{operation}(rvASOddTategen[C,{degree},{i}]) == zero[ASOddTate(C,{degree})];").unwrap();
            }
        }
    }
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "rvASOddTateJordan[C] = [4, 4, 4]\n",
            "rvASOddTatedim[C, 0] = 3\n",
            "rvASOddTatedim[C, 1] = 3\n",
            "rvASOddTateJordan[Half] = [4]\n",
            "card = 5 / 5\n",
            "rvASOddTateJordan[Line] = []\n",
            "rvASOddTatedim[Line, 0] = 0\n",
            "card = 1 / 1\n"
        )
    );
}

#[test]
fn cycles_degrees_covers_and_basis_indices_are_checked() {
    let setup = "rvasodd A[3] {(7,1)}; rvasodd B[3] {(7,1)};";
    for (statement, reason) in [
        (
            "show rvASOddTateClass[0](rvASOddDeRhamGen[A,4]);",
            "input is not a Tate cycle in degree 0",
        ),
        (
            "show rvASOddTatePair(zero[ASOddTate(A,0)],zero[ASOddTate(A,2)]);",
            "opposite Tate parities",
        ),
        (
            "show rvASOddTatePair(zero[ASOddTate(A,0)],zero[ASOddTate(B,1)]);",
            "same cover",
        ),
        (
            "show add(zero[ASOddTate(A,0)],zero[ASOddTate(A,2)]);",
            "expected type ASOddTate(A, 0)",
        ),
        ("show rvASOddTategen[A,0,6];", "out of range"),
        (
            "show rvASOddTateF(rvASOddDeRhamGen[A,0]);",
            "an ASOddTate class",
        ),
        (
            "show rvASOddTateClass[9223372036854775808](rvASOddDeRhamGen[A,0]);",
            "Tate degree is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup} {statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
}
