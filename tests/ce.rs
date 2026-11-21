mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn central_extension_transgressions_leave_the_periodicity_class() {
    assert_eq!(
        stdout(example("cartan-eilenberg")),
        concat!(
            "CEAn(1, 2, 2, 3, 2){g0}\n",
            "CEAn(1, 3, 3, 6, 3){g0}\n",
            "ExtAn(1, 4, 12){g0}\n",
        )
    );
}

#[test]
fn associated_graded_and_first_differential_count_top_generator_occurrences() {
    let program = r"
        let mixed = resA[1](cword{[1, 1]});
        let a = ceClass[0, 0](mixed);
        assert ceD(a) == zero[CEAn(1, 0, 2, 4, 0)];
        let first = ceNext(a);
        let expected = resA[1](add(cword{[1], [0, 1]}, cword{[0, 1], [1]}));
        assert ceD(first) == ceClass[1, 1](expected);
        assert ceD(ceD(first)) == zero[CEAn(1, 1, 3, 4, 2)];
        assert cePrev(first) == a;
        let lift = \x : CEAn(1, 1, 1, 4, 0) => ceClass[1, 0](ceRep(x));
        assert lift == (\x : CEAn(1, 1, 1, 4, 0) => x);
    ";
    assert!(stdout(interpret(program)).is_empty());
}

#[test]
fn a_zero_collapses_and_a_two_detects_the_top_reduced_coproduct() {
    let program = r"
        let u = ceClass[0, 0](resA[0](cword{}));
        let x = ceClass[0, 0](resA[0](cword{[1]}));
        assert ceMul(u, x) == x;
        let h = ceNext(ceNext(x));
        assert ceNext(cePrev(h)) == h;
        let square = ceMul(h, h);
        assert square == ceangen[0, 2, 2, 2, 0, 0];
        assert add(square, square) == zero[CEAn(0, 2, 2, 2, 0)];
        assert extClass(ceRep(square)) == extangen[0, 2, 2, 0];
        ceandim[0, 15, 4, 4, 0];
        ceandim[0, 15, 4, 4, 1];
        betti CEAn(0, 2, 2, 2, 0);

        let v = ceClass[2, 0](resA[2](cword{[0, 0, 1]}));
        let coproduct = resA[2](add(cword{[0, 2], [1]}, cword{[4], [0, 1]}));
        assert cobarD(ceRep(v)) == coproduct;
        assert ceD(v) == ceClass[2, 2](coproduct);
        assert ceD(v) == ceangen[2, 2, 2, 7, 2, 0];
        assert ceNext(ceD(v)) == zero[CEAn(2, 3, 2, 7, 2)];
        let vSquare = ceNext(ceMul(v, v));
        assert ceD(vSquare) == ceangen[2, 3, 3, 14, 3, 0];
        assert ceNext(ceD(vSquare)) == zero[CEAn(2, 4, 3, 14, 3)];
        ceandim[2, 2, 2, 7, 2];
        ceandim[2, 3, 2, 7, 2];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "ceandim[0, 15, 4, 4, 0] = 1\n",
            "ceandim[0, 15, 4, 4, 1] = 0\n",
            "betti(CEAn(0, 2, 2, 2, 0)) = [2]\n",
            "ceandim[2, 2, 2, 7, 2] = 1\n",
            "ceandim[2, 3, 2, 7, 2] = 0\n",
        )
    );
}

#[test]
fn cartan_eilenberg_operations_check_families_profiles_and_cycles() {
    for (term, reason) in [
        ("ceClass[2, 0](cword{[1]})", "requires a CobarAn cochain"),
        (
            "ceClass[2, 0](resA[1](cword{[1, 1]}))",
            "Cartan–Eilenberg page cycle condition",
        ),
        (
            "ceClass[0, 1](resA[1](cword{[0, 1]}))",
            "below the requested Cartan–Eilenberg filtration",
        ),
        (
            "ceNext(ceangen[1, 2, 1, 3, 0, 0])",
            "ceNext requires a cycle for d_2",
        ),
        (
            "ceNext(mayangen[1, 2, 1, 1, 0, 0])",
            "requires a CEAn class",
        ),
        (
            "mayNext(ceangen[1, 2, 1, 3, 0, 0])",
            "requires a May or MapMay class",
        ),
        (
            "cePrev(zero[CEAn(1, 0, 0, 0, 0)])",
            "requires a page greater than zero",
        ),
        (
            "ceMul(zero[CEAn(1, 0, 0, 0, 0)], zero[CEAn(2, 0, 0, 0, 0)])",
            "same Steenrod algebra profile",
        ),
        (
            "ceMul(zero[CEAn(1, 0, 0, 0, 0)], zero[CEAn(1, 1, 0, 0, 0)])",
            "same page",
        ),
        ("resA[0](ceangen[1, 2, 1, 3, 0, 0])", "resA requires"),
    ] {
        let error = diagnostic(interpret(&format!("show {term};")));
        assert!(error.contains(reason), "{term}\n{error}");
    }
    let max = usize::MAX;
    for term in [
        format!("ceNext(zero[CEAn(1, {max}, 0, 0, 0)])"),
        format!("ceD(zero[CEAn(1, {max}, 0, 0, 1)])"),
    ] {
        let error = diagnostic(interpret(&format!("show {term};")));
        assert!(error.contains("out of range"), "{error}");
    }
    let program = format!("ceandim[{max}, 2, {max}, {max}, 0];");
    assert_eq!(
        stdout(interpret(&program)),
        format!("ceandim[{max}, 2, {max}, {max}, 0] = 0\n")
    );
}
