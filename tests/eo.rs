mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn a_genus_five_zip_specialization_is_not_ordinary_bruhat_order() {
    assert_eq!(
        stdout(example("rank-eo-closure")),
        concat!(
            "EOType(5){0, 0, 1, 2, 2}\n",
            "EOType(5){1, 1, 1, 1, 2}\n",
            "eoDim = 5\n",
            "eoDim = 6\n",
            "S.lo\n",
            "S.hi\n",
            "eoClosureWitness = [2, 1, 3, 4, 5]\n",
            "eoWeyl = [1, 2, 6, 7, 3, 8, 4, 5, 9, 10]\n",
            "eoWeyl = [6, 1, 2, 3, 7, 4, 8, 9, 10, 5]\n",
        )
    );
}

#[test]
fn type_values_from_curves_and_classification_of_small_genus() {
    let input = r"
        rvascover C[[0,1]] {rvmon[[0,1]]{0,5},rvmon[[0,1]]{5,0}};
        let x = rvASEOType[C];
        let sup = eoType{0,0};
        let ord = eoType{1,2};
        assert x == eoType{0,1};
        assert eoBruhat(sup,x) == S.hi;
        assert eoClosure(sup,x) == S.hi;
        assert eoBruhat(x,ord) == S.hi;
        assert eoClosure(x,ord) == S.hi;
        assert eoClosure(ord,x) == S.lo;
        assert eoClosure(x,x) == S.hi;
        eoClosureWitness(ord,x);
        eoClosureWitness(x,x);
        eoDim(eoType{});
        eoWeyl(eoType{});
        eoClosureWitness(eoType{},eoType{});
        card (\v : EOType(4) => S.hi);
        card (\v : EOType(0) => S.hi);
    ";
    assert_eq!(
        stdout(interpret(input)),
        concat!(
            "eoClosureWitness = none\n",
            "eoClosureWitness = [1, 2]\n",
            "eoDim = 0\n",
            "eoWeyl = []\n",
            "eoClosureWitness = []\n",
            "card = 16 / 16\n",
            "card = 1 / 1\n",
        )
    );
}

#[test]
fn impossible_orderings_and_unfinished_search_have_distinct_results() {
    let setup = r"
        let a = eoType{0,0,1,2,2,2,2,2,2};
        let b = eoType{1,1,1,1,2,2,2,2,2};
        let ord = eoType{1,2,3,4,5,6,7,8,9};
    ";
    assert_eq!(
        stdout(interpret(&format!(
            "{setup} show eoBruhat(a,b); show eoClosure(a,ord); show eoClosure(ord,a);"
        ))),
        "S.lo\nS.hi\nS.lo\n"
    );
    for (expression, reason) in [
        ("eoClosure(a,b)", "Levi search supports genus at most 8"),
        (
            "eoClosure(eoType{0},eoType{0,1})",
            "expected type EOType(1)",
        ),
        (
            "eoBruhat(eoType{0,2},eoType{1,2})",
            "increase by zero or one at position 2",
        ),
        ("eoType{1,0}", "increase by zero or one at position 2"),
        ("eoType{2}", "increase by zero or one at position 1"),
        ("eoClosure(S.lo,S.hi)", "an EO final type"),
    ] {
        let error = diagnostic(interpret(&format!("{setup} show {expression};")));
        assert!(error.contains(reason), "{expression}\n{error}");
    }
}
