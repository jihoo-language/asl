mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn fourfold_zero_products_give_the_nonlinear_rank_one_loci() {
    assert_eq!(stdout(example("map-massey-four")), "card = 10 / 16\n");
    let program = r"
        space Three { a, b, c };
        let f = zero[MapExt(Three, One, 1, 0)];
        let g = zero[MapExt(One, One, 0, 0)];
        let h = zero[MapExt(One, One, 1, 0)];
        let k = zero[MapExt(One, Three, 0, 0)];
        let rankOne = emap(\x : Three => Three.a);
        let locus = eBracket(k, h, g, f);
        assert locus(rankOne) == S.hi;
        assert locus(eid[Three]) == S.lo;
        assert eInBracket(k, h, g, f, rankOne) == S.hi;
        assert eInBracket(k, h, g, f, eid[Three]) == S.lo;
        card locus;
    ";
    assert_eq!(stdout(interpret(program)), "card = 50 / 512\n");
}

#[test]
fn a_nonexact_adjacent_product_makes_a_fourfold_product_undefined() {
    let program = "
        let h0 = mapextgen[One, One, 1, 1, 0];
        let h1 = mapextgen[One, One, 1, 2, 0];
        let z = zero[MapExt(One, One, 1, 1)];
        let locus = eBracket(z, z, h0, h0);
        assert locus(zero[MapExt(One, One, 2, 4)]) == S.lo;
        assert eInBracket(z, z, h0, h0, ecomp(h1, h1)) == S.lo;
        card locus;
    ";
    assert_eq!(stdout(interpret(program)), "card = 0 / 2\n");
}

#[test]
fn proper_obstructions_constrain_the_compatible_choices() {
    let program = "
        space Two { a, b };
        let first = zero[MapExt(One, One, 1, 0)];
        let second = eClass(mword[One, 0, 0, Two, 0, 0](cword{}));
        let third = eClass(mword[Two, 0, 1, One, 0, 0](cword{}));
        let fourth = zero[MapExt(One, One, 1, 0)];
        let zeroMap = zero[MapExt(One, One, 0, 0)];
        assert ecomp(third, second) == zeroMap;
        let locus = eBracket(fourth, third, second, first);
        assert locus(zeroMap) == S.hi;
        assert locus(eid[One]) == S.lo;
        assert eInBracket(fourth, third, second, first, eid[One]) == S.lo;
        card locus;
    ";
    assert_eq!(stdout(interpret(program)), "card = 1 / 2\n");

    let undefined = "
        let h0 = mapextgen[One, One, 1, 1, 0];
        let h1 = mapextgen[One, One, 1, 2, 0];
        assert ecomp(h0, h1) == zero[MapExt(One, One, 2, 3)];
        assert ecomp(h1, h0) == zero[MapExt(One, One, 2, 3)];
        let locus = eBracket(h1, h0, h1, h0);
        assert locus(zero[MapExt(One, One, 2, 6)]) == S.lo;
        assert eInBracket(h1, h0, h1, h0, zero[MapExt(One, One, 2, 6)]) == S.lo;
    ";
    assert!(stdout(interpret(undefined)).is_empty());
}

#[test]
fn fivefold_choices_use_both_pair_and_triple_witnesses() {
    let program = r"
        space Two { a, b };
        space Void {};
        let a = zero[MapExt(Two, One, 1, 0)];
        let b = zero[MapExt(One, One, 0, 0)];
        let c = zero[MapExt(One, Void, 1, 0)];
        let d = zero[MapExt(Void, Void, 0, 0)];
        let e = zero[MapExt(Void, Two, 1, 0)];
        let rankOne = emap(\x : Two => Two.a);
        let locus = eBracket(e, d, c, b, a);
        let direct = \x : MapExt(Two, Two, 0, 0) => eInBracket(e, d, c, b, a, x);
        assert locus == direct;
        assert direct(rankOne) == S.hi;
        assert direct(eid[Two]) == S.lo;
        card locus;
    ";
    assert_eq!(stdout(interpret(program)), "card = 10 / 16\n");
}

#[test]
fn direct_membership_avoids_target_enumeration_and_unused_extreme_grades() {
    let program = r"
        space Four { a, b, c, d };
        let f = zero[MapExt(Four, One, 1, 0)];
        let g = zero[MapExt(One, One, 0, 0)];
        let h = zero[MapExt(One, One, 1, 0)];
        let k = zero[MapExt(One, Four, 0, 0)];
        let rankOne = emap(\x : Four => Four.a);
        assert eInBracket(k, h, g, f, rankOne) == S.hi;
        assert eInBracket(k, h, g, f, eid[Four]) == S.lo;
    ";
    assert!(stdout(interpret(program)).is_empty());

    let max = i64::MAX;
    let signed = format!(
        "
        let a = zero[MapExt(One, One, 1, {max})];
        let b = zero[MapExt(One, One, 0, 1)];
        let c = zero[MapExt(One, One, 1, -{max})];
        let d = zero[MapExt(One, One, 0, -1)];
        let e = zero[MapExt(One, One, 1, 0)];
        assert eInBracket(e, d, c, b, a, zero[MapExt(One, One, 0, 0)]) == S.hi;
        "
    );
    assert!(stdout(interpret(&signed)).is_empty());
}

#[test]
fn longer_products_check_composability_output_degree_and_search_limits() {
    for source in [
        "show eBracket(eid[One], eid[One], eid[One], eid[One]);",
        "show eBracket(eid[One], eid[S], eid[One], zero[MapExt(One, One, 2, 0)]);",
        "show eInBracket(eid[One], eid[One], eid[One], zero[MapExt(One, One, 2, 0)], eid[S]);",
    ] {
        diagnostic(interpret(source));
    }
    let inputs = vec!["z"; 400].join(", ");
    let source = format!(
        "let z = zero[MapExt(One, One, 1, 0)]; \
         show eInBracket({inputs}, zero[MapExt(One, One, 2, 0)]);"
    );
    let error = diagnostic(interpret(&source));
    assert!(error.contains("limit"), "{error}");
}
