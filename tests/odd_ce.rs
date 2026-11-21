mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn the_three_primary_extension_has_the_signed_third_differential() {
    assert_eq!(
        stdout(example("odd-ce")),
        concat!(
            "PCEAn(3, 1, 3, 4, 14, 4){2*g0}\n",
            "pceandim[3, 1, 4, 3, 14, 1] = 0\n",
            "pceandim[3, 1, 4, 4, 14, 4] = 0\n",
        )
    );
}

#[test]
fn the_five_primary_extension_has_a_negative_fifth_differential_and_a_permanent_fifth_power() {
    let source = "
        let a = resA[1](pword[5]{([], [0])});
        let x = resA[1](pword[5]{([1], [])});
        let v = resA[1](pword[5]{([], [1])});
        let v2 = concat(v, v);
        let v4 = concat(v2, v2);
        let e2 = ceClass[2, 1](concat(v4, x));
        let e3 = ceNext(e2);
        let e4 = ceNext(e3);
        let e5 = ceNext(e4);
        assert e5 == pceangen[5, 1, 5, 5, 44, 1, 0];
        assert ceNext(cePrev(pscale[3](e5))) == pscale[3](e5);
        assert ceClass[5, 1](ceRep(e5)) == e5;
        let beta = add(add(pword[5]{([1], []), ([4], [])}, pword[5]{([4], []), ([1], [])}),
                       pscale[2](add(pword[5]{([2], []), ([3], [])}, pword[5]{([3], []), ([2], [])})));
        let a5 = ceClass[5, 1](a);
        let a2 = ceMul(a5, a5);
        let target = ceMul(ceMul(a2, a2), ceClass[5, 2](resA[1](beta)));
        assert target == pceangen[5, 1, 5, 6, 44, 6, 0];
        assert ceD(e5) == pscale[4](target);
        assert ceD(pscale[2](e5)) == pscale[3](target);
        assert ceNext(target) == zero[PCEAn(5, 1, 6, 6, 44, 6)];
        pceandim[5, 1, 6, 5, 44, 1];
        pceandim[5, 1, 6, 6, 44, 6];
        let vp2 = ceClass[2, 0](concat(v4, v));
        let vp3 = ceNext(vp2);
        let vp4 = ceNext(vp3);
        let vp5 = ceNext(vp4);
        let vp6 = ceNext(vp5);
        let permanent = ceNext(vp6);
        assert extClass(ceRep(permanent)) == pscale[4](pextangen[5, 1, 5, 45, 0]);
    ";
    assert_eq!(
        stdout(interpret(source)),
        "pceandim[5, 1, 6, 5, 44, 1] = 0\npceandim[5, 1, 6, 6, 44, 6] = 0\n"
    );
}

#[test]
fn the_associated_graded_has_zero_differential_and_uses_all_field_digits() {
    let source = "
        let a = resA[0](pword[5]{([], [0])});
        let a0 = ceClass[0, 0](a);
        assert a0 == pceangen[5, 0, 0, 1, 1, 0, 0];
        assert ceD(a0) == zero[PCEAn(5, 0, 0, 2, 1, 0)];
        assert cePrev(ceNext(a0)) == a0;
        assert ceMul(pscale[2](a0), pscale[3](a0)) == ceMul(a0, a0);
        betti PCEAn(5, 0, 1, 2, 2, 0);
        pceandim[5, 100, 0, 1, 4, 0];
    ";
    assert_eq!(
        stdout(interpret(source)),
        "betti(PCEAn(5, 0, 1, 2, 2, 0)) = [5]\npceandim[5, 100, 0, 1, 4, 0] = 0\n"
    );
}

#[test]
fn odd_central_extension_operations_check_profiles_primes_and_differentials() {
    for (expression, expected) in [
        ("ceClass[0, 0](pword[3]{})", "PCobarAn"),
        ("ceNext(pceangen[3, 1, 2, 1, 5, 0, 0])", "cycle for d_2"),
        ("cePrev(pceangen[3, 1, 0, 1, 5, 0, 0])", "greater than zero"),
        ("mayNext(pceangen[3, 1, 0, 1, 5, 0, 0])", "May"),
        (
            "ceMul(zero[PCEAn(3, 1, 0, 0, 0, 0)], zero[PCEAn(3, 0, 0, 0, 0, 0)])",
            "same Steenrod algebra profile",
        ),
        (
            "ceMul(zero[PCEAn(3, 1, 0, 0, 0, 0)], zero[PCEAn(5, 1, 0, 0, 0, 0)])",
            "same prime",
        ),
    ] {
        let error = diagnostic(interpret(&format!("show {expression};")));
        assert!(error.contains(expected), "{expression}\n{error}");
    }
}
