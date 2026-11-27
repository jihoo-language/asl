mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn lens_spaces_with_same_homology_order_have_distinct_values() {
    assert_eq!(
        stdout(example("rv-lens-tv")),
        concat!(
            "RVCycloQ(16){[] / 1}\nRVCycloQ(16){[1] / 1}\n",
            "RVCycloQ(16){[2, 0, -1, 0, 0, 0, 1] / 4}\n",
            "RVCycloQ(16){[1] / 1}\nRVCycloQ(12){[1] / 2}\n",
            "RVCycloQ(12){[] / 1}\nRVCycloQ(12){[1] / 1}\n",
            "RVCycloQ(48){[1] / 3}\nRVCycloQ(48){[1] / 1}\n",
        )
    );
}

#[test]
fn gluing_framing_and_canonical_rationals_are_exact() {
    let program = r"
        let u = rvCycloRat[1](rvCyclo[16,0]);
        assert rvLensTV[2,2,0,1] == u;
        assert rvLensTV[2,2,0,-1] == u;
        assert rvWzwTV[2,2]{} == u;
        assert rvWzwTV[2,2]{S,S} == u;
        assert rvLensTV[2,2,1,0] == rvCycloRat[4](rvCyclo[16,0]);
        let a = rvLensTV[2,2,5,2];
        assert a == rvLensTV[2,2,5,3];
        assert a == rvLensTV[2,2,5,-2];
        assert a == rvLensTV[2,2,5,12];
        assert a == rvWzwTV[2,2]{T[0],S,T[-3],S,T[-2],S};
        let b = rvLensTV[2,2,5,1];
        assert b == rvWzwTV[2,2]{T[-9],S,T[5],S,T[23]};
        let half = rvCycloRat[8](rvCycloScale[4](rvCyclo[16,0]));
        assert add(half,half) == u;
        assert add(a,neg(a)) == zero[RVCycloQ(16)];
        assert rvCycloRat[170141183460469231731687303715884105728](rvCycloScale[-170141183460469231731687303715884105728](rvCyclo[16,0])) == neg(u);
        show rvCycloRat[340282366920938463463374607431768211455](rvCyclo[16,0]);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "RVCycloQ(16){[1] / 340282366920938463463374607431768211455}\n"
    );
}

#[test]
fn long_words_cancel_normalization_before_integer_powers_overflow() {
    let word = std::iter::repeat_n("S,T[0]", 64)
        .collect::<Vec<_>>()
        .join(",");
    let program = format!("assert rvWzwTV[3,2]{{{word}}} == rvCycloRat[1](rvCyclo[60,0]);");
    assert_eq!(stdout(interpret(&program)), "");
    if usize::BITS >= 64 {
        assert_eq!(
            stdout(interpret(
                "show rvLensTV[3,1,12200160415121876738,7540113804746346429];"
            )),
            "RVCycloQ(48){[1] / 3}\n"
        );
    }
}

#[test]
fn lens_and_rational_boundaries_are_diagnosed() {
    let oversized = std::iter::repeat_n("S", 257).collect::<Vec<_>>().join(",");
    let oversized = format!("show rvWzwTV[2,2]{{{oversized}}};");
    for (program, expected) in [
        ("show rvLensTV[2,2,0,0];", "p=0 requires q=1 or q=-1"),
        ("show rvLensTV[2,2,0,3];", "p=0 requires q=1 or q=-1"),
        ("show rvLensTV[2,2,8,4];", "must be coprime"),
        ("show rvLensTV[2,2,3,0];", "must be coprime"),
        (
            "show rvCycloRat[0](rvCyclo[16,0]);",
            "denominator must be positive",
        ),
        (
            "show rvCycloRat[340282366920938463463374607431768211456](rvCyclo[16,0]);",
            "denominator exceeds the u128 range",
        ),
        ("show rvWzwTV[2,2]{U[1]};", "S or T[power]"),
        (oversized.as_str(), "limit of 256 letters"),
        (
            r"show \x: RVCycloQ(16) => x;",
            "cannot enumerate a rational cyclotomic field",
        ),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}
