mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn divided_connecting_cochain_detects_the_nonzero_page_differential() {
    assert_eq!(
        stdout(example("bp-vadic")),
        concat!(
            "BPVPage(2, 2, 2, 0, 1, 12, 0){g0}\n",
            "BPVPage(2, 2, 2, 1, 1, 12, 2){g1}\n",
            "BPI(2, 2, 1, 1){v1^2*t1_1^4}\n",
            "bpVdim[2, 2, 3, 0, 1, 12, 0] = 0\n",
            "BPVPage(3, 2, 9, 1, 1, 144, 9){g0}\n",
            "BPI(3, 2, 1, 1){v1^9*t1_1^27}\n",
            "bpVdim[3, 2, 10, 0, 1, 144, 0] = 0\n",
        )
    );
}

#[test]
fn prime_power_differentials_come_from_normalized_quotient_cocycles() {
    for (p, power) in [(2, 4), (2, 8), (2, 16), (3, 3), (3, 9), (5, 5), (7, 7)] {
        let t = 2 * (p * p - 1) * power;
        let lower_page = power - 1;
        let next_page = power + 1;
        let t1_power = p * power;
        let output_grade = format!("BPVPage({p},2,{power},1,1,{t},{power})");
        let program = format!(
            r"
            let v = bpReduce[1](bpv[{p},2,0,2]);
            let w = bpPow[{power}](v);
            let early = bpVClass[{lower_page},{t},0](w);
            let x = bpVNext(early);
            assert x == bpVClass[{power},{t},0](w);
            assert bpVPrev(x) == early;
            let y = bpVD(x);
            let v1 = bpReduce[1](bpv[{p},2,1,1]);
            let t1 = bpReduce[1](bpt[{p},2,1,1,1]);
            let leading = bpMul(bpPow[{power}](v1),bpPow[{t1_power}](t1));
            assert y == bpVClass[{power},{t},{power}](leading);
            assert bpVD(y) == zero[BPVPage({p},2,{power},2,1,{t},{twice_power})];
            show y;
            bpVdim[{p},2,{next_page},0,1,{t},0];
        ",
            twice_power = 2 * power,
        );
        let result = stdout(interpret(&program));
        assert!(
            result.starts_with(&format!("{}{{g0}}\n", output_grade.replace(',', ", "))),
            "{result}"
        );
        assert!(result.ends_with(" = 0\n"), "{result}");
    }
}

#[test]
fn corrected_lifts_and_a_higher_invariant_generator() {
    let program = r"
        let lower = bpReduce[1](bpMul(bpv[2,2,1,2],bpt[2,2,1,1,1]));
        let v1 = bpReduce[1](bpv[2,2,1,1]);
        let t1 = bpReduce[1](bpt[2,2,1,1,1]);
        let t2 = bpReduce[1](bpt[2,2,1,1,2]);
        let corrected = add(lower,bpMul(v1,add(bpPow[3](t1),t2)));
        let x = bpVNext(bpVClass[1,8,0](lower));
        assert bpVRep(x) == corrected;
        assert bpVClass[2,8,0](corrected) == x;
        assert bpVNext(bpVPrev(x)) == x;
        let top = bpPow[2](bpReduce[2](bpv[2,3,0,3]));
        let y = bpVD(bpVClass[2,28,0](top));
        let leading = bpMul(bpPow[2](bpReduce[2](bpv[2,3,1,2])),bpPow[8](bpReduce[2](bpt[2,3,1,1,1])));
        assert y == bpVClass[2,28,2](leading);
        show y;
        bpVdim[2,3,3,0,2,28,0];
        bpVdim[2,2,1,1,1,12,2];
        bpIExtdim[2,2,1,2,8];
        bpVdim[2,2,10,1,1,12,0];
        bpVdim[2,2,10,1,1,12,1];
        bpVdim[2,2,10,1,1,12,2];
        bpVdim[2,2,10,1,1,12,3];
        bpVdim[2,2,10,1,1,12,4];
        bpVdim[2,2,10,1,1,12,5];
        bpVdim[2,2,10,1,1,12,6];
        bpIExtdim[2,2,1,1,12];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "BPVPage(2, 3, 2, 1, 2, 28, 2){g0}\n",
            "bpVdim[2, 3, 3, 0, 2, 28, 0] = 0\n",
            "bpVdim[2, 2, 1, 1, 1, 12, 2] = 2\n",
            "bpIExtdim[2, 2, 1, 2, 8] = 2\n",
            "bpVdim[2, 2, 10, 1, 1, 12, 0] = 0\n",
            "bpVdim[2, 2, 10, 1, 1, 12, 1] = 0\n",
            "bpVdim[2, 2, 10, 1, 1, 12, 2] = 1\n",
            "bpVdim[2, 2, 10, 1, 1, 12, 3] = 0\n",
            "bpVdim[2, 2, 10, 1, 1, 12, 4] = 0\n",
            "bpVdim[2, 2, 10, 1, 1, 12, 5] = 1\n",
            "bpVdim[2, 2, 10, 1, 1, 12, 6] = 0\n",
            "bpIExtdim[2, 2, 1, 1, 12] = 2\n",
        )
    );
}

#[test]
fn odd_signs_scalars_and_filtered_products() {
    let program = r"
        let v2 = bpReduce[1](bpv[3,2,0,2]);
        let raw = bpPow[3](v2);
        let v1 = bpReduce[1](bpv[3,2,1,1]);
        let t = bpReduce[1](bpt[3,2,1,1,1]);
        let a = bpMul(bpPow[3](v1),bpPow[9](t));
        let b = bpMul(bpPow[9](v1),bpPow[3](t));
        assert bpD(raw) == add(a,bpNeg(b));
        let x = bpVClass[3,48,0](raw);
        let y = bpVD(x);
        assert bpVD(pscale[2](x)) == pscale[2](y);
        assert add(y,pscale[2](y)) == zero[BPVPage(3,2,3,1,1,48,3)];
        let shift = bpVClass[3,4,1](bpReduce[1](bpv[3,2,0,1]));
        let shifted = bpCup(shift,x);
        assert bpVD(shifted) == bpCup(shift,y);
        let check = bpVClass[3,52,4](bpMul(v1,a));
        assert bpVD(shifted) == check;
        let z = bpVClass[0,4,0](t);
        let z2 = bpVClass[0,8,0](bpPow[2](t));
        let positive = bpCup(z,bpVD(z2));
        let d = bpVD(bpCup(z,z2));
        assert d == pscale[2](positive);
        assert bpVD(d) == zero[BPVPage(3,2,0,4,1,12,0)];
        show d;
        let u = bpVClass[0,4,1](bpReduce[1](bpv[3,2,0,1]));
        assert bpVNext(bpCup(u,z)) == bpCup(bpVNext(u),bpVNext(z));
        assert bpCup(bpCup(u,u),z) == bpCup(u,bpCup(u,z));
        betti BPVPage(3,2,3,0,1,48,0);
    ";
    let result = stdout(interpret(program));
    assert!(
        result.contains("BPVPage(3, 2, 0, 3, 1, 12, 0){2*g0}"),
        "{result}"
    );
    assert!(
        result.ends_with("betti(BPVPage(3, 2, 3, 0, 1, 48, 0)) = [3]\n"),
        "{result}"
    );
}

#[test]
fn validation_rejects_wrong_invariant_weight_filtration_or_lift() {
    for (term, reason) in [
        (
            "zero[BPVPage(3,2,1,0,3,0,0)]",
            "requires a BPI quotient with an invariant",
        ),
        (
            "bpVClass[1,12,0](bpv[3,2,0,2])",
            "a BPI cochain with an invariant",
        ),
        (
            "bpVClass[1,16,0](zero[BPThick(3,2,0,2,2)])",
            "a BPI cochain with an invariant",
        ),
        (
            "bpVClass[1,16,0](bpReduce[1](bpv[3,2,0,1]))",
            "homogeneous internal degree",
        ),
        (
            "bpVClass[1,16,1](bpReduce[1](bpv[3,2,0,2]))",
            "below the requested",
        ),
        (
            "bpVClass[2,16,0](bpReduce[1](bpv[3,2,0,2]))",
            "page cycle condition",
        ),
        (
            "bpVNext(bpVgen[3,2,1,0,1,16,0,0])",
            "requires a cycle for d_1",
        ),
        (
            "bpVPrev(zero[BPVPage(3,2,0,0,1,16,0)])",
            "zero has no previous page",
        ),
    ] {
        let error = diagnostic(interpret(&format!("show {term};")));
        assert!(error.contains(reason), "{term}\n{error}");
    }
}
