mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn powers_of_tau_support_genuine_differentials_through_page_sixteen() {
    assert_eq!(
        stdout(example("rho-bockstein")),
        concat!(
            "RhoBock(1, 1, 0, -1, 1){g0}\n",
            "RhoBock(2, 1, 0, -2, 2){g0}\n",
            "RhoBock(4, 1, 0, -4, 4){g0}\n",
            "RhoBock(8, 1, 0, -8, 8){g0}\n",
            "RhoBock(16, 1, 0, -16, 16){g0}\n",
        ),
    );
}

#[test]
fn first_page_agrees_with_shifted_complex_motivic_ext() {
    let source = r"
        let h = rhoClass[1,0](rmcword[0,0]{([1],[])});
        assert h == rhobockgen[1,1,2,1,0,0];
        assert extClass(realComplex(rhoRep(h))) == mextgen[1,2,1,0];
        let shifted = rhoClass[1,3](rmcword[3,0]{([1],[])});
        assert shifted == rhobockgen[1,1,-1,-2,3,0];
        assert rhoClass[1,3](rhoRep(shifted)) == shifted;
        card(\u : RhoBock(1,1,-1,-2,3) => S.hi);
        rhobockdim[1,1,-1,-2,3];
        mextdim[1,2,1];
        rhobockdim[1,1,-1,0,0];
        rhobockdim[0,3,-1,2,9];
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "card = 2 / 2\n",
            "rhobockdim[1, 1, -1, -2, 3] = 1\n",
            "mextdim[1, 2, 1] = 1\n",
            "rhobockdim[1, 1, -1, 0, 0] = 0\n",
            "rhobockdim[0, 3, -1, 2, 9] = 0\n",
        ),
    );
}

#[test]
fn differential_product_and_transitions_use_real_cochain_lifts() {
    let source = r"
        let z = rhoClass[1,0](rmcword[0,0]{});
        let t = rhoClass[1,0](rmcword[0,1]{});
        let h = rhoClass[1,0](rmcword[0,0]{([],[1])});
        assert rhoMul(z,h) == h;
        assert rhoMul(h,z) == h;
        assert rhoMul(z,add(h,h)) == zero[RhoBock(1,1,1,0,0)];
        let x = rhoMul(t,h);
        assert x == rhoClass[1,0](concat(rhoRep(t),rhoRep(h)));
        assert rhoD(x) == add(rhoMul(rhoD(t),h),rhoMul(t,rhoD(h)));
        assert rhoMul(h,rhoMul(t,t)) == rhoMul(rhoMul(h,t),t);
        let square = rhoMul(t,t);
        assert rhoD(square) == zero[RhoBock(1,1,0,-2,1)];
        assert rhoNext(square) == rhoClass[2,0](rmcword[0,2]{});

        let a = rhoClass[0,2](rmcword[2,0]{});
        let b = rhoClass[0,0](rmcword[0,0]{([],[1,1])});
        assert rhoD(rhoMul(a,b)) == rhoMul(a,rhoD(b));
        assert rhoClass[0,2](rhoRep(rhoMul(a,b))) == rhoMul(a,b);
    ";
    assert_eq!(stdout(interpret(source)), "");
}

#[test]
fn enormous_page_indices_stabilize_the_finite_filtration() {
    let max = usize::MAX;
    let source = format!(
        "let x = rhoClass[{max},2](rmcword[2,0]{{}});\n\
         assert x == rhobockgen[{max},0,-2,-2,2,0];\n\
         assert rhoNext(rhoPrev(x)) == x;\n\
         rhobockdim[{max},1,0,-1,1];\n\
         rhobockdim[{max},0,-2,-2,3];"
    );
    assert_eq!(
        stdout(interpret(&source)),
        format!("rhobockdim[{max}, 1, 0, -1, 1] = 0\nrhobockdim[{max}, 0, -2, -2, 3] = 0\n"),
    );
}

#[test]
fn page_operations_reject_wrong_filtrations_noncycles_and_mixed_pages() {
    for (source, reason) in [
        (
            "show rhoClass[0,1](rmcword[0,0]{});",
            "lies below the requested rho-Bockstein filtration",
        ),
        (
            "show rhoClass[2,0](rmcword[0,1]{});",
            "does not satisfy the rho-Bockstein page cycle condition",
        ),
        (
            "show rhoNext(rhoClass[1,0](rmcword[0,1]{}));",
            "rhoNext requires a cycle for d_1",
        ),
        (
            "show rhoPrev(zero[RhoBock(0,0,0,0,0)]);",
            "rhoPrev requires a positive page number",
        ),
        (
            "show rhoMul(zero[RhoBock(0,0,0,0,0)],zero[RhoBock(1,0,0,0,0)]);",
            "rhoMul requires RhoBock classes on the same page",
        ),
        (
            "show rhoClass[1,0](mcword[0]{});",
            "rhoClass requires an RMCobar cochain",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(reason), "{source}\n{error}");
    }
}
