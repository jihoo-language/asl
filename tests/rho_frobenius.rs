mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn algebraic_frobenius_respects_chain_maps_and_twisted_products() {
    let source = r"
        let a = rmcword[0,1]{([],[1,1])};
        let b = rmcword[1,0]{([1],[1]),([],[1])};
        let p = concat(a,b);
        assert cobarD(rhoFrob(p)) == rhoFrob(cobarD(p));
        assert rhoFrob(p) == concat(rhoFrob(a),rhoFrob(b));
        assert rhoFrob(rhoShift[3](p)) == rhoShift[6](rhoFrob(p));
        assert rhoFrob(rhoReduce[3](p)) == rhoReduce[6](rhoFrob(p));
        assert rhoFrob(rhoReduce[2](rhoReduce[3](p))) == rhoReduce[4](rhoFrob(rhoReduce[3](p)));
        assert realComplex(rhoFrob(a)) == rhoFrob(realComplex(a));
        let x = rmcword[0,0]{([],[0,1])};
        let y = rmcword[0,0]{([1],[1])};
        assert rhoFrob(add(x,y)) == add(rhoFrob(x),rhoFrob(y));

        let h = extClass(rmcword[0,0]{([],[1])});
        let k = extClass(rmcword[0,0]{([1],[])});
        assert rhoFrob(yoneda(h,h)) == yoneda(rhoFrob(h),rhoFrob(h));
        assert realComplex(rhoFrob(h)) == rhoFrob(realComplex(h));
        assert realize(rhoFrob(realComplex(k))) == extSq[0](realize(realComplex(k)));
        assert rhoFrob(rhoReduce[3](h)) == rhoReduce[6](rhoFrob(h));
        let q = rhoReduce[3](k);
        assert rhoFrob(rhoShift[1](q)) == rhoShift[2](rhoFrob(q));
        assert rhoFrob(yoneda(q,q)) == yoneda(rhoFrob(q),rhoFrob(q));
        assert rhoFrob(extClass(realComplex(cobarD(a)))) == zero[MExt(2,8,0)];
    ";
    assert_eq!(stdout(interpret(source)), "");
}

#[test]
fn frobenius_doubles_page_filtration_and_commutes_with_the_differential() {
    let source = r"
        let t = rhoClass[1,0](rmcword[0,1]{});
        let t2 = rhoFrob(t);
        let t4 = rhoFrob(t2);
        assert t2 == rhobockgen[2,0,0,-2,0,0];
        assert t4 == rhobockgen[4,0,0,-4,0,0];
        assert rhoD(t2) == rhoFrob(rhoD(t));
        assert rhoD(t4) == rhoFrob(rhoD(t2));
        assert rhoD(t4) == rhobockgen[4,1,0,-4,4,0];
        assert rhoFrob(rhoMul(t,t)) == rhoMul(t2,t2);
        let x = rhoClass[0,2](rmcword[2,0]{([],[1,1])});
        assert rhoD(rhoFrob(x)) == rhoFrob(rhoD(x));
    ";
    assert_eq!(stdout(interpret(source)), "");
}

#[test]
fn frobenius_can_kill_a_nonzero_class_in_ext_and_on_a_page() {
    let source = r"
        let h = extClass(rmcword[0,0]{([1],[])});
        let h2 = yoneda(h,h);
        let h4 = yoneda(h2,h2);
        assert h4 == rmextgen[4,8,4,0];
        assert rhoFrob(h4) == zero[RMExt(4,16,8)];
        let a = rhoClass[1,0](extRep(h4));
        assert a == rhobockgen[1,4,8,4,0,0];
        assert rhoFrob(a) == zero[RhoBock(2,4,16,8,0)];
    ";
    assert_eq!(stdout(interpret(source)), "");
}

#[test]
fn the_two_slot_frobenius_has_a_nonzero_raw_connecting_homotopy() {
    let output = stdout(example("rho-frobenius"));
    assert_eq!(
        output,
        "RhoBock(4, 1, 0, -4, 4){g0}\nRMExt(3, 6, 2){0}\nrhoorder = 64\n"
    );
    let source = r"
        let c = rhoReduce[1](rmcword[0,0]{([],[1]),([],[1])});
        let fc = rhoFrob(c);
        show rhoFrobConnectNull(c);
        show rhoConnectRep(fc);
    ";
    let output = stdout(interpret(source));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(lines.len(), 2);
    assert_eq!(lines[0].matches("rmcword").count(), 8);
    assert_eq!(lines[1].matches("rmcword").count(), 6);
}

#[test]
fn connecting_witness_corrects_the_boundary_between_two_canonical_choices() {
    let source = r"
        let c = rhoReduce[1](rmcword[0,0]{([],[1]),([],[1])});
        let x = extClass(c);
        let y = rhoFrob(x);
        let can = rhoReduce[2](add(rmcword[0,2]{([1],[]),([1],[])},rmcword[1,1]{([],[1]),([2],[])}));
        assert extRep(y) == can;
        assert extClass(rhoFrob(c)) == y;
        let change = add(add(rmcword[0,1]{([1],[]),([1],[]),([1],[])},
            rmcword[1,0]{([1],[1]),([1],[]),([1],[])}),
            add(rmcword[1,0]{([],[0,1]),([1],[]),([1],[])},rmcword[0,0]{([],[1]),([],[1]),([2],[])}));
        assert rhoConnectRep(y) == change;
        assert cobarD(rhoFrobConnectNull(x)) == add(rhoFrob(rhoConnectRep(x)),change);
        assert rhoConnect(y) == rhoFrob(rhoConnect(x));
        assert rhoBeta(y) == rhoFrob(rhoBeta(x));
    ";
    assert_eq!(stdout(interpret(source)), "");
}

#[test]
fn frobenius_doubles_the_typed_connecting_order_from_32_to_64() {
    let source = r"
        let a = extClass(rhoReduce[32](rmcword[0,32]{}));
        let fa = rhoFrob(a);
        assert fa == extClass(rhoReduce[64](rmcword[0,64]{}));
        let x = rhoFrob(rhoConnect(a));
        assert x == rhoConnect(fa);
        assert rhoFrob(rhoBeta(a)) == rhoBeta(fa);
        assert cobarD(rhoFrobConnectNull(a)) == add(rhoFrob(rhoConnectRep(a)),rhoConnectRep(fa));
        let witness = rhoNull[64](x);
        assert witness == rmcword[0,64]{};
        assert cobarD(witness) == rhoShift[64](extRep(x));
        rhoorder(x);
    ";
    assert_eq!(stdout(interpret(source)), "rhoorder = 64\n");
}

#[test]
fn frobenius_operations_diagnose_wrong_types_noncycles_and_overflow() {
    let cutoff = usize::MAX;
    for (source, reason) in [
        ("show rhoFrob(S.lo);".to_owned(), "rhoFrob requires"),
        (
            "show rhoFrobConnectNull(rmcword[0,0]{});".to_owned(),
            "an RMRhoCobar cochain or RMRhoExt class",
        ),
        (
            "show rhoFrobConnectNull(rhoReduce[2](rmcword[0,1]{}));".to_owned(),
            "requires a closed rho-truncated cochain",
        ),
        (
            "show rhoConnectRep(rhoReduce[2](rmcword[0,1]{}));".to_owned(),
            "connecting lift is not closed modulo rho^2",
        ),
        (
            format!("show rhoFrob(zero[RhoBock({cutoff},0,0,0,0)]);"),
            "doubled grade is out of range",
        ),
        (
            format!("show rhoFrob(zero[RMRhoExt({cutoff},0,0,0)]);"),
            "doubled grade is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&source));
        assert!(error.contains(reason), "{source}\n{error}");
    }
}
