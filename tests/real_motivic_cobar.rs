mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn unit_transport_and_specialization_are_exact() {
    assert_eq!(
        stdout(example("real-motivic-cobar")),
        concat!(
            "RMCobar(1, 0, -1){rmcword[1, 0]{([], [1])}}\n",
            "RMCobar(1, 0, -2){rmcword[3, 0]{([1], [1])} + rmcword[3, 0]{([], [0, 1])} + rmcword[2, 1]{([1], [])}}\n",
        ),
    );
}

#[test]
fn frobenius_powers_of_tau_have_their_full_normalized_differentials() {
    for power in [1, 2, 4, 8, 16] {
        let source = format!(
            "let x = rmcword[0,{power}]{{}};\n\
             let d = cobarD(x);\n\
             assert d == rhoShift[{power}](rmcword[0,0]{{([], [{power}])}});\n\
             assert cobarD(d) == zero[RMCobar(2,0,-{power})];\n\
             show d;"
        );
        let output = stdout(interpret(&source));
        assert_eq!(
            output.matches("rmcword[").count(),
            2 * power - 1,
            "{output}"
        );
    }
}

#[test]
fn multiplication_is_associative_rho_linear_and_obeys_leibniz() {
    let source = r"
        let x = rmcword[0,1]{([],[1,1])};
        let y = rmcword[0,2]{([1],[1])};
        let z = rmcword[0,1]{([],[0,1])};
        assert concat(concat(x,y),z) == concat(x,concat(y,z));
        assert cobarD(concat(x,y)) == add(concat(cobarD(x),y),concat(x,cobarD(y)));
        assert cobarD(cobarD(concat(x,y))) == zero[RMCobar(4,7,-1)];
        assert rhoShift[2](concat(x,y)) == concat(x,rhoShift[2](y));
        assert rhoShift[2](concat(x,y)) == concat(rhoShift[2](x),y);
        assert cobarD(rhoShift[3](z)) == rhoShift[3](cobarD(z));
        assert realComplex(concat(x,y)) == concat(realComplex(x),realComplex(y));
        assert cobarD(realComplex(x)) == realComplex(cobarD(x));

        let both = rmcword[0,0]{([],[1]),([],[2])};
        let shift = rhoShift[1](rmcword[0,0]{([],[2]),([1],[])});
        let parts = add(add(rmcword[0,1]{([],[1]),([1],[])},rmcword[1,0]{([],[1]),([1],[1])}),
                        rmcword[1,0]{([],[1]),([],[0,1])});
        assert both == add(parts,shift);
    ";
    assert_eq!(stdout(interpret(source)), "");
}

#[test]
fn signed_bases_and_specialization_preserve_degree() {
    let source = r"
        let a = rmcword[2,3]{};
        assert realComplex(a) == zero[MCobar(0,-2,-5)];
        assert realComplex(cobarD(a)) == cobarD(realComplex(a));
        assert concat(a,rmcword[0,0]{}) == a;
        card(\x : RMCobar(0,-2,-5) => S.hi);
        card(\x : RMCobar(0,-2,-1) => S.hi);
        card(\x : RMCobar(1,0,-2) => S.hi);
        card(\x : MCobar(1,-2,-5) => S.hi);
        mextdim[1,-2,-5];
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "card = 2 / 2\n",
            "card = 1 / 1\n",
            "card = 32 / 32\n",
            "card = 1 / 1\n",
            "mextdim[1, -2, -5] = 0\n",
        ),
    );
}

#[test]
fn real_cobar_checks_generators_types_and_signed_bounds() {
    for (source, reason) in [
        (
            "show rmcword[0,0]{([],[])};",
            "real-motivic word factors must have positive degree",
        ),
        (
            "show rhoShift[1](mcword[0]{});",
            "rhoShift requires an RMCobar cochain",
        ),
        (
            "show tauShift[1](rmcword[0,0]{});",
            "tauShift requires an MCobar cochain or MExt class",
        ),
        (
            "show concat(rmcword[0,0]{}, mcword[0]{});",
            "concat requires an RMCobar cochain",
        ),
        (
            "show realize(realComplex(rmcword[1,0]{}));",
            "realize target is outside the supported classical topological degree range",
        ),
        (
            "show rhoShift[1](zero[RMCobar(0,-170141183460469231731687303715884105728,0)]);",
            "real-motivic topological degree is out of range",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(reason), "{source}\n{error}");
    }
}
