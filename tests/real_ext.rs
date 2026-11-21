mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn localized_and_complex_maps_distinguish_torsion_and_free_classes() {
    assert_eq!(
        stdout(example("real-motivic-ext")),
        concat!(
            "rmextdim[1, 1, 0] = 2\n",
            "rhoorder = 1\n",
            "rhoorder = infinite\n",
            "rhoorder = 1\n",
            "rhoorder = 32\n",
            "RMCobar(0, 0, -32){rmcword[0, 32]{}}\n",
        ),
    );
}

#[test]
fn true_real_ext_classes_and_maps_preserve_products_and_degrees() {
    let source = r"
        let one = rmextgen[0,0,0,0];
        let rho = extClass(rmcword[1,0]{});
        let h1 = extClass(rmcword[0,0]{([1],[])});
        let h2 = extClass(rmcword[0,0]{([2],[])});
        assert extClass(extRep(h2)) == h2;
        assert yoneda(one,h1) == h1;
        assert yoneda(h2,one) == h2;
        assert yoneda(rho,h1) == rhoShift[1](h1);
        assert yoneda(h1,rho) == rhoShift[1](h1);
        assert rhoShift[5](yoneda(h1,h2)) == yoneda(rhoShift[2](h1),rhoShift[3](h2));
        assert realComplex(yoneda(h1,h2)) == yoneda(realComplex(h1),realComplex(h2));
        assert rhoProject(yoneda(h1,h2)) == yoneda(rhoProject(h1),rhoProject(h2));
        assert rhoProject(h2) == extgen[1,2,0];
        assert rhoProject(rhoShift[6](h2)) == extgen[1,2,0];
        assert realComplex(rhoShift[6](h2)) == zero[MExt(1,-2,-4)];
        assert rhoProject(rho) == extgen[0,0,0];
        let g = \x : RMExt(1,2,1) => rhoProject(x);
        assert g(h1) == extgen[1,1,0];
        card(\x : RMExt(1,1,0) => S.hi);
        rmextdim[1,-5,-6];
        rmextdim[1,-5,-4];
        rmextdim[2,0,1];
        rhoorder(rhoShift[6](h2));
        rhoorder(zero[RMExt(65,130,0)]);
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "card = 4 / 4\n",
            "rmextdim[1, -5, -6] = 1\n",
            "rmextdim[1, -5, -4] = 0\n",
            "rmextdim[2, 0, 1] = 0\n",
            "rhoorder = infinite\n",
            "rhoorder = 0\n",
        ),
    );
}

#[test]
fn powers_of_tau_zero_have_exact_rho_exponents_and_witnesses() {
    for power in [1, 2, 4, 8, 16, 32] {
        let source = format!(
            "let z = rmcword[0,0]{{([],[{power}])}};\n\
             rhoorder(z);\n\
             let w = rhoNull[{power}](z);\n\
             assert w == rmcword[0,{power}]{{}};\n\
             assert cobarD(w) == rhoShift[{power}](z);\n\
             assert extPrimitive(cobarD(w)) == w;\n\
             rhoorder(cobarD(w));"
        );
        assert_eq!(
            stdout(interpret(&source)),
            format!("rhoorder = {power}\nrhoorder = 0\n"),
        );
    }
}

#[test]
fn a_long_raw_tower_has_an_exact_sparse_null() {
    let factors = std::iter::repeat_n("([],[1])", 30)
        .collect::<Vec<_>>()
        .join(",");
    let source = format!(
        "let z = rmcword[0,0]{{{factors}}};
         rhoorder(z);
         assert cobarD(rhoNull[1](z)) == rhoShift[1](z);"
    );
    assert_eq!(
        stdout(interpret(&source)),
        "rhoorder = 1
"
    );
}

#[test]
fn exact_boundary_witnesses_are_typed_and_zero_is_killed_at_power_zero() {
    let source = r"
        let z = zero[RMExt(2,3,1)];
        assert rhoNull[0](z) == zero[RMCobar(1,3,1)];
        let t = rmcword[0,1]{};
        let d = cobarD(t);
        assert extClass(d) == zero[RMExt(1,0,-1)];
        assert rhoNull[0](d) == t;
        let h0 = extClass(rmcword[0,0]{([],[1])});
        let h00 = yoneda(h0,h0);
        assert extClass(extRep(h00)) == h00;
        assert cobarD(rhoNull[1](h00)) == rhoShift[1](extRep(h00));
        assert rhoShift[1](h00) == zero[RMExt(2,1,-1)];
        rhoorder(rhoShift[1](h00));
    ";
    assert_eq!(stdout(interpret(source)), "rhoorder = 0\n");
}

#[test]
fn witnesses_reject_noncocycles_nontorsion_and_nonminimal_powers() {
    for (source, reason) in [
        (
            "rhoorder(rmcword[0,1]{});",
            "requires a closed real-motivic cobar cochain",
        ),
        (
            "show extClass(rmcword[0,1]{});",
            "requires a closed real-motivic cobar cochain",
        ),
        (
            "show rhoNull[0](rmcword[0,0]{([],[1])});",
            "rhoNull power does not kill",
        ),
        (
            "show rhoNull[2](rmcword[0,0]{([],[1])});",
            "rhoNull power is not the minimal exponent",
        ),
        (
            "show rhoNull[6](rmextgen[1,2,1,0]);",
            "rhoNull power does not kill",
        ),
        (
            "show rhoNull[1](zero[RMExt(1,2,1)]);",
            "rhoNull power is not the minimal exponent",
        ),
        (
            "show rhoNull[0](zero[RMExt(0,0,0)]);",
            "rhoNull requires positive cobar length",
        ),
        (
            "show extPrimitive(rmcword[0,0]{([1],[])});",
            "requires a real-motivic cobar boundary",
        ),
        (
            "show rhoProject(zero[RMExt(0,-2,-1)]);",
            "rho-localized degree t-w is outside the supported classical range",
        ),
        (
            "show rhoShift[-1](rmextgen[1,2,1,0]);",
            "rhoShift on RMExt requires an unsigned exponent",
        ),
        (
            "show yoneda(rmextgen[1,2,1,0],mextgen[1,2,1,0]);",
            "yoneda requires an RMExt class",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(reason), "{source}\n{error}");
    }
}
