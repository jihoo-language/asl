mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn explicit_localization_contraction_recovers_classical_ext_classes() {
    assert_eq!(
        stdout(example("rho-localization")),
        concat!(
            "RhoLoc(0, 1, 0){rholocword[1, 1]{}}\n",
            "Ext(1, 1){g0}\n",
            "Ext(1, 2){g0}\n",
        ),
    );
}

#[test]
fn localized_homotopy_checks_the_retraction_and_side_identities() {
    let source = r"
        let x = rholocword[-2,0]{(0,[1]),(3,[]),(1,[1])};
        let a = rholocword[-2,3]{(0,[1]),(1,[1])};
        let b = rholocword[-2,2]{(1,[1]),(1,[1])};
        let c = rholocword[-2,1]{(2,[1]),(1,[1])};
        let d = rholocword[-2,0]{(3,[1]),(1,[1])};
        assert rhoHomotopy(x) == add(add(a,b),add(c,d));
        assert add(cobarD(rhoHomotopy(x)),rhoHomotopy(cobarD(x))) == add(x,rhoSection[-2](rhoProject(x)));
        assert rhoProject(rhoHomotopy(x)) == zero[Cobar(2,6)];
        assert rhoHomotopy(rhoHomotopy(x)) == zero[RhoLoc(1,-2,-8)];
        assert cobarD(cobarD(x)) == zero[RhoLoc(5,-2,-8)];

        let t = rholocword[0,3]{};
        assert rhoHomotopy(cobarD(t)) == t;
        let y = rhoSection[-2](cword{[0,1]});
        assert rhoProject(y) == cword{[0,1]};
        assert rhoHomotopy(y) == zero[RhoLoc(0,-2,-5)];
        assert rhoProject(cobarD(y)) == cobarD(rhoProject(y));
        assert cobarD(y) == rhoSection[-2](cobarD(cword{[0,1]}));
        assert rhoProject(rhoSection[-8](cword{[1],[0,1]})) == cword{[1],[0,1]};
    ";
    assert_eq!(stdout(interpret(source)), "");
}

#[test]
fn real_generators_and_balanced_products_localize_as_chain_maps() {
    let source = r"
        let tauone = rmcword[0,0]{([],[0,1])};
        let u2 = rholocword[3,0]{(2,[])};
        let txi = rholocword[3,1]{(0,[1])};
        let uxi = rholocword[3,0]{(1,[1])};
        assert rhoLocalize(tauone) == add(u2,add(txi,uxi));
        let a = rmcword[0,1]{([],[1,1]),([1],[0,0,1])};
        let b = rmcword[2,2]{([],[1])};
        let c = rmcword[0,1]{([],[0,1])};
        let x = rhoLocalize(a);
        let y = rhoLocalize(b);
        let z = rhoLocalize(c);
        assert rhoLocalize(cobarD(a)) == cobarD(x);
        assert rhoLocalize(concat(a,b)) == concat(x,y);
        assert rhoLocalize(rhoShift[2](b)) == rhoShift[2](y);
        assert rhoShift[-2](rhoShift[2](x)) == x;
        assert concat(concat(x,y),z) == concat(x,concat(y,z));
        assert cobarD(concat(x,y)) == add(concat(cobarD(x),y),concat(x,cobarD(y)));
        assert rhoProject(concat(x,y)) == concat(rhoProject(x),rhoProject(y));
        assert concat(rhoSection[-2](cword{[1]}),rhoSection[1](cword{[0,1]})) == rhoSection[-1](cword{[1],[0,1]});

        assert rhoProject(rhoLocalize(rmcword[0,0]{([1],[])})) == cword{[1]};
        assert rhoProject(rhoLocalize(rmcword[0,0]{([2],[])})) == cword{[2]};
        assert rhoProject(rhoLocalize(rmcword[0,0]{([4],[])})) == cword{[4]};
        assert rhoProject(rhoLocalize(rmcword[0,0]{([8],[])})) == cword{[8]};
    ";
    assert_eq!(stdout(interpret(source)), "");
}

#[test]
fn fixed_signed_localized_grades_have_finite_bases() {
    let source = r"
        card(\x : RhoLoc(0,-4,-6) => S.hi);
        card(\x : RhoLoc(1,-3,-4) => S.hi);
        card(\x : RhoLoc(1,0,-2) => S.hi);
        card(\x : RhoLoc(1,-2,-1) => S.hi);
        card(\x : RhoLoc(3,4,2) => S.hi);
        let q = \x : RhoLoc(1,-3,-4) => rhoProject(x);
        assert q(rholocword[-3,0]{(0,[1])}) == cword{[1]};
        assert q(rholocword[-3,0]{(1,[])}) == zero[Cobar(1,1)];
    ";
    assert_eq!(
        stdout(interpret(source)),
        "card = 2 / 2\ncard = 4 / 4\ncard = 32 / 32\ncard = 1 / 1\ncard = 1 / 1\n"
    );
}

#[test]
fn localized_operators_check_types_and_supported_degrees() {
    for (source, reason) in [
        (
            "show rholocword[0,1]{(0,[])};",
            "rho-localized word factors must have positive degree",
        ),
        (
            "show rhoProject(zero[RhoLoc(0,-2,-1)]);",
            "rho-localized degree t-w is outside the supported classical range",
        ),
        (
            "show rhoHomotopy(rholocword[0,1]{});",
            "rhoHomotopy requires positive cobar length",
        ),
        (
            "show rhoSection[0](resA[1](cword{[1]}));",
            "rhoSection requires a classical Cobar cochain",
        ),
        (
            "show rhoShift[-1](rmcword[0,0]{});",
            "rhoShift on RMCobar requires an unsigned exponent",
        ),
        (
            "show concat(rholocword[0,0]{},cword{});",
            "concat requires a RhoLoc cochain",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(reason), "{source}\n{error}");
    }
}
