mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn level_two_framing_scalar_and_verlinde_entry() {
    assert_eq!(
        stdout(example("rv-wzw-modular")),
        concat!(
            "RVCyclo(60){0, 0, 0, 0, 375, 0, 0, 0, 0, 0, 0, 0, 0, 0, 375}\n",
            "RVCyclo(60){3, 0, -1, 0, -1, 0, -3, 0, -3, 0, -1, 0, 4, 0, 3}\n",
            "RVCyclo(60){-1, 0, -1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, -1}\n",
        )
    );
}

#[test]
fn modular_words_act_linearly_on_nonselfdual_states() {
    let program = r"
        let u = rvWzwBasis[3,1,[0,0,0]];
        let f = rvWzwBasis[3,1,[0,0,1]];
        let g = rvWzwBasis[3,1,[0,1,1]];
        let scalar = add(rvCyclo[48,3],neg(rvCyclo[48,2]));
        let v = add(f,rvWzwScale(scalar,g));
        assert add(v,neg(v)) == zero[RVWzw(3,1)];
        assert rvWzwCharge(v) == add(g,rvWzwScale(scalar,f));
        assert rvWzwS(v) == add(rvWzwS(f),rvWzwScale(scalar,rvWzwS(g)));
        assert rvWzwT[2](v) == add(rvWzwT[2](f),rvWzwScale(scalar,rvWzwT[2](g)));
        let x = rvWzwS(rvWzwT[1](v));
        let y = rvWzwS(rvWzwT[1](x));
        let z = rvWzwS(rvWzwT[1](y));
        let b = add(rvCycloScale[384](rvCyclo[48,4]),rvCycloScale[-192](rvCyclo[48,12]));
        assert z == rvWzwScale(b,rvWzwCharge(v));
        assert rvWzwSAdjoint(rvWzwS(v)) == rvWzwScale(rvCycloScale[48](rvCyclo[48,0]),v);
        assert rvWzwCoord[[0,0,0]](v) == zero[RVCyclo(48)];
        assert rvWzwS(rvWzwS(u)) == rvWzwScale(rvCycloScale[-48](rvCyclo[48,0]),u);
    ";
    assert_eq!(stdout(interpret(program)), "");
}

#[test]
fn roots_embeddings_signed_twists_and_level_zero() {
    let program = r"
        let f = rvWzwBasis[3,2,[-7,-7,-6]];
        assert f == rvWzwBasis[3,2,[0,0,1]];
        let twistmin = rvWzwT[-170141183460469231731687303715884105728](f);
        assert twistmin == rvWzwT[52](f);
        assert rvWzwTEntry[3,2,[0,0,1],-1] == rvCycloConj(rvWzwTEntry[3,2,[0,0,1]]);
        assert rvWzwT[-1](f) == rvWzwScale(rvWzwTEntry[3,2,[0,0,1],-1],f);
        assert rvCycloEmbed[60](rvCycloEmbed[15](rvCyclo[5,-2])) == rvCycloEmbed[60](rvCyclo[5,-2]);
        assert rvCycloEmbed[60](rvCyclo[15,-3]) == rvCyclo[60,-12];
        let a = rvCyclo[15,2];
        let b = rvCyclo[15,4];
        assert rvCycloEmbed[60](rvCycloMul(a,b)) == rvCycloMul(rvCycloEmbed[60](a),rvCycloEmbed[60](b));
        let v = rvWzwBasis[4,0,[15,15,15,15]];
        assert rvWzwT[3](v) == v;
        assert rvWzwCharge(v) == v;
        assert rvWzwS(rvWzwS(v)) == rvWzwScale(rvCycloScale[256](rvCyclo[32,0]),v);
    ";
    assert_eq!(stdout(interpret(program)), "");
}

#[test]
fn modular_inputs_and_storage_limits_are_checked() {
    for (program, expected) in [
        (
            "show rvWzwBasis[1,2,[0]];",
            "fusion rank must be between 2 and 24",
        ),
        ("show rvWzwBasis[3,2,[0,1]];", "must match the ring rank"),
        ("show rvWzwTEntry[3,2,[0,3,3]];", "outside its level alcove"),
        (
            "show rvWzwSEntry[3,2,[1,0,0],[0,0,0]];",
            "must be weakly increasing",
        ),
        (
            "show rvWzwBasis[3,1000,[0,0,0]];",
            "cyclotomic order must be between 1 and 4096",
        ),
        (
            "show rvWzwS(rvWzwBasis[3,100,[0,0,0]]);",
            "limit of 4096 integrable weights",
        ),
        ("show rvWzwS(S.hi);", "requires a WZW torus state"),
        (
            "show rvWzwScale(rvCyclo[15,0],rvWzwBasis[3,2,[0,0,0]]);",
            "RVCyclo(60)",
        ),
        (
            "show rvCycloEmbed[15](rvCyclo[6,1]);",
            "source order to divide",
        ),
        (
            r"show \v : RVWzw(3,2) => v;",
            "cannot enumerate a cyclotomic torus state module",
        ),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}\n{error}");
    }
}
