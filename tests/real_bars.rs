mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn non_dyadic_bars_compute_actual_spectral_differentials_and_nulls() {
    assert_eq!(
        stdout(example("real-motivic-bars")),
        concat!(
            "rmbars[2, 4] = [(2, 1), (4, 2), (8, infinite)]\n",
            "rmbarRank[2, 4, 4, 1] = 2\n",
            "rmbarRank[2, 4, 4, 2] = 1\n",
            "rhoorder = infinite\n",
            "rmbars[2, 5] = [(5, 1), (8, 6), (9, 1), (10, infinite)]\n",
            "rhoorder = 6\n",
            "RhoBock(6, 2, 2, -3, 6){g0}\n",
            "rmbars[3, 6] = [(6, 1), (10, 1), (11, 7), (12, infinite)]\n",
            "rhoorder = 7\n",
            "RhoBock(7, 3, 4, -2, 7){g0}\n",
        ),
    );
}

#[test]
fn shifted_bars_span_their_actual_real_motivic_components() {
    let source = r"
        let a = rmbarGen[1,1,0];
        let b = rmbarGen[1,1,1];
        let shifted = rhoShift[1](b);
        let z = add(a,shifted);
        assert a == extClass(rmcword[0,0]{([],[1])});
        assert rhoProject(shifted) == extgen[1,1,0];
        assert rhoProject(a) == zero[Ext(1,1)];
        assert rhoShift[1](z) == rhoShift[1](shifted);
        assert cobarD(rmbarNull[1,1,0]) == rhoShift[1](rmbarRep[1,1,0]);
        card(\u : RMExt(1,1,0) => S.hi);
        rmbarRank[1,1,1,0];
        rmbarRank[1,1,1,1];
        rmbarRank[1,1,-8,90];
        rmbarRank[1,1,3,8];
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "card = 4 / 4\n",
            "rmbarRank[1, 1, 1, 0] = 2\n",
            "rmbarRank[1, 1, 1, 1] = 1\n",
            "rmbarRank[1, 1, -8, 90] = 1\n",
            "rmbarRank[1, 1, 3, 8] = 0\n",
        ),
    );
}

#[test]
fn the_reducer_handles_high_filtration_and_empty_groups() {
    let source = r"
        rmbars[0,0];
        rmbars[0,1];
        rmbars[1,8];
        rmbars[4,8];
        rmbars[8,0];
        rmbars[1,-2];
        rmbarRank[1,-2,0,0];
        assert rmbarGen[0,0,0] == rmextgen[0,0,0,0];
        assert rmbarRep[0,0,0] == rmcword[0,0]{};
        let q = rmbarGen[4,8,1];
        assert rhoShift[0](q) == q;
        assert extClass(rmbarRep[4,8,1]) == q;
        assert cobarD(rmbarNull[4,8,1]) == rhoShift[2](rmbarRep[4,8,1]);
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "rmbars[0, 0] = [(0, infinite)]\n",
            "rmbars[0, 1] = []\n",
            "rmbars[1, 8] = [(8, 8), (16, infinite)]\n",
            "rmbars[4, 8] = [(4, 1), (13, 2)]\n",
            "rmbars[8, 0] = []\n",
            "rmbars[1, -2] = []\n",
            "rmbarRank[1, -2, 0, 0] = 0\n",
        ),
    );
}

#[test]
fn bar_operations_check_indices_signed_shifts_and_infinite_intervals() {
    for (source, reason) in [
        (
            "show rmbarGen[2,4,3];",
            "rho bar index 3 is out of range for 3 bars",
        ),
        (
            "show rmbarNull[2,4,2];",
            "rmbarNull requires a finite rho bar",
        ),
        (
            "show rmbarRep[0,1,0];",
            "rho bar index 0 is out of range for 0 bars",
        ),
        (
            "rmbarRank[0,0,-170141183460469231731687303715884105728,1];",
            "rmbarRank target degree is out of range",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(reason), "{source}\n{error}");
    }
}
