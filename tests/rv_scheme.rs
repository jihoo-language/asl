mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn projective_closure_preserves_the_double_point_and_cusp_singular_scheme() {
    assert_eq!(
        stdout(example("rank-schemes")),
        concat!(
            "rvgroebner[Projective] = [x0, x1^2]\n",
            "rvgroebner[ColonX] = [x2, x1, x0]\n",
            "rvgroebner[ChartZ] = [x0, x1^2]\n",
            "rvhilbert[Thick, 1] = 3\n",
            "rvhilbert[Projective, 1] = 2\n",
            "rvhilbert[Projective, 1000] = 2\n",
            "rvgroebner[Corner] = [x1*x2, x0*x2]\n",
            "rvgroebner[Singular] = [x1^2, x0^2]\n",
            "rvhilbert[Singular, 1000] = 4\n",
        )
    );
}

#[test]
fn nonmonomial_elimination_distinguishes_one_division_from_saturation() {
    let program = r"
        let x = rvmon[[0, 1, 2]]{1,0,0};
        let y = rvmon[[0, 1, 2]]{0,1,0};
        let z = rvmon[[0, 1, 2]]{0,0,1};
        let f = add(x, y);
        let g = add(y, z);
        let h = add(x, z);
        let fg = rvMul(f, g);
        let ffg = rvMul(f, fg);
        rvideal A[[0, 1, 2]] { ffg };
        rvideal B[[0, 1, 2]] { rvMul(f, h) };
        rvcolon Once[A](f);
        rvsaturate All[A](f);
        rvintersect Both[A, B];
        assert rvInIdeal[Once](fg) == S.hi;
        assert rvInIdeal[Once](g) == S.lo;
        assert rvInIdeal[All](g) == S.hi;
        assert rvInIdeal[All](h) == S.lo;
        let lcm = rvMul(ffg, h);
        assert rvInIdeal[Both](lcm) == S.hi;
        assert rvInIdeal[Both](ffg) == S.lo;
        assert rvWitness[Both, 0](lcm) == rvmon[[0, 1, 2]]{0,0,0};
        let y2 = rvMul(y, y);
        rvideal Thick[[0, 1, 2]] { rvMul(f, f), rvMul(f, y), rvMul(f, z), y2 };
        rvprojective Projective[Thick];
        assert rvInIdeal[Projective](f) == S.hi;
        assert rvInIdeal[Projective](x) == S.lo;
        let reduce = \p : RVPoly([0,1,2],1) => rvNormal[Projective](p);
        assert reduce(x) == y;
        rvgroebner[Once];
        rvgroebner[All];
        rvgroebner[Projective];
        rvhilbert[Projective,1000];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvgroebner[Once] = [x1*x2 + x1^2 + x0*x2 + x0*x1]\n",
            "rvgroebner[All] = [x2 + x1]\n",
            "rvgroebner[Projective] = [x1 + x0, x1^2]\n",
            "rvhilbert[Projective, 1000] = 2\n",
        )
    );
}

#[test]
fn empty_and_unit_ideal_rules_and_large_exponents_are_exact() {
    let program = format!(
        "rvideal Empty[[]] {{}};
        rvprojective Irrelevant[Empty];
        assert rvInIdeal[Irrelevant](rvmon[[]]{{}}) == S.hi;
        rvideal Z[[0,1]] {{}};
        let x = rvmon[[0,1]]{{1,0}};
        let one = rvmon[[0,1]]{{0,0}};
        let zeroX = zero[RVPoly([0,1],1)];
        rvcolon ZeroColon[Z](zeroX);
        rvsaturate ZeroSat[Z](zeroX);
        assert rvInIdeal[ZeroColon](one) == S.hi;
        assert rvInIdeal[ZeroSat](one) == S.hi;
        rvprojective ProjectiveZero[Z];
        rvintersect Bottom[Z, ZeroColon];
        assert rvInIdeal[ProjectiveZero](x) == S.lo;
        assert rvInIdeal[Bottom](x) == S.lo;
        rvideal Large[[0,1]] {{ rvmon[[0,1]]{{{},0}}, rvmon[[0,1]]{{0,2}} }};
        rvcolon ColonOne[Large](one);
        rvsaturate SatOne[Large](one);
        assert rvInIdeal[ColonOne](x) == S.lo;
        assert rvInIdeal[SatOne](x) == S.lo;
        rvsaturate SatPower[Large](x);
        assert rvInIdeal[SatPower](one) == S.hi;
        rvprojective Artin[Large];
        assert rvInIdeal[Artin](one) == S.hi;",
        usize::MAX,
    );
    assert_eq!(stdout(interpret(&program)), "");
}

#[test]
fn declarations_check_redefinitions_and_polynomial_rings() {
    let setup = "rvideal J[[0,1]] {}; rvideal H[[0,2]] {};";
    for (statement, reason) in [
        ("rvintersect K[J,H];", "same primitive selection"),
        (
            "rvcolon K[J](rvmon[[0,2]]{1,0});",
            "different primitive selection",
        ),
        ("rvsaturate K[J](S.hi);", "an RVPoly polynomial"),
        ("rvprojective J[H];", "already defined"),
        ("rvprojective K[Missing];", "unknown rank ideal"),
    ] {
        let error = diagnostic(interpret(&format!("{setup}\n{statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
}

#[test]
fn a_derived_ideal_supports_quotient_koszul_tor() {
    let program = r"
        let one = rvmon[[0,1,2]]{0,0,0};
        let y = rvmon[[0,1,2]]{0,1,0};
        rvideal Thick[[0,1,2]] {
            rvmon[[0,1,2]]{2,0,0}, rvmon[[0,1,2]]{1,1,0},
            rvmon[[0,1,2]]{1,0,1}, rvmon[[0,1,2]]{0,2,0}
        };
        rvprojective DoublePoint[Thick];
        let a = rvTorClass(rvkterm[DoublePoint,[0]](one));
        let b = rvTorClass(rvkterm[DoublePoint,[1]](y));
        assert a == rvtorgen[DoublePoint,1,1,0];
        assert b == rvtorgen[DoublePoint,1,2,0];
        assert rvTorProd(a,b) == rvtorgen[DoublePoint,2,3,0];
        rvbetti[DoublePoint,2,3];
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvbetti[DoublePoint, 2, 3] = 1\n"
    );
}
