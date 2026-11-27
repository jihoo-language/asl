mod common;

use common::{diagnostic, example, interpret, stdout};

fn setup() -> &'static str {
    include_str!("../examples/hyperelliptic-og6-monodromy.asl")
        .split_once("rvOG6Orbit(delta)")
        .unwrap()
        .0
}

#[test]
fn example_distinguishes_stable_and_full_monodromy_and_wall_types() {
    assert_eq!(
        stdout(example("hyperelliptic-og6-monodromy")),
        concat!(
            "rvOG6Orbit = content=1, square=-2, div=2, residue=[1,0], discq=-1/2\n",
            "rvOG6Orbit = content=1, square=-2, div=2, residue=[0,1], discq=-1/2\n",
            "rvOG6SameOrbit = stable=false, full=true\n",
            "rvOG6SameOrbit = stable=false, full=false\n",
            "rvOG6Wall = content=1, square=-4, div=2, residue=[1,1], discq=-1, kind=divisorial\n",
            "rvOG6Wall = content=2, square=-2, div=2, residue=[0,1], discq=-1/2, kind=divisorial\n",
            "rvOG6Wall = content=1, square=-2, div=1, residue=[0,0], discq=0, kind=small\n",
            "rvOG6Orbit = content=1, square=0, div=1, residue=[0,0], discq=0\n",
        )
    );
}

#[test]
fn residue_is_read_from_the_actual_perpendicular_mod_two() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{}
            let other=rvAbMukai[B,2,3]{{2,3,0,0,5,4}};
            rvog6 N(other);
            // g = 2b - w for b = (-1, wedge23, 0).
            let g=rvOG6[N,0](rvAbMukai[B,-4,-3]{{-2,-3,0,0,-5,-2}});
            let err=rvOG6ExcHalf[N];
            let s=add(g,err);
            let u=rvOG6[N,0](rvAbMukai[B,0,0]{{0,0,1,0,0,0}});
            let v=rvOG6[N,0](rvAbMukai[B,0,0]{{0,0,0,-1,0,0}});
            let iso=add(add(u,v),err);
            rvOG6Orbit(g);
            rvOG6Orbit(s);
            rvOG6Orbit(add(g,u));
            rvOG6Orbit(iso);
            rvOG6Orbit(add(iso,iso));
            rvOG6SameOrbit(iso,u);
            rvOG6SameOrbit(g,rvOG6Reflect(s,g));
            rvOG6SameOrbit(s,rvOG6Reflect(g,s));
            rvOG6Wall(add(g,neg(err)));
            assert rvOG6Reflect(s,g)==neg(err);
            assert rvOG6Reflect(s,err)==neg(g);
            ",
            setup()
        ))),
        concat!(
            "rvOG6Orbit = content=1, square=-2, div=2, residue=[1,0], discq=-1/2\n",
            "rvOG6Orbit = content=1, square=-4, div=2, residue=[1,1], discq=-1\n",
            "rvOG6Orbit = content=1, square=-2, div=1, residue=[0,0], discq=0\n",
            "rvOG6Orbit = content=1, square=0, div=1, residue=[0,0], discq=0\n",
            "rvOG6Orbit = content=2, square=0, div=1, residue=[0,0], discq=0\n",
            "rvOG6SameOrbit = stable=true, full=true\n",
            "rvOG6SameOrbit = stable=false, full=true\n",
            "rvOG6SameOrbit = stable=true, full=true\n",
            "rvOG6Wall = content=1, square=-4, div=2, residue=[1,1], discq=-1, kind=divisorial\n",
        )
    );
}

#[test]
fn reflection_is_an_integral_involution_and_preserves_pairing() {
    let output = stdout(interpret(&format!(
        r"{}
        let theta=rvOG6[M,0](rvAbMukai[B,0,0]{{0,1,0,0,1,0}});
        let x=add(add(e,theta),delta);
        let y=add(add(e,e),add(delta,theta));
        let rx=rvOG6Reflect(add(exchange,exchange),x);
        let ry=rvOG6Reflect(neg(exchange),y);
        assert rvOG6Reflect(exchange,rx)==x;
        assert rvOG6Reflect(exchange,ry)==y;
        assert rvOG6Reflect(delta,delta)==neg(delta);
        assert rvOG6Reflect(exchange,theta)==theta;
        rvOG6Pair(x,y);
        rvOG6Pair(rx,ry);
        rvOG6Pair(x,x);
        rvOG6Pair(rx,rx);
        rvOG6SameOrbit(x,rx);
        rvOG6SameOrbit(e,rvOG6Reflect(e,e));
        rvOG6Wall(rvOG6[M,0](rvAbMukai[B,0,0]{{0,1,0,0,-2,0}}));
        ",
        setup()
    )));
    assert_eq!(
        output,
        concat!(
            "rvOG6Pair = -4\n",
            "rvOG6Pair = -4\n",
            "rvOG6Pair = -2\n",
            "rvOG6Pair = -2\n",
            "rvOG6SameOrbit = stable=true, full=true\n",
            "rvOG6SameOrbit = stable=true, full=true\n",
            "rvOG6Wall = content=1, square=-4, div=1, residue=[0,0], discq=0, kind=regular\n",
        )
    );
}

#[test]
fn wide_norm_and_content_and_reflection_cancel_without_narrowing() {
    let m = i128::MAX;
    let min = i128::MIN;
    let square = "57896044618658097711785492504343953925614145232057466629338668181661260185604";
    assert_eq!(
        stdout(interpret(&format!(
            r"{}
            let x=rvOG6[M,{m}](rvAbMukai[B,{m},{m}]{{0,0,0,0,0,0}});
            let x2=rvOG6Reflect(exchange,x);
            assert x2==neg(x);
            rvOG6Orbit(x);
            rvOG6SameOrbit(x,x2);
            let giant=rvOG6[M,0](rvAbMukai[B,0,0]{{0,{m},0,0,{},0}});
            let negativeSquare=rvOG6[M,0](rvAbMukai[B,0,0]{{0,{m},0,0,{},0}});
            rvOG6Orbit(giant);
            rvOG6SameOrbit(giant,neg(giant));
            rvOG6SameOrbit(giant,negativeSquare);
            let minimum=rvOG6[M,{min}](zero[RVAbMukai(B)]);
            rvOG6Orbit(minimum);
            rvOG6Orbit(zero[RVOG6(M)]);
            rvOG6SameOrbit(zero[RVOG6(M)],zero[RVOG6(M)]);
            ",
            setup(),
            m - 1,
            1 - m
        ))),
        format!(
            concat!(
                "rvOG6Orbit = content={m}, square=-4, div=2, residue=[1,1], discq=-1\n",
                "rvOG6SameOrbit = stable=true, full=true\n",
                "rvOG6Orbit = content=1, square={square}, div=1, residue=[0,0], discq=0\n",
                "rvOG6SameOrbit = stable=true, full=true\n",
                "rvOG6SameOrbit = stable=false, full=false\n",
                "rvOG6Orbit = content={absmin}, square=-2, div=2, residue=[0,1], discq=-1/2\n",
                "rvOG6Orbit = content=0, square=0, div=0, residue=[0,0], discq=0\n",
                "rvOG6SameOrbit = stable=true, full=true\n",
            ),
            m = m,
            square = square,
            absmin = min.unsigned_abs()
        )
    );
}

#[test]
fn nominal_wrong_wall_and_true_output_range_errors_are_diagnosed() {
    for (source, expected) in [
        ("rvOG6Orbit(w);", "an RVOG6 class"),
        ("rvOG6Wall(zero[RVOG6(M)]);", "requires a nonzero ray"),
        ("show rvOG6Reflect(small,e);", "requires a divisorial OG6 wall ray"),
        ("show rvOG6Reflect(zero[RVOG6(M)],e);", "requires a divisorial OG6 wall ray"),
        ("rvog6 Twin(w); rvOG6SameOrbit(e,rvOG6ExcHalf[Twin]);", "expected type RVOG6(M)"),
        ("rvog6 Twin(w); show rvOG6Reflect(e,rvOG6ExcHalf[Twin]);", "expected type RVOG6(M)"),
        (
            "let minimum=rvOG6[M,-170141183460469231731687303715884105728](zero[RVAbMukai(B)]); show rvOG6Reflect(e,minimum);",
            "half-exceptional coefficient exceeds",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{}{source}", setup())));
        assert!(error.contains(expected), "{error}");
    }
}
