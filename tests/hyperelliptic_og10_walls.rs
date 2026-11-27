mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn setup() -> &'static str {
    include_str!("../examples/hyperelliptic-og10-walls.asl")
        .split_once("// The numerical classes")
        .unwrap()
        .0
}

const INTEGRAL_SEGMENT: &str = r"
    let h0=add(add(H,H),sigma);
    let h1=add(add(H,H),neg(sigma));
";

#[test]
fn example_matches_the_erratum_and_finds_a_divisorial_crossing() {
    assert_eq!(
        stdout(example("hyperelliptic-og10-walls")),
        concat!(
            "rvOG10Wall = content=1, square=-10, div=1, kind=regular\n",
            "rvOG10Wall = content=1, square=-4, div=1, kind=small codimension 5\n",
            "rvOG10Wall = content=1, square=-24, div=3, kind=small codimension 3\n",
            "rvOG10Wall = content=2, square=-6, div=3, kind=divisorial\n",
            "rvOG10Walls = count=5, crossings=[{t=0, square=-6, div=3, kind=divisorial}, {t=1/6, square=-4, div=1, kind=small codimension 5}, {t=1/2, square=-6, div=3, kind=divisorial}, {t=5/6, square=-4, div=1, kind=small codimension 5}, {t=1, square=-6, div=3, kind=divisorial}]\n",
            "rvOG10Walls = count=3, crossings=[{t=0, square=-6, div=3, kind=divisorial}, {t=1/2, square=-4, div=1, kind=small codimension 5}, {t=1, square=-6, div=3, kind=divisorial}]\n",
        )
    );
}

#[test]
fn primitive_rays_in_the_saturated_plane_and_endpoints_are_oriented_and_sorted() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{}{INTEGRAL_SEGMENT}
            let h0x=add(h0,h0);
            let h1x=add(h1,add(h1,h1));
            let initial=rvOG10WallAt[0](h0x,h1x);
            let a=rvOG10WallAt[1](h0x,h1x);
            let b=rvOG10WallAt[2](h0x,h1x);
            let c=rvOG10WallAt[3](h0x,h1x);
            let final=rvOG10WallAt[4](h0x,h1x);
            let h3=add(H,add(H,H));
            let s2=add(sigma,sigma);
            assert initial==neg(add(h3,s2));
            assert final==add(h3,neg(s2));
            assert a==neg(add(H,sigma));
            assert b==neg(sigma);
            assert c==add(H,neg(sigma));
            assert initial==rvOG10WallAt[0](h0,h1);
            assert a==rvOG10WallAt[1](h0,h1);
            assert b==rvOG10WallAt[2](h0,h1);
            assert c==rvOG10WallAt[3](h0,h1);
            assert final==rvOG10WallAt[4](h0,h1);
            assert neg(initial)==rvOG10WallAt[4](h1,h0);
            assert neg(a)==rvOG10WallAt[3](h1,h0);
            assert neg(b)==rvOG10WallAt[2](h1,h0);
            assert neg(c)==rvOG10WallAt[1](h1,h0);
            assert neg(final)==rvOG10WallAt[0](h1,h0);
            assert rvOG10WallAt[0](H,h0)==sigma;
            assert rvOG10WallAt[2](h0,H)==neg(sigma);
            rvOG10Walls(h0x,h1x);
            rvOG10Walls(h1x,h0x);
            let h2=add(add(H,H),h0);
            rvOG10Walls(h2,add(h2,H));
            ",
            setup()
        ))),
        concat!(
            "rvOG10Walls = count=5, crossings=[{t=0, square=-6, div=3, kind=divisorial}, {t=2/17, square=-4, div=1, kind=small codimension 5}, {t=2/5, square=-6, div=3, kind=divisorial}, {t=10/13, square=-4, div=1, kind=small codimension 5}, {t=1, square=-6, div=3, kind=divisorial}]\n",
            "rvOG10Walls = count=5, crossings=[{t=0, square=-6, div=3, kind=divisorial}, {t=3/13, square=-4, div=1, kind=small codimension 5}, {t=3/5, square=-6, div=3, kind=divisorial}, {t=15/17, square=-4, div=1, kind=small codimension 5}, {t=1, square=-6, div=3, kind=divisorial}]\n",
            "rvOG10Walls = count=0, crossings=[]\n",
        )
    );
}

#[test]
fn the_small_codimension_three_wall_can_be_an_actual_segment_crossing() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{}
            let E=rvOG10[M,-1](rvMukai[2,2](Hk3));
            let F=add(add(E,add(E,E)),sigma);
            let J=rvOG10[M,0](rvMukai[0,0](rvK3Lattice[B]{{[0,0,2,-2,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}}));
            let j3=add(J,add(J,J));
            let h0=add(j3,F);
            let h1=add(j3,neg(F));
            rvOG10Walls(h0,h1);
            assert rvOG10WallAt[0](h0,h1)==neg(F);
            assert rvOG10WallAt[0](h1,h0)==F;
            ",
            setup()
        ))),
        "rvOG10Walls = count=1, crossings=[{t=1/2, square=-24, div=3, kind=small codimension 3}]\n"
    );
}

#[test]
fn divisorial_reflections_use_primitive_content_and_preserve_the_form() {
    let output = stdout(interpret(&format!(
        r"{}{INTEGRAL_SEGMENT}
        let root=rvOG10[M,0](rvMukai[0,0](rvK3FromKummerDiv(rvKummerNodeDiv[B,0])));
        let h=rvOG10Reflect(root,H);
        let t=rvOG10Reflect(add(sigma,sigma),h0);
        assert rvOG10Reflect(root,h)==H;
        assert h==add(H,neg(add(root,root)));
        assert rvOG10Reflect(neg(sigma),t)==h0;
        assert t==h1;
        rvOG10Wall(root);
        rvOG10Wall(neg(add(root,root)));
        rvOG10Pair(h,h);
        rvOG10Pair(t,h0);
        rvOG10Pair(h1,h0);
        ",
        setup()
    )));
    assert_eq!(
        output,
        concat!(
            "rvOG10Wall = content=1, square=-2, div=1, kind=divisorial\n",
            "rvOG10Wall = content=2, square=-2, div=1, kind=divisorial\n",
            "rvOG10Pair = 2\n",
            "rvOG10Pair = 14\n",
            "rvOG10Pair = 14\n",
        )
    );
}

#[test]
fn exact_wide_square_roots_and_segment_parameters_survive_large_coordinates() {
    let b = 1i128 << 124;
    let a = 2 * b + 1;
    let output = stdout(interpret(&format!(
        r"{}
        let far=rvOG10[M,{}](rvMukai[-{b},-{b}](rvK3Lattice[B]{{[0,{},0,0,{},0],[{},0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}}));
        rvOG10Wall(far);
        rvOG10Walls(H,far);
        assert rvOG10WallAt[0](H,far)==sigma;
        assert rvOG10WallAt[1](H,far)==add(H,sigma);
        rvOG10Walls(far,H);
        assert rvOG10WallAt[0](far,H)==neg(add(H,sigma));
        assert rvOG10WallAt[1](far,H)==neg(sigma);
        ",
        setup(),
        2 * b,
        2 * a,
        2 * a,
        2 * a
    )));
    let square = "904625697166532776746648320380374280273812938660786138289949678777705431042";
    assert_eq!(
        output,
        format!(
            concat!(
                "rvOG10Wall = content=1, square={square}, div=3, kind=regular\n",
                "rvOG10Walls = count=2, crossings=[{{t=0, square=-6, div=3, kind=divisorial}}, {{t=1/{b}, square=-4, div=1, kind=small codimension 5}}]\n",
                "rvOG10Walls = count=2, crossings=[{{t={before_end}/{b}, square=-4, div=1, kind=small codimension 5}}, {{t=1, square=-6, div=3, kind=divisorial}}]\n",
            ),
            square = square,
            b = b,
            before_end = b - 1,

        )
    );
}

#[test]
fn invalid_cones_small_reflections_and_search_bounds_are_diagnosed() {
    let setup = format!("{}{INTEGRAL_SEGMENT}", setup());
    for (source, expected) in [
        ("rvOG10Wall(zero[RVOG10(M)]);", "require a nonzero class"),
        ("rvOG10Walls(sigma,H);", "requires positive endpoints in the same cone component"),
        ("rvOG10Walls(H,neg(h0));", "requires positive endpoints in the same cone component"),
        ("rvOG10Walls(H,add(H,H));", "independent endpoints spanning a plane"),
        ("let j=rvOG10[M,0](rvMukai[0,0](rvK3Lattice[B]{[0,0,2,-2,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]})); rvOG10Walls(H,add(H,j));", "independent endpoints spanning a plane"),
        ("show rvOG10WallAt[5](h0,h1);", "index 5 is out of range"),
        ("show rvOG10Reflect(add(H,sigma),h0);", "requires a divisorial wall ray"),
        ("let D=rvOG10[M,0](rvMukai[3,3](add(Hk3,Hk3))); show rvOG10Reflect(D,h0);", "requires a divisorial wall ray"),
        ("rvmukaimoduli Other(rvMukaiIdeal[B,2]); rvOG10Walls(H,rvOG10Exc[Other]);", "expected type RVOG10(M)"),
        ("rvmukaimoduli Other(rvMukaiIdeal[B,2]); show rvOG10Reflect(sigma,rvOG10Exc[Other]);", "expected type RVOG10(M)"),
    ] {
        let output = diagnostic(interpret(&format!("{setup}{source}")));
        assert!(output.contains(expected), "{output}");
    }

    let mut source = setup;
    let (mut a, mut b) = (1i128, 0i128);
    for _ in 0..21 {
        (a, b) = (2 * a + 3 * b, a + 2 * b);
    }
    writeln!(
        source,
        "let far=rvOG10[M,{}](rvMukai[-{b},-{b}](rvK3Lattice[B]{{[0,{},0,0,{},0],[{},0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}})); rvOG10Walls(H,far);",
        2*b, 2*a, 2*a, 2*a
    ).unwrap();
    let error = diagnostic(interpret(&source));
    assert!(
        error.contains("OG10 wall segment search exceeds the work limit"),
        "{error}"
    );
}
