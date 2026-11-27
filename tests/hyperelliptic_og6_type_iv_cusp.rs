mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn setup() -> &'static str {
    include_str!("../examples/hyperelliptic-og6-type-iv-cusp.asl")
        .split_once("// The quotient")
        .unwrap()
        .0
}

#[test]
fn runnable_example_computes_the_cusp_gram_and_nontrivial_commutator() {
    assert_eq!(
        stdout(example("hyperelliptic-og6-type-iv-cusp")),
        concat!(
            "rvOG6CuspGram[P] = rank=3, signature=(0,3), determinant=-8, gram=[[-2, 0, 0], [0, -2, 0], [0, 0, -2]]\n",
            "rvOG6Pair = 0\n",
            "rvOG6Pair = -2\n",
        )
    );
}

#[test]
fn saturated_quotient_handles_a_divisibility_two_polarization_and_plane_shear() {
    let mut script = format!(
        r"{}
        let h2=add(add(h,h),delta);
        rvOG6Div(h2);
        rvog6cusp Q(h2,add(e,f),f);
        rvOG6CuspGram[Q];
        rvog6cusp R(add(hE,add(hF,hF)),e,f);
        rvOG6CuspGram[R];
        ",
        setup()
    );
    for i in 0..3 {
        writeln!(
            script,
            "let k{i}=rvOG6CuspBasis[Q,{i}]; rvOG6Pair(k{i},h2); rvOG6Pair(k{i},e); rvOG6Pair(k{i},f);"
        )
        .unwrap();
    }
    for i in 0..3 {
        for j in i..3 {
            writeln!(script, "rvOG6Pair(k{i},k{j});").unwrap();
        }
    }
    let output = stdout(interpret(&script));
    let mut lines = output.lines();
    assert_eq!(lines.next(), Some("rvOG6Div = 2"));
    assert!(lines.next().unwrap().contains("determinant=-6"));
    assert!(lines.next().unwrap().contains("determinant=-16"));
    for _ in 0..9 {
        assert_eq!(lines.next(), Some("rvOG6Pair = 0"));
    }
    let g: Vec<i128> = lines
        .map(|line| line.strip_prefix("rvOG6Pair = ").unwrap().parse().unwrap())
        .collect();
    let [a, b, c, d, e, f] = g[..] else {
        panic!("unexpected Gram entries")
    };
    assert_eq!(
        a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d),
        -6
    );
}

#[test]
fn integral_transvections_preserve_every_basis_pairing_and_fix_polarization() {
    let mut script = format!(
        r"{}
        let u=add(delta,f);
        let v=add(rvOG6ExcHalf[M],add(hE,neg(hF)));
        assert rvOG6CuspTransvection[P,0](u,h)==h;
        assert rvOG6CuspTransvection[P,1](v,h)==h;
        assert rvOG6CuspCentral[P,7](h)==h;
        assert rvOG6CuspTransvection[P,0](f,eDual)==add(eDual,f);
        assert rvOG6CuspTransvection[P,0](e,eDual)==eDual;
        assert rvOG6CuspTransvection[P,0](u,eDual)
            ==rvOG6CuspCentral[P,1](rvOG6CuspTransvection[P,0](delta,eDual));
        assert rvOG6CuspCentral[P,-3](fDual)==add(fDual,add(e,add(e,e)));
        ",
        setup()
    );
    for i in 0..8 {
        writeln!(
            script,
            r"
            let a{i}=rvOG6Basis[M,{i}];
            let t{i}=rvOG6CuspTransvection[P,0](u,a{i});
            let s{i}=rvOG6CuspTransvection[P,1](v,a{i});
            let c{i}=rvOG6CuspCentral[P,7](a{i});
            assert rvOG6CuspTransvection[P,0](neg(u),t{i})==a{i};
            assert rvOG6CuspTransvection[P,1](neg(v),s{i})==a{i};
            assert rvOG6CuspCentral[P,-7](c{i})==a{i};
            "
        )
        .unwrap();
    }
    for i in 0..8 {
        for j in i..8 {
            writeln!(script,
                "rvOG6Pair(a{i},a{j}); rvOG6Pair(t{i},t{j}); rvOG6Pair(s{i},s{j}); rvOG6Pair(c{i},c{j});"
            ).unwrap();
        }
    }
    let output = stdout(interpret(&script));
    let pairings: Vec<_> = output.lines().collect();
    assert_eq!(pairings.len(), 36 * 4);
    assert!(pairings
        .chunks_exact(4)
        .all(|row| row.iter().all(|x| *x == row[0])));
}

#[test]
fn commutator_and_center_work_in_the_declared_orientation_on_the_whole_basis() {
    let mut script = format!(
        r"{}
        let u=add(delta,f);
        let v=add(delta,rvOG6ExcHalf[M]);
        rvOG6Pair(u,v);
        rvog6cusp Swap(h,f,e);
        assert rvOG6CuspCentral[P,-2](eDual)==add(eDual,neg(add(f,f)));
        ",
        setup()
    );
    for i in 0..8 {
        writeln!(
            script,
            r"
            let x{i}=rvOG6Basis[M,{i}];
            let a{i}=rvOG6CuspTransvection[P,1](neg(v),x{i});
            let b{i}=rvOG6CuspTransvection[P,0](neg(u),a{i});
            let c{i}=rvOG6CuspTransvection[P,1](v,b{i});
            assert rvOG6CuspTransvection[P,0](u,c{i})==rvOG6CuspCentral[P,-2](x{i});
            assert rvOG6CuspCentral[Swap,7](x{i})==rvOG6CuspCentral[P,-7](x{i});
            assert rvOG6CuspCentral[P,5](rvOG6CuspCentral[P,-3](x{i}))==rvOG6CuspCentral[P,2](x{i});
            assert rvOG6CuspCentral[P,9](rvOG6CuspTransvection[P,0](u,x{i}))
                ==rvOG6CuspTransvection[P,0](u,rvOG6CuspCentral[P,9](x{i}));
            assert rvOG6CuspCentral[P,9](rvOG6CuspTransvection[P,1](v,x{i}))
                ==rvOG6CuspTransvection[P,1](v,rvOG6CuspCentral[P,9](x{i}));
            "
        )
        .unwrap();
    }
    assert_eq!(stdout(interpret(&script)), "rvOG6Pair = -2\n");
}

#[test]
fn plucker_minors_and_transvection_intermediates_allow_wide_cancellation() {
    let n = i128::MAX;
    assert_eq!(
        stdout(interpret(&format!(
            r"{}
            let bigE=rvOG6[M,0](rvAbMukai[B,0,0]{{{n},0,0,0,0,0}});
            let bigDual=rvOG6[M,0](rvAbMukai[B,0,0]{{0,0,0,0,0,{negative_n}}});
            assert rvOG6CuspTransvection[P,0](bigE,bigDual)==bigDual;
            let wideE=rvOG6[M,0](rvAbMukai[B,0,0]{{{n},{m},0,0,0,0}});
            let wideF=rvOG6[M,0](rvAbMukai[B,0,0]{{{m},{l},0,0,0,0}});
            rvog6cusp Wide(h,wideE,wideF);
            rvOG6CuspGram[Wide];
            assert rvOG6CuspCentral[Wide,1](eDual)==rvOG6CuspCentral[P,-1](eDual);
            assert rvOG6CuspCentral[P,{min}](e)==e;
            ",
            setup(),
            n = n, m = n - 1, l = n - 2, min = i128::MIN, negative_n = -n
        ))),
        "rvOG6CuspGram[Wide] = rank=3, signature=(0,3), determinant=-8, gram=[[-2, 0, 0], [0, -2, 0], [0, 0, -2]]\n"
    );
}

#[test]
fn rejects_nonprimitive_planes_incompatible_directions_and_real_overflow() {
    for (source, expected) in [
        ("rvog6cusp Bad(add(h,h),e,f);", "requires a primitive positive polarization"),
        ("rvog6cusp Bad(rvOG6ExcHalf[M],e,f);", "requires a primitive positive polarization"),
        ("rvog6cusp Bad(h,e,e);", "requires two independent isotropic generators"),
        ("rvog6cusp Bad(h,e,add(e,add(f,f)));", "requires a primitive isotropic plane"),
        ("rvog6cusp Bad(h,e,eDual);", "orthogonal isotropic generators"),
        ("rvog6cusp Bad(h,hE,f);", "perpendicular to the polarization"),
        ("rvog6cusp P(h,e,f);", "is already defined"),
        ("rvOG6CuspGram[Missing];", "unknown OG6 cusp"),
        ("show rvOG6CuspBasis[P,3];", "between 0 and 2"),
        ("show rvOG6CuspTransvection[P,2](delta,eDual);", "axis must be 0 or 1"),
        ("show rvOG6CuspTransvection[P,0](h,eDual);", "requires a direction perpendicular"),
        ("show rvOG6CuspTransvection[P,1](fDual,f);", "requires a direction perpendicular"),
        ("rvog6 N(rvAbMukai[B,1,-1]{0,0,0,0,0,0}); show rvOG6CuspCentral[P,1](zero[RVOG6(N)]);", "expected type RVOG6(M)"),
        ("rvog6 N(rvAbMukai[B,1,-1]{0,0,0,0,0,0}); rvog6cusp Bad(h,e,zero[RVOG6(N)]);", "expected type RVOG6(M)"),
        ("let x=add(eDual,f); show rvOG6CuspCentral[P,170141183460469231731687303715884105727](x);", "abelian Mukai coordinate exceeds"),
    ] {
        let error = diagnostic(interpret(&format!("{}{source}", setup())));
        assert!(error.contains(expected), "{error}");
    }
}
