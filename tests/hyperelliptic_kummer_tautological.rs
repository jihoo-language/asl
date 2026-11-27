mod common;

use common::{diagnostic, example, interpret, stdout};

const MARKED: &str = r"
    rvunram W[3,1]{1,0,1};
    let z=rvUnramElem[W]{};
    let o=rvUnramElem[W]{1};
    let i=rvUnramElem[W]{0,1};
    rvhyperelliptic C[W]{z,neg(o),z,z,z,o};
    let e1=rvJacPoint[C](z,z);
    let f1=rvJacPoint[C](o,z);
    let e2=add(add(e1,f1),rvJacPoint[C](neg(o),z));
    let f2=add(add(e1,f1),rvJacPoint[C](i,z));
    rvjactorsionbasis B[C,2](neg(o)){e1,e2,f1,f2};
";

#[test]
fn example_computes_weight_three_top_segre_independently_by_wreath_cup() {
    assert_eq!(
        stdout(example("hyperelliptic-kummer-tautological")),
        "rvHilbTautChernDen[3, 2] = 3\nrvHilbTautSegre[3] = -160\n"
    );
}

#[test]
fn c0_c1_and_c2_agree_with_the_vacuum_surface_and_hilbert_square() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
        let h=rvK3LatticeBasis[B,0];
        let empty=rvHilbTautChern[0,0,1](h);
        assert empty==rvNakajimaVac[B];
        assert rvHilbTautChern[0,1,1](h)==zero[RVNakajima(B,0,2)];
        assert rvHilbTautChern[1,0,1](h)==rvNakajimaMonomial[B]{{(1,0,1)}};
        assert rvHilbTautChern[1,1,1](h)==rvNakajimaMonomial[B]{{(1,1,1)}};
        assert rvHilbTautChern[2,0,2](h)==rvNakajimaMonomial[B]{{(1,0,2)}};
        let x=rvHilbTautChern[2,1,2](h);
        let e=rvHilbSquareBasis[B,2,0];
        let delta=rvHilbSquareBasis[B,2,22];
        let expected=add(rvHilbSquareScale[2,1](e),rvHilbSquareScale[-2,1](delta));
        assert rvNakajimaHilbSquare(x)==expected;
        let y=rvHilbTautChern[2,2,2](h);
        assert y==add(rvNakajimaMonomial[B]{{(1,1,2)}},neg(rvNakajimaMonomial[B]{{(2,1,1)}}));
        rvHilbTautChernDen[0,0](h);
        rvHilbTautChernDen[1,1](h);
        rvHilbTautChernDen[2,1](h);
        rvHilbTautChernDen[3,0](h);
        rvHilbTautSegre[0](h);
        rvHilbTautSegre[1](h);
    "
        ))),
        concat!(
            "rvHilbTautChernDen[0, 0] = 1\n",
            "rvHilbTautChernDen[1, 1] = 1\n",
            "rvHilbTautChernDen[2, 1] = 2\n",
            "rvHilbTautChernDen[3, 0] = 6\n",
            "rvHilbTautSegre[0] = 1\n",
            "rvHilbTautSegre[1] = -2\n"
        )
    );
}

#[test]
fn scaled_and_mixed_colors_have_their_reduced_denominator() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
        let h=rvK3LatticeBasis[B,0];
        let g=rvK3LatticeBasis[B,1];
        let twice=add(h,h);
        let mixed=add(h,neg(g));
        assert rvHilbTautChern[2,2,1](twice)==add(rvNakajimaScale[2](rvNakajimaMonomial[B]{{(1,1,2)}}),neg(rvNakajimaMonomial[B]{{(2,1,1)}}));
        let a=rvNakajimaMonomial[B]{{(1,1,2)}};
        let b=rvNakajimaMonomial[B]{{(1,2,2)}};
        let c=rvNakajimaScale[-2](rvNakajimaMonomial[B]{{(1,1,1),(1,2,1)}});
        let d=neg(rvNakajimaMonomial[B]{{(2,1,1)}});
        let e=rvNakajimaMonomial[B]{{(2,2,1)}};
        assert rvHilbTautChern[2,2,2](mixed)==add(add(add(a,b),add(c,d)),e);
        assert rvHilbTautChern[2,2,-2](mixed)==neg(rvHilbTautChern[2,2,2](mixed));
        rvHilbTautChernDen[2,2](twice);
        rvHilbTautChernDen[2,2](mixed);
        rvHilbTautChernDen[3,3](h);
        rvHilbTautChernDen[3,3](zero[RVK3Lattice(B)]);
        rvHilbTautChernDen[3,4](h);
    "
        ))),
        concat!(
            "rvHilbTautChernDen[2, 2] = 1\n",
            "rvHilbTautChernDen[2, 2] = 2\n",
            "rvHilbTautChernDen[3, 3] = 6\n",
            "rvHilbTautChernDen[3, 3] = 1\n",
            "rvHilbTautChernDen[3, 4] = 1\n"
        )
    );
}

#[test]
fn weight_two_top_segre_of_two_distinct_lines_matches_wreath_cup() {
    let mut source =
        format!("{MARKED} let trivial=zero[RVK3Lattice(B)]; let h=rvK3LatticeBasis[B,0];");
    for (i, (line, segre)) in [("trivial", 12i128), ("h", 24)].into_iter().enumerate() {
        source.push_str(&format!(
            r"
            let x{i}=rvHilbTautChern[2,1,2]({line});
            let y{i}=rvHilbTautChern[2,2,2]({line});
            let x2n{i}=rvNakajimaCup(x{i},x{i});
            let x4n{i}=rvNakajimaCup(x2n{i},x2n{i});
            let mixed{i}=rvNakajimaCup(x2n{i},y{i});
            let y2n{i}=rvNakajimaCup(y{i},y{i});
            let scaled{i}=add(x4n{i},add(rvNakajimaScale[-6](mixed{i}),rvNakajimaScale[4](y2n{i})));
            assert scaled{i}==rvNakajimaScale[{}](rvNakajimaMonomial[B]{{(1,23,2)}});
            rvHilbTautSegre[2]({line});
        ",
            16 * segre
        ));
    }
    assert_eq!(
        stdout(interpret(&source)),
        "rvHilbTautSegre[2] = 12\nrvHilbTautSegre[2] = 24\n"
    );
}

#[test]
fn four_point_partitions_and_large_cancelled_numerators_have_exact_denominators() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
        let trivial=zero[RVK3Lattice(B)];
        let c2=rvHilbTautChern[4,2,24](trivial);
        let a=rvNakajimaScale[8](rvNakajimaMonomial[B]{{(1,0,1),(3,0,1)}});
        let b=rvNakajimaScale[3](rvNakajimaMonomial[B]{{(2,0,2)}});
        assert c2==add(a,b);
        rvHilbTautChernDen[4,2](trivial);
        let huge=rvK3Lattice[B]{{[0,12000000000000000000000000000000,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}};
        rvHilbTautChernDen[3,3](huge);
    "
        ))),
        "rvHilbTautChernDen[4, 2] = 24\nrvHilbTautChernDen[3, 3] = 1\n"
    );
}

#[test]
fn weight_three_nonzero_line_uses_all_three_chern_classes_in_segre_cup() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{MARKED}
        let h=rvK3LatticeBasis[B,0];
        let x=rvHilbTautChern[3,1,2](h);
        let y=rvHilbTautChern[3,2,6](h);
        let q=rvHilbTautChern[3,3,6](h);
        let x2=rvNakajimaCup(x,x);
        let x3=rvNakajimaCup(x2,x);
        let x4=rvNakajimaCup(x2,x2);
        let x6=rvNakajimaCup(x4,x2);
        let y2=rvNakajimaCup(y,y);
        let y3=rvNakajimaCup(y2,y);
        let q2=rvNakajimaCup(q,q);
        let part1=add(rvNakajimaScale[27](x6),rvNakajimaScale[-90](rvNakajimaCup(x4,y)));
        let part2=add(rvNakajimaScale[72](rvNakajimaCup(x2,y2)),rvNakajimaScale[-8](y3));
        let part3=add(rvNakajimaScale[144](rvNakajimaCup(x3,q)),rvNakajimaScale[-144](rvNakajimaCup(x,rvNakajimaCup(y,q))));
        let scaled=add(add(part1,part2),add(part3,rvNakajimaScale[48](q2)));
        assert scaled==rvNakajimaScale[-483840](rvNakajimaMonomial[B]{{(1,23,3)}});
        rvHilbTautSegre[3](h);
    "
        ))),
        "rvHilbTautSegre[3] = -280\n"
    );
}

#[test]
fn enormous_formal_grades_and_signed_segre_edge_are_exact() {
    let n = 1_000_000_000_000_000_000_000_000u128;
    assert_eq!(stdout(interpret(&format!(r"{MARKED}
        let l=zero[RVK3Lattice(B)];
        let h=rvK3LatticeBasis[B,0];
        assert rvHilbTautChern[{n},{j},{n}](l)==neg(rvNakajimaMonomial[B]{{({n},0,1)}});
        assert rvHilbTautChern[{n},{high},1](h)==zero[RVNakajima(B,{n},{degree})];
        assert rvHilbTautChern[{n},0,0](h)==zero[RVNakajima(B,{n},0)];
        rvHilbTautChernDen[{n},{j}](l);
        rvHilbTautChernDen[{n},{n}](l);
        let positive=rvK3Lattice[B]{{[0,2,0,0,{w},0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}};
        rvHilbTautSegre[{n}](positive);
        let edge=rvK3Lattice[B]{{[0,2,0,0,252,0],[2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}};
        rvHilbTautSegre[127](edge);
    ", j=n-1, high=n+1, degree=2*(n+1), w=2*(n-1)))), format!("rvHilbTautChernDen[{n}, {}] = {n}\nrvHilbTautChernDen[{n}, {n}] = 1\nrvHilbTautSegre[{n}] = 0\nrvHilbTautSegre[127] = {}\n",n-1,i128::MIN));
}

#[test]
fn limits_and_insufficient_scale_are_reported() {
    for (input, expected) in [
        (
            "show rvHilbTautChern[2,1,1](rvK3LatticeBasis[B,0]);",
            "scale does not clear the denominator",
        ),
        (
            "show rvHilbTautChern[2,2,1](rvK3LatticeBasis[B,0]);",
            "scale does not clear the denominator",
        ),
        (
            "rvHilbTautChernDen[1000000000000000000000,0](zero[RVK3Lattice(B)]);",
            "denominator exceeds the unsigned 128-bit range",
        ),
        (
            "rvHilbTautChernDen[1000000000,1000000000](rvK3LatticeBasis[B,0]);",
            "partition enumeration exceeds the computation limit",
        ),
        (
            "rvHilbTautSegre[1000000000](zero[RVK3Lattice(B)]);",
            "Segre number exceeds the signed 128-bit range",
        ),
        ("rvHilbTautSegre[1](o);", "an RVK3Lattice vector"),
    ] {
        let error = diagnostic(interpret(&format!("{MARKED}{input}")));
        assert!(error.contains(expected), "{error}");
    }
}
