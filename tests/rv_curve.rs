mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn sparse_nonclassical_orders_and_osculating_images_are_exact() {
    assert_eq!(
        stdout(example("rank-osculation")),
        concat!(
            "rvorders[Twist] = [0, 1, 4, 5]\n",
            "rvgroebner[Tangents] = [x6, x5 + x4, x3, x5^2 + x2*x7]\n",
            "rvgroebner[Planes] = [x9*x10 + x8*x11, x10^4 + x9*x11^3, x9^2*x11^2 + x8*x10^3, x9^3*x11 + x8^2*x10^2, x9^4 + x8^3*x10]\n",
            "rvorders[Large] = [0, 1, 1024, 1025]\n",
            "rvgroebner[LargeTangents] = [x6, x5 + x4, x3, x5^2 + x2*x7]\n",
        )
    );
}

#[test]
fn mobius_change_preserves_hasse_orders_and_transforms_the_minors() {
    let program = r"
        let a = rvmon[[0,1]]{1025,0};
        let b = rvmon[[0,1]]{1024,1};
        let c = rvmon[[0,1]]{1,1024};
        let d = rvmon[[0,1]]{0,1025};
        // Substitute s+t for s, keeping t unchanged.
        rvcurve Moved[[0,1]] { add(add(a,b),add(c,d)), add(b,d), add(c,d), d };
        rvorders[Moved];
        assert rvOscMinor[Moved,1,[0,1]] == add(rvmon[[0,1]]{2048,0},rvmon[[0,1]]{0,2048});
        assert rvOscMinor[Moved,1,[3,0]] == add(rvmon[[0,1]]{1024,1024},rvmon[[0,1]]{0,2048});
        assert rvOscMinor[Moved,3,[3,0,2,1]] == rvmon[[0,1]]{0,0};
        rvosculating Tangents[Moved,1,[2,3,4,5,6,7]];
        rvgroebner[Tangents];
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvorders[Moved] = [0, 1, 1024, 1025]\nrvgroebner[Tangents] = [x6, x5 + x4, x3, x5^2 + x2*x7]\n"
    );
}

#[test]
fn minors_keep_ramification_and_degenerate_curves_report_actual_rank() {
    let program = r"
        rvcurve Cusp[[0,1]] { rvmon[[0,1]]{3,0}, rvmon[[0,1]]{1,2}, rvmon[[0,1]]{0,3} };
        rvorders[Cusp];
        assert rvOscMinor[Cusp,1,[0,2]] == rvmon[[0,1]]{2,2};
        assert rvOscMinor[Cusp,1,[1,2]] == rvmon[[0,1]]{0,4};
        assert rvOscMinor[Cusp,2,[0,1,2]] == rvmon[[0,1]]{1,2};
        rvosculating Tangent[Cusp,1,[2,3,4]];
        rvgroebner[Tangent];
        rvosculating Whole[Cusp,2,[2]];
        rvgroebner[Whole];

        let s = rvmon[[0,1]]{1,0};
        let t = rvmon[[0,1]]{0,1};
        rvcurve Degenerate[[0,1]] { s,t,add(s,t) };
        rvorders[Degenerate];
        assert rvOscMinor[Degenerate,0,[2]] == add(s,t);
        rvosculating Line[Degenerate,0,[2,3,4]];
        rvgroebner[Line];
        rvosculating Point[Degenerate,1,[2,3,4]];
        rvgroebner[Point];
        rvhilbert[Point,1000];
    ";
    assert_eq!(
        stdout(interpret(program)),
        concat!(
            "rvorders[Cusp] = [0, 1, 2]\n",
            "rvgroebner[Tangent] = [x2]\n",
            "rvgroebner[Whole] = []\n",
            "rvorders[Degenerate] = [0, 1]\n",
            "rvgroebner[Line] = [x4 + x3 + x2]\n",
            "rvgroebner[Point] = [x4 + x3, x4 + x2]\n",
            "rvhilbert[Point, 1000] = 1\n",
        )
    );
}

#[test]
fn invalid_coordinates_basepoints_and_orders_are_diagnosed() {
    let setup = r"
        let s = rvmon[[0,1]]{1,0};
        let t = rvmon[[0,1]]{0,1};
        let g = add(add(rvmon[[0,1]]{2,0},rvmon[[0,1]]{1,1}),rvmon[[0,1]]{0,2});
        rvcurve C[[0,1]] {s,t};
    ";
    for (statement, reason) in [
        ("rvcurve E[[0,1]] {};", "at least one"),
        ("rvcurve E[[0]] {rvmon[[0]]{1}};", "exactly two"),
        (
            "rvcurve E[[0,1]] {rvmon[[0,1]]{0,0}};",
            "common positive degree",
        ),
        (
            "rvcurve E[[0,1]] {rvmon[[0,2]]{1,0}};",
            "different primitive selection",
        ),
        ("rvcurve E[[0,1]] {s,rvMul(s,t)};", "common positive degree"),
        (
            "rvcurve E[[0,1]] {rvMul(s,s),rvMul(s,t)};",
            "basepoint at [0:1]",
        ),
        (
            "rvcurve E[[0,1]] {rvMul(s,g),rvMul(t,g)};",
            "basepoints on the chart",
        ),
        ("rvcurve C[[0,1]] {s,t};", "already defined"),
        ("rvorders[Missing];", "unknown rational curve"),
        (
            "show rvOscMinor[C,2,[0,1,2]];",
            "out of range for coordinate rank 2",
        ),
        ("show rvOscMinor[C,1,[0]];", "exactly k+1"),
        (
            "show rvOscMinor[C,1,[0,0]];",
            "distinct coordinate positions",
        ),
        (
            "show rvOscMinor[C,1,[0,2]];",
            "coordinate position is out of range",
        ),
        ("rvosculating D[C,1,[2,3]];", "one new variable"),
        ("rvosculating D[C,1,[0]];", "fresh coordinate labels"),
    ] {
        let error = diagnostic(interpret(&format!("{setup}\n{statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
}

#[test]
fn huge_sparse_gaps_and_plucker_degree_checks_do_not_enumerate_degrees() {
    let q = 1usize << (usize::BITS / 2);
    let degree = usize::MAX;
    let program = format!(
        "rvcurve Huge[[0,1]] {{rvmon[[0,1]]{{{},0}},rvmon[[0,1]]{{{q},1}},
          rvmon[[0,1]]{{1,{q}}},rvmon[[0,1]]{{0,{}}}}};
         rvorders[Huge];
         assert rvOscMinor[Huge,3,[0,1,2,3]] == rvmon[[0,1]]{{0,0}};",
        q + 1,
        q + 1,
    );
    assert_eq!(
        stdout(interpret(&program)),
        format!("rvorders[Huge] = [0, 1, {q}, {}]\n", q + 1)
    );

    let overflow = format!(
        "rvcurve Huge[[0,1]] {{rvmon[[0,1]]{{{degree},0}},
          add(rvmon[[0,1]]{{0,{degree}}},rvmon[[0,1]]{{{},1}})}};
         show rvOscMinor[Huge,1,[0,1]];",
        degree - 1,
    );
    let error = diagnostic(interpret(&overflow));
    assert!(
        error.contains("osculating Plucker degree is out of range"),
        "{error}"
    );
}
