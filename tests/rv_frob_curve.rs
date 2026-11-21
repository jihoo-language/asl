mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn ballico_hefez_has_the_exact_frobenius_gap_and_ramification_divisor() {
    assert_eq!(
        stdout(example("rank-frobenius")),
        concat!(
            "rvorders[BH] = [0, 1, 4]\n",
            "rvFrobOrders[BH, 1] = [0, 1]\n",
            "rvFrobOrders[BH, 2] = [0, 4]\n",
            "rvFrobOrders[BH, 4] = [0, 1]\n",
            "RVPoly([0, 1], 22){x0^2*x1^20 + x0^5*x1^17 + x0^17*x1^5 + x0^20*x1^2}\n",
            "rvOrdAt[7, [1, 0]] = 2\n",
            "rvOrdAt[7, [1, 1]] = 2\n",
            "rvOrdAt[7, [1, 2]] = 2\n",
            "rvOrdAt[7, [1, 3]] = 2\n",
            "rvOrdAt[7, [0, 1]] = 2\n",
            "rvOrdAt[19, [1, 2]] = 1\n",
        )
    );
}

#[test]
fn vanishing_orders_cover_each_f16_point_and_are_scaling_invariant() {
    let mut program = include_str!("../examples/rank-frobenius.asl").to_owned();
    let mut expected = stdout(example("rank-frobenius"));
    for t in 0..16 {
        let scaled = if t < 8 { t << 1 } else { (t << 1) ^ 19 };
        writeln!(program, "rvOrdAt[19,[2,{scaled}]](w);").unwrap();
        let order = if [0, 1, 6, 7].contains(&t) { 2 } else { 1 };
        writeln!(expected, "rvOrdAt[19, [2, {scaled}]] = {order}").unwrap();
    }
    program.push_str("rvOrdAt[19,[0,11]](w);\n");
    expected.push_str("rvOrdAt[19, [0, 11]] = 2\n");
    assert_eq!(stdout(interpret(&program)), expected);
}

#[test]
fn a_large_segre_wronskian_stays_sparse_and_degenerate_coordinates_are_skipped() {
    let h = 1024;
    let e = 18;
    let q = 1 << e;
    let program = format!(
        "let a = rvmon[[0,1]]{{{},0}};
         let b = rvmon[[0,1]]{{{h},1}};
         let c = rvmon[[0,1]]{{1,{h}}};
         let d = rvmon[[0,1]]{{0,{}}};
         rvcurve Segre[[0,1]] {{a,b,c,d}};
         rvcurve Dependent[[0,1]] {{zero[RVPoly([0,1],{})], add(a,c), d, add(a,c), b, c}};
         rvFrobOrders[Segre,{e}];
         rvFrobOrders[Dependent,{e}];
         let f = add(rvmon[[0,1]]{{1,{q}}},rvmon[[0,1]]{{{q},1}});
         let f_h = add(rvmon[[0,1]]{{{h},{}}},rvmon[[0,1]]{{{},{h}}});
         let w = rvFrobWronskian[Segre,{e}];
         assert w == rvMul(f,f_h);
         assert rvFrobWronskian[Dependent,{e}] == w;
         rvOrdAt[7,[0,1]](w);
         rvOrdAt[7,[1,2]](w);",
        h + 1,
        h + 1,
        h + 1,
        h * q,
        h * q,
    );
    assert_eq!(
        stdout(interpret(&program)),
        concat!(
            "rvFrobOrders[Segre, 18] = [0, 1, 1024]\n",
            "rvFrobOrders[Dependent, 18] = [0, 1, 1024]\n",
            "rvOrdAt[7, [0, 1]] = 1025\n",
            "rvOrdAt[7, [1, 2]] = 1025\n",
        )
    );
}

#[test]
fn the_local_binary_split_skips_huge_cancelled_orders() {
    let h = 1usize << (usize::BITS / 2);
    let program = format!(
        "let f = add(rvmon[[0,1]]{{{h},{}}},rvmon[[0,1]]{{0,{}}});
         rvOrdAt[3,[1,1]](f);
         rvOrdAt[3,[1,0]](f);
         rvOrdAt[3,[0,1]](f);
         rvOrdAt[7,[1,2]](f);
         rvOrdAt[7,[1,2]](zero[RVPoly([0,1],12)]);
         rvOrdAt[7,[1,2]](rvmon[[0,1]]{{0,0}});",
        h - 1,
        2 * h - 1,
    );
    assert_eq!(
        stdout(interpret(&program)),
        format!("rvOrdAt[3, [1, 1]] = {h}\nrvOrdAt[3, [1, 0]] = {}\nrvOrdAt[3, [0, 1]] = 0\nrvOrdAt[7, [1, 2]] = 0\nrvOrdAt[7, [1, 2]] = infinite\nrvOrdAt[7, [1, 2]] = 0\n", h - 1)
    );
}

#[test]
fn bad_curves_exponents_fields_and_projective_points_report_diagnostics() {
    let setup = r"
        let s = rvmon[[0,1]]{1,0};
        let t = rvmon[[0,1]]{0,1};
        rvcurve C[[0,1]] {s,t};
    ";
    for (statement, reason) in [
        ("rvFrobOrders[Missing,2];", "unknown rational curve"),
        ("rvFrobOrders[C,0];", "Frobenius extension exponent"),
        ("show rvFrobWronskian[C,0];", "Frobenius extension exponent"),
        (
            "rvcurve Basepoint[[0,1]] {rvMul(s,s),rvMul(s,t)};",
            "common projective basepoint",
        ),
        ("rvOrdAt[5,[1,0]](s);", "reducible"),
        ("rvOrdAt[3,[0,0]](zero[RVPoly([0,1],1)]);", "zero vector"),
        ("rvOrdAt[3,[1]](s);", "two field coordinates"),
        ("rvOrdAt[7,[1,4]](s);", "coefficient 4 is out of range"),
        (
            "rvOrdAt[3,[1,0]](rvmon[[0,1,2]]{1,0,0});",
            "requires a binary form",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}\n{statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
    for (statement, reason) in [
        (
            format!("rvFrobOrders[C,{}];", usize::BITS),
            "Frobenius extension exponent",
        ),
        (
            format!(
                "rvcurve D[[0,1]] {{rvMul(s,s),rvMul(t,t)}}; rvFrobOrders[D,{}];",
                usize::BITS - 1
            ),
            "Frobenius row degree is out of range",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}\n{statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
}
