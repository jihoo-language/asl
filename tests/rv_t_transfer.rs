#[allow(dead_code)]
mod common;

use std::fmt::Write;

use common::{diagnostic, interpret, stdout};

#[test]
fn unary_lift_binary_product_and_morphism_identity_hold_over_a_nonmonomial_ideal() {
    let source = format!(
        "{}
        let left_q = rvTorClass(u);
        let right_q = rvTorClass(v);
        assert rvTLift(left_q) == rvTorRep(left_q);
        assert rvTHigher(left_q,right_q) == rvTorProd(left_q,right_q);
        assert rvkD(rvTLift(left_q,right_q)) ==
            add(rvkProd(rvTLift(left_q),rvTLift(right_q)),rvTorRep(rvTHigher(left_q,right_q)));
        ",
        include_str!("../examples/rank-syzygies.asl")
    );
    assert!(stdout(interpret(&source)).ends_with("RVTor(DoublePoint, 2, 4){0}\n"));
}

#[test]
fn a_noncycle_binary_sum_uses_the_full_homotopy_and_nonzero_morphism_corrections() {
    let source = format!(
        "{}
        let raw3 = add(rvkProd(rvTLift(a),rvTLift(c,b)),rvkProd(rvTLift(a,c),rvTLift(b)));
        let nonzeroD = rvkterm[Flag,[1,3,5]](rvmon[[1,2,3,4,5,6]]{{0,1,0,1,0,1}});
        assert rvkD(raw3) == nonzeroD;
        let chosen = rvTHigher(a,c,b);
        assert chosen == rvtorgen[Flag,4,6,0];
        let correction = add(rvTLift(rvTHigher(a,c),b),rvTLift(a,rvTHigher(c,b)));
        assert rvkD(correction) == nonzeroD;
        assert rvkD(rvTLift(a,c,b)) == add(add(raw3,rvTorRep(chosen)),correction);
        assert rvTHigher(rvTHigher(a,c),b) == rvTHigher(a,rvTHigher(c,b));
        let unit = rvtorgen[Flag,0,0,0];
        assert rvTHigher(unit,c) == c;
        assert rvTHigher(unit,b,c) == zero[RVTor(Flag,3,4)];
        assert rvTLift(unit,b,c) == zero[RVKoszul(Flag,4,4)];
        ",
        include_str!("../examples/rank-massey.asl")
    );
    assert!(stdout(interpret(&source)).ends_with("RVTor(Flag, 4, 6){g0}\n"));
}

#[test]
fn the_arity_four_stasheff_terms_cancel_with_a_nonzero_summand() {
    let pairs = [(1, 2), (3, 4), (5, 6), (1, 3), (4, 5), (7, 7)];
    let mut program = ideal(7, &pairs);
    for (name, i, j) in [('a', 1, 2), ('b', 3, 4), ('c', 5, 6), ('d', 7, 7)] {
        writeln!(program, "let {name} = rvTorClass({});", term(7, &[j], &[i])).unwrap();
    }
    program.push_str(
        "
        let u=rvTHigher(rvTHigher(a,c,b),d);
        let v=rvTHigher(a,rvTHigher(c,b,d));
        let w=rvTHigher(rvTHigher(a,c),b,d);
        let y=rvTHigher(a,rvTHigher(c,b),d);
        let z=rvTHigher(a,c,rvTHigher(b,d));
        assert u == rvtorgen[J,5,8,0];
        assert add(add(u,v),add(w,add(y,z))) == zero[RVTor(J,5,8)];
        ",
    );
    assert!(stdout(interpret(&program)).is_empty());
}

#[test]
fn staircase_operations_first_appear_at_arities_three_through_six() {
    for n in 3..=6 {
        let pairs = (1..=n)
            .flat_map(|j| {
                (1..=j)
                    .filter(move |&k| (j, k) != (n, 1))
                    .map(move |k| (j, n + k))
            })
            .collect::<Vec<_>>();
        let mut program = ideal(2 * n, &pairs);
        for j in 1..=n {
            writeln!(
                program,
                "let a{j} = rvTorClass({});",
                term(2 * n, &[n + j], &[j])
            )
            .unwrap();
        }
        let inputs = |start, end| {
            (start..=end)
                .map(|j| format!("a{j}"))
                .collect::<Vec<_>>()
                .join(",")
        };
        for span in 2..n {
            for start in 1..=n - span + 1 {
                let args = inputs(start, start + span - 1);
                writeln!(
                    program,
                    "assert rvTHigher({args}) == zero[RVTor(J,{},{})];",
                    2 * span - 2,
                    2 * span
                )
                .unwrap();
            }
        }
        let wedge = (1..=2 * n)
            .filter(|&k| k != n && k != n + 1)
            .collect::<Vec<_>>();
        let args = inputs(1, n);
        writeln!(
            program,
            "let target = rvTorClass({});",
            term(2 * n, &wedge, &[n, n + 1])
        )
        .unwrap();
        writeln!(
            program,
            "assert target == rvtorgen[J,{},{},0];",
            2 * n - 2,
            2 * n
        )
        .unwrap();
        writeln!(program, "assert rvTHigher({args}) == target;").unwrap();

        if n == 4 {
            for j in 1..=n {
                let local = term(2 * n, &[n + j], &[j]).replacen("rvkterm", "rvmkterm", 1);
                writeln!(program, "let l{j} = rvTorClass({local});").unwrap();
            }
            program.push_str(
                "
                assert rvMForget(rvTHigher(l1,l2,l3,l4)) == target;
                assert rvMForget(rvTLift(l1,l2,l3,l4)) == rvTLift(a1,a2,a3,a4);
                let u = rvkProd(rvTLift(a1),rvTLift(a2,a3,a4));
                let v = rvkProd(rvTLift(a1,a2),rvTLift(a3,a4));
                let w = rvkProd(rvTLift(a1,a2,a3),rvTLift(a4));
                assert rvkD(rvTLift(a1,a2,a3,a4)) == add(add(u,v),add(w,rvTorRep(target)));
                ",
            );
        }
        assert!(stdout(interpret(&program)).is_empty(), "arity {n}");
    }
}

#[test]
fn transfer_checks_arity_input_family_and_degree_overflow() {
    let program = "rvideal J[[0]] {}; rvideal K[[0]] {}; let a = rvtorgen[J,0,0,0];";
    for (expr, message) in [
        ("rvTHigher(a)", "at least two"),
        ("rvTLift()", "at least one"),
        ("rvTLift(One.pt)", "requires RVTor"),
        ("rvTHigher(a,zero[RVTor(K,0,0)])", "same rank ideal"),
    ] {
        let error = diagnostic(interpret(&format!("{program} show {expr};")));
        assert!(error.contains(message), "{expr}: {error}");
    }
    let overflow = format!("{program} show rvTLift(a,zero[RVTor(J,{},0)]);", usize::MAX);
    assert!(diagnostic(interpret(&overflow)).contains("out of range"));
}

fn monomial(width: usize, factors: &[usize]) -> String {
    let selected: Vec<_> = (1..=width).collect();
    let mut powers = vec![0; width];
    for &index in factors {
        powers[index - 1] += 1;
    }
    let powers = powers
        .iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join(",");
    format!("rvmon[{selected:?}]{{{powers}}}")
}

fn term(width: usize, wedge: &[usize], powers: &[usize]) -> String {
    format!("rvkterm[J,{wedge:?}]({})", monomial(width, powers))
}

fn ideal(width: usize, pairs: &[(usize, usize)]) -> String {
    let selected: Vec<_> = (1..=width).collect();
    let generators = pairs
        .iter()
        .map(|&(i, j)| monomial(width, &[i, j]))
        .collect::<Vec<_>>()
        .join(",");
    format!("rvideal J[{selected:?}] {{{generators}}};\n")
}
