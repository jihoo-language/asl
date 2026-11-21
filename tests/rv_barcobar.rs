mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn bar_and_split_terms_give_the_exact_differential_and_leibniz_rule() {
    assert_eq!(
        stdout(example("rank-barcobar")),
        concat!(
            "RVBarCobar(P, 0, 3, 4){",
            "[(0, 1, 1):g0] * [(1, 2, 1):g0 | (2, 3, 1):g0] + ",
            "[(0, 1, 1):g0 | (1, 2, 1):g0] * [(2, 3, 1):g0] + ",
            "[(0, 1, 1):g0 | (1, 3, 2):g0] + ",
            "[(0, 2, 2):g0 | (2, 3, 1):g0]}\n"
        )
    );
}

#[test]
fn the_retract_identity_holds_on_every_chain_in_three_finite_slices() {
    let program = include_str!("../examples/rank-barcobar.asl");
    for (b, size) in [(3, 16), (4, 16), (5, 2)] {
        let source = format!(
            r"{program}
            let identity = \x : RVBarCobar(P,0,3,{b}) => x;
            let retract = \x : RVBarCobar(P,0,3,{b}) =>
                add(add(rvBCD(rvBCHomotopy(x)),rvBCHomotopy(rvBCD(x))),
                    rvBCDesuspend(rvBarWord[P,[0,3]](rvBCProject(x))));
            assert retract == identity;
            let checkd2 = \x : RVBarCobar(P,0,3,{b}) => rvBCD(rvBCD(x));
            assert checkd2 == (\x : RVBarCobar(P,0,3,{b}) => zero[RVBarCobar(P,0,3,{out})]);
            card \x : RVBarCobar(P,0,3,{b}) => S.hi;
            ",
            out = b - 2,
        );
        let output = stdout(interpret(&source));
        assert!(
            output.ends_with(&format!("card = {size} / {size}\n")),
            "{output}"
        );
    }
}

#[test]
fn the_first_staircase_split_is_nonclosed_and_projects_to_the_nonzero_top() {
    let source = format!(
        r"{}
        rvbarpath P[Stair] {{
            [1,0,0,0,1,0,0,0], [0,1,0,0,0,1,0,0],
            [0,0,1,0,0,0,1,0], [0,0,0,1,0,0,0,1]
        }};
        let raw = rvBarWord[P,[0,1,2,3,4]](
            rvMComponent[[1,0,0,0,1,0,0,0]](a),
            rvMComponent[[0,1,0,0,0,1,0,0]](b),
            rvMComponent[[0,0,1,0,0,0,1,0]](c),
            rvMComponent[[0,0,0,1,0,0,0,1]](d));
        let ADJ = rvBCDesuspend(rvBarWord[P,[0,1]](rvMComponent[[1,0,0,0,1,0,0,0]](a)));
        let BCD = rvBCDesuspend(rvBarWord[P,[1,2,3,4]](
            rvMComponent[[0,1,0,0,0,1,0,0]](b),
            rvMComponent[[0,0,1,0,0,0,1,0]](c),
            rvMComponent[[0,0,0,1,0,0,0,1]](d)));
        let w = rvBCMul(ADJ,BCD);
        let second = rvBCDesuspend(rvBarWord[P,[1,2]](rvMComponent[[0,1,0,0,0,1,0,0]](b)));
        let third = rvBCDesuspend(rvBarWord[P,[3,4]](rvMComponent[[0,0,0,1,0,0,0,1]](d)));
        let tail = rvBCDesuspend(rvBarWord[P,[2,3,4]](
            rvMComponent[[0,0,1,0,0,0,1,0]](c),
            rvMComponent[[0,0,0,1,0,0,0,1]](d)));
        let head = rvBCDesuspend(rvBarWord[P,[1,2,3]](
            rvMComponent[[0,1,0,0,0,1,0,0]](b),
            rvMComponent[[0,0,1,0,0,0,1,0]](c)));
        let left_gap = rvBCMul(second,tail);
        let right_gap = rvBCMul(head,third);
        assert rvBCD(w) == rvBCMul(ADJ,add(left_gap,right_gap));
        let long = rvBCDesuspend(raw);
        let top = rvBarWord[P,[0,4]](rvMComponent[[1,1,1,1,1,1,1,1]](value));
        let incl = rvBCDesuspend(top);
        assert rvBCProject(w) == rvmtorgen[Stair,6,[1,1,1,1,1,1,1,1],0];
        assert rvBCHomotopy(w) == long;
        assert add(rvBCD(rvBCHomotopy(w)),rvBCHomotopy(rvBCD(w))) == add(w,incl);
        assert rvBCD(rvBCD(w)) == zero[RVBarCobar(P,0,4,4)];
        assert rvBCProject(rvBCD(w)) == zero[RVMTor(Stair,5,[1,1,1,1,1,1,1,1])];
        assert rvBCHomotopy(rvBCHomotopy(w)) == zero[RVBarCobar(P,0,4,8)];
        assert rvBCHomotopy(incl) == zero[RVBarCobar(P,0,4,7)];
        let retra = \x : RVBarCobar(P,0,4,6) =>
            add(add(rvBCD(rvBCHomotopy(x)),rvBCHomotopy(rvBCD(x))),
                rvBCDesuspend(rvBarWord[P,[0,4]](rvBCProject(x))));
        assert retra == (\x : RVBarCobar(P,0,4,6) => x);
        ",
        include_str!("../examples/rank-massey-four.asl")
    );
    assert!(stdout(interpret(&source)).ends_with("RVTor(Stair, 6, 8){g0}\n"));
}

#[test]
fn rectification_operations_check_paths_degrees_and_input_kind() {
    let src = include_str!("../examples/rank-barcobar.asl");
    for (expr, expected) in [
        ("rvBCMul(A,A)", "sharing an endpoint"),
        ("rvBCMul(A,a)", "requires an RVBarCobar"),
        ("rvBCProject(a)", "requires an RVBarCobar"),
        ("rvBCDesuspend(A)", "requires an RVBar chain"),
        ("rvBCD(zero[RVBarCobar(P,0,1,0)])", "positive"),
    ] {
        let error = diagnostic(interpret(&format!("{src} show {expr};")));
        assert!(error.contains(expected), "{expr}: {error}");
    }
    let error = diagnostic(interpret(&format!(
        "{src} show rvBCHomotopy(zero[RVBarCobar(P,0,1,{})]);",
        usize::MAX
    )));
    assert!(error.contains("degree is out of range"), "{error}");
}
