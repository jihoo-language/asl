#[allow(dead_code)]
mod common;

use common::{diagnostic, interpret, stdout};

const QUADRICS: &str = "
    rvideal Q[[0,1]] {rvmon[[0,1]]{2,0}, rvmon[[0,1]]{0,2}};
    let A=rvmkterm[Q,[0]](rvmon[[0,1]]{1,0});
    let B=rvmkterm[Q,[1]](rvmon[[0,1]]{0,1});
    let a=rvTorClass(A);
    let b=rvTorClass(B);
";

#[test]
fn local_products_non_squarefree_classes_and_forget_are_compatible() {
    let source = format!(
        "{QUADRICS}
        rvmbetti[Q,1,[2,0]];
        rvmbetti[Q,1,[3,0]];
        rvmbetti[Q,2,[2,2]];
        assert a == rvmtorgen[Q,1,[2,0],0];
        assert rvkD(A) == zero[RVMKoszul(Q,0,[2,0])];
        assert rvTorRep(a) == A;
        let m = rvTorProd(a,b);
        assert m == rvmtorgen[Q,2,[2,2],0];
        assert rvMForget(m) == rvTorProd(rvMForget(a),rvMForget(b));
        assert rvMForget(rvkProd(A,B)) == rvkProd(rvMForget(A),rvMForget(B));
        assert rvMComponent[[2,2]](rvMForget(m)) == m;
        assert rvMComponent[[1,3]](rvMForget(m)) == zero[RVMTor(Q,2,[1,3])];
        show m;
        "
    );
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "rvmbetti[Q, 1, [2, 0]] = 1\n",
            "rvmbetti[Q, 1, [3, 0]] = 0\n",
            "rvmbetti[Q, 2, [2, 2]] = 1\n",
            "RVMTor(Q, 2, [2, 2]){g0}\n"
        )
    );
}

#[test]
fn raw_components_split_total_chains_and_keep_differential_and_preimages() {
    let source = format!(
        r"{QUADRICS}
        let edge=rvmkterm[Q,[0,1]](rvmon[[0,1]]{{0,0}});
        let bd=rvkD(edge);
        let Z=add(rvMForget(A),rvMForget(bd));
        assert rvMComponent[[2,0]](Z) == A;
        assert rvMComponent[[1,1]](Z) == bd;
        assert rvMComponent[[0,2]](Z) == zero[RVMKoszul(Q,1,[0,2])];
        assert rvMForget(bd) == rvkD(rvMForget(edge));
        assert rvkPrimitive(bd) == edge;
        assert rvMForget(rvkPrimitive(bd)) == rvkPrimitive(rvMForget(bd));
        let sum=add(rvMForget(a),rvMForget(b));
        assert rvMComponent[[2,0]](sum)==a;
        assert rvMComponent[[0,2]](sum)==b;
        let identity=\x:RVMKoszul(Q,1,[1,1]) => x;
        assert identity(bd) == bd;
        card \x:RVMKoszul(Q,1,[1,1]) => S.hi;
        "
    );
    assert_eq!(stdout(interpret(&source)), "card = 4 / 4\n");
}

#[test]
fn one_local_fiber_does_not_require_a_large_total_basis() {
    let n = usize::MAX;
    let source = format!(
        r"
        rvideal Unit[[0,1,2]] {{rvmon[[0,1,2]]{{0,0,0}}}};
        rvideal Empty[[]] {{}};
        rvideal J[[0,1,2]] {{rvmon[[0,1,2]]{{0,1,1}}}};
        rvmbetti[Unit,0,[0,0,0]];
        rvmbetti[Empty,0,[]];
        rvmbetti[J,0,[{n},0,0]];
        let m=rvmkterm[J,[]](rvmon[[0,1,2]]{{{n},0,0}});
        assert rvkD(rvkPrimitive(m))==m;
        card \x:RVMKoszul(J,0,[{n},0,0]) => S.hi;
        "
    );
    assert_eq!(
        stdout(interpret(&source)),
        format!(
            "rvmbetti[Unit, 0, [0, 0, 0]] = 0\nrvmbetti[Empty, 0, []] = 1\nrvmbetti[J, 0, [{n}, 0, 0]] = 0\ncard = 2 / 2\n"
        )
    );
}

#[test]
fn local_grade_and_conversion_errors_are_clear() {
    for (input, expected) in [
        ("rvmbetti[Q,0,[1]];", "one exponent"),
        (
            "rvideal Binomial[[0,1]] {add(rvmon[[0,1]]{1,0},rvmon[[0,1]]{0,1})}; rvmbetti[Binomial,0,[0,0]];",
            "requires a monomial",
        ),
        ("show rvMComponent[[0,0]](rvMForget(a));", "same total internal degree"),
        ("show rvMComponent[[2,0]](a);", "requires a total"),
        ("show rvMForget(rvMForget(a));", "requires a multigraded"),
        ("show rvTorProd(a,rvMForget(b));", "two total grades or two multidegrees"),
        (
            "let c=rvmon[[0,1]]{1,0}; show rvmkterm[Q,[0]](c);",
            "requires an rvmon literal",
        ),
        ("show rvmkterm[Q,[0]](rvmon[[0,2]]{1,0});", "different primitive selection"),
        ("show rvmkterm[Q,[0]](rvmon[[0,1]]{1});", "one exponent"),
    ] {
        let error = diagnostic(interpret(&format!("{QUADRICS}\n{input}")));
        assert!(error.contains(expected), "{input}: {error}");
    }
    let bad = format!("{QUADRICS} rvmbetti[Q,0,[{},1]];", usize::MAX);
    assert!(diagnostic(interpret(&bad)).contains("out of range"));
}
