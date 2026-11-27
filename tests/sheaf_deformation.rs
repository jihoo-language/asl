mod common;

use common::{diagnostic, example, interpret, stdout};

const PARALLEL: &str = include_str!("../examples/sheaf-deformation.asl");

#[test]
fn the_parallel_arms_first_obstruct_in_order_five() {
    let result = stdout(example("sheaf-deformation"));
    assert_eq!(result, "ShExt(F, F, 2){g0}\n");
    let wrong_coefficient = diagnostic(interpret(&format!(
        "{PARALLEL} show shDefAppend(j1,zero[ShHomC(F,F,1)]);"
    )));
    assert!(wrong_coefficient.contains("does not satisfy dA + curvature = 0"));
    let blocked = diagnostic(interpret(&format!("{PARALLEL} show shDefLift(j4);")));
    assert!(
        blocked.contains("nonzero obstruction in sheaf Ext degree two"),
        "{blocked}"
    );
}

#[test]
fn valid_jets_remain_functorial_under_noncommuting_matrix_gauges() {
    let source = r"
        space C {a<b<c};
        sheaf F[C,5] {a:2,b:2,c:2} {a<b:[1,1,0,1],b<c:[1,0,2,1]};
        let j0 = shDefStart[F];
        let z = shDefLift(shDefLift(shDefLift(j0)));
        let qa = shHomTerm[F,F,[a],1,1];
        let qb = shHomTerm[F,F,[b],1,0];
        let qc = add(shHomTerm[F,F,[c],0,0],shHomTerm[F,F,[c],0,1]);
        let q = add(qa,add(qb,qc));
        let k = shDefGauge[1](z,q);
        let c1 = shD(q);
        let c2 = shScale[-1](shExtCompose(c1,q));
        let c3 = shScale[-1](shExtCompose(c2,q));
        assert shDefTerm[1](k) == c1;
        assert shDefTerm[2](k) == c2;
        assert shDefTerm[3](k) == c3;
        assert shDefTerm[0](k) == shDefTerm[0](z);
        let rab = add(shHomTerm[F,F,[a,b],0,0],
            add(shHomTerm[F,F,[a,b],1,0],shHomTerm[F,F,[a,b],1,1]));
        let rbc = add(shHomTerm[F,F,[b,c],0,0],
            add(shHomTerm[F,F,[b,c],0,1,2],shHomTerm[F,F,[b,c],1,1]));
        let rac = add(shHomTerm[F,F,[a,c],0,0],
            add(shHomTerm[F,F,[a,c],1,0],
            add(shHomTerm[F,F,[a,c],0,1,2],shHomTerm[F,F,[a,c],1,1,3])));
        assert shDefTerm[0](z) == add(rab,add(rbc,rac));
        let w = shDefGaugeWitness[1](z,q);
        let curvatureDifference = add(shDefCurvature(z),shScale[-1](shDefCurvature(k)));
        assert shD(w) == curvatureDifference;
        assert shDefObstruction(k) == shDefObstruction(z);
        let k1 = shDefAppend(j0,shDefTerm[1](k));
        let k2 = shDefAppend(k1,shDefTerm[2](k));
        assert shDefAppend(k2,shDefTerm[3](k)) == k;
        assert shD(shDefTerm[4](shDefLift(k))) == shScale[-1](shDefCurvature(k));
        let l = shDefGauge[2](k,qc);
        assert shDefTerm[1](l) == shDefTerm[1](k);
        assert shDefTerm[2](l) == add(shDefTerm[2](k),shD(qc));
        assert shD(shDefGaugeWitness[2](k,qc)) ==
            add(shDefCurvature(k),shScale[-1](shDefCurvature(l)));
        assert shDefCurvature(k) == shHomTerm[F,F,[a,b,c],1,1];
        show shDefCurvature(k);
    ";
    let got = stdout(interpret(source));
    assert!(got.starts_with("ShHomC(F, F, 2){"), "{got}");
    assert_ne!(got, "ShHomC(F, F, 2){0}\n");
}

#[test]
fn an_obstructed_nonzero_class_is_preserved_under_a_change_of_basis() {
    let source = format!(
        r"{PARALLEL}
        let q = shHomTerm[F,F,[x0],0,0,2];
        let k4 = shDefGauge[1](j4,q);
        let w = shDefGaugeWitness[1](j4,q);
        assert shD(w) == add(shDefCurvature(j4),shScale[-1](shDefCurvature(k4)));
        assert shDefObstruction(k4) == shDefObstruction(j4);
        assert shDefTerm[1](k4) == a1;
        "
    );
    assert_ne!(stdout(interpret(&source)), "ShExt(F, F, 2){0}\n");
}

#[test]
fn jets_are_first_class_finite_values_with_checked_domains_and_coefficients() {
    let source = r"
        space B {a<b};
        sheaf F[B,3] {a:1,b:1} {a<b:[0]};
        card (\j : ShDef(F,0) => S.hi);
        card (\j : ShDef(F,2) => S.hi);
        let one = shHomTerm[F,F,[a,b],0,0];
        let j1 = shDefAppend(shDefStart[F],one);
        let j2 = shDefAppend(j1,shScale[-1](one));
        let coeff = (\j : ShDef(F,2) => shDefTerm[2](j));
        assert coeff(j2) == shScale[-1](one);
        assert (\j : ShDef(F,2) => j)(j2) == j2;
        show j2;
    ";
    assert_eq!(
        stdout(interpret(source)),
        "card = 1 / 1\ncard = 9 / 9\nShDef(F, 2){t^1=ShHomC(F, F, 1){[a<b]#0->0}, t^2=ShHomC(F, F, 1){2*[a<b]#0->0}}\n"
    );
}

#[test]
fn invalid_prefixes_and_types_are_rejected() {
    let source = r"
        space C {a<b<c};
        sheaf E[C,5] {a:1,b:1,c:1} {a<b:[0],b<c:[0]};
        sheaf V[C,5] {a:1,b:1,c:1} {a<b:[1],b<c:[1]};
        let j0 = shDefStart[E];
        let j1 = shDefLift(j0);
        let b = shHomTerm[E,E,[a,c],0,0];
        let h = shHomTerm[E,E,[a],0,0];
    ";
    for (statement, expected) in [
        (
            "show shDefAppend(j0,b);",
            "does not satisfy dA + curvature = 0",
        ),
        ("show shDefAppend(j0,h);", "ShHomC(E, E, 1)"),
        (
            "show shDefAppend(j0,shHomTerm[V,V,[a,b],0,0]);",
            "ShHomC(E, E, 1)",
        ),
        ("show shDefTerm[2](j1);", "beyond the known truncation"),
        (
            "show shDefGauge[0](j1,h);",
            "between 1 and the truncation order",
        ),
        (
            "show shDefGauge[2](j1,h);",
            "between 1 and the truncation order",
        ),
        ("show shDefGauge[1](j1,b);", "ShHomC(E, E, 0)"),
        ("show shDefCurvature(h);", "a ShDef jet"),
        ("show shDefStart[Bad];", "unknown sheaf"),
        (
            r"card (\j : ShDef(V,3) => S.hi);",
            "sheaf deformations exceed the finite type enumeration limit",
        ),
        (
            r"show (\j : ShDef(E,18446744073709551615) => shDefLift(j));",
            "order is out of range",
        ),
    ] {
        let blocked = diagnostic(interpret(&format!("{source}{statement}")));
        assert!(blocked.contains(expected), "{statement}: {blocked}");
    }
}
