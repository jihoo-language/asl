mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn setup() -> &'static str {
    include_str!("../examples/hyperelliptic-og6-spin-monodromy.asl")
        .split_once("// Left regular")
        .unwrap()
        .0
}

fn spin_data(name: &str, rank: usize, factor: &str) -> String {
    let factors = vec![factor; rank].join(", ");
    let middle = 64 - 2 * rank;
    format!("rvOG6SpinData[{name}] = rank=64, abelian_dimension=32, log_rank={rank}, graded=[{rank}, {middle}, {rank}], smith=[{factors}]\n")
}

#[test]
fn example_reports_type_ii_and_type_iii_weight_one_and_smith() {
    assert_eq!(
        stdout(example("hyperelliptic-og6-spin-monodromy")),
        spin_data("A", 16, "3") + &spin_data("B3", 32, "1")
    );
}

#[test]
fn all_polarized_anticommutators_and_pbw_order_are_checked_on_a_shear() {
    let mut source = format!(
        r"{}
        let h2=add(e,add(h,h));
        rvog6cusp Q(h2,add(e,f),f);
        let unit=rvOG6CliffordBasis[Q,0];
        ",
        setup()
    );
    for i in 0..7 {
        writeln!(
            source,
            "let x{i}=rvOG6PolarizedBasis[Q,{i}]; let a{i}=rvOG6CliffordVector[Q](x{i}); assert a{i}==rvOG6CliffordBasis[Q,{}];",
            1 << i
        ).unwrap();
    }
    let mut pairings = String::new();
    for i in 0..7 {
        for j in i..7 {
            writeln!(pairings, "rvOG6Pair(x{i},x{j});").unwrap();
        }
    }
    let actual = stdout(interpret(&format!("{source}{pairings}")));
    let mut entries = actual.lines().map(|line| {
        line.strip_prefix("rvOG6Pair = ")
            .unwrap()
            .parse::<i128>()
            .unwrap()
    });
    let mut nonorthogonal = false;
    for i in 0..7 {
        for j in i..7 {
            let gram = entries.next().unwrap();
            writeln!(
                source,
                "assert add(rvOG6CliffordMul(a{i},a{j}),rvOG6CliffordMul(a{j},a{i}))==rvOG6CliffordScale[{gram}](unit);"
            ).unwrap();
            if i < j {
                nonorthogonal |= gram != 0;
                writeln!(
                    source,
                    "assert rvOG6CliffordMul(a{i},a{j})==rvOG6CliffordBasis[Q,{}];",
                    (1 << i) | (1 << j)
                )
                .unwrap();
            }
        }
    }
    assert!(nonorthogonal);
    assert_eq!(entries.next(), None);
    source.push_str(
        r"
        let b=rvOG6CliffordBasis[Q,127];
        let c=rvOG6CliffordBasis[Q,29];
        assert rvOG6CliffordMul(rvOG6CliffordMul(b,c),a3)
            ==rvOG6CliffordMul(b,rvOG6CliffordMul(c,a3));
        ",
    );
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn spin_lifts_square_to_zero_and_conjugate_the_whole_polarized_basis() {
    let mut source = format!(
        r"{}
        let h2=add(e,add(h,h));
        rvog6cusp Q(h2,add(e,f),f);
        rvog6nilpotentii Q2[Q,-3];
        let qeDual=add(add(eDual,eDual),neg(hF));
        rvog6nilpotentiii Q3[Q,1](add(e,qeDual));
        rvog6nilpotentiii Shift[P,0](add(u,e));
        assert rvOG6SpinLog[Shift]==rvOG6SpinLog[B3];
        ",
        setup()
    );
    for (name, cusp) in [("A", "P"), ("B3", "P"), ("Q2", "Q"), ("Q3", "Q")] {
        writeln!(
            source,
            r"
            let log{name}=rvOG6SpinLog[{name}];
            let exp{name}=rvOG6SpinExp[{name},-3];
            let inv{name}=rvOG6SpinExp[{name},3];
            let unit{name}=rvOG6CliffordBasis[{cusp},0];
            assert rvOG6CliffordMul(log{name},log{name})==zero[RVOG6Clifford({cusp},0)];
            assert rvOG6CliffordMul(exp{name},inv{name})==unit{name};
            assert rvOG6CliffordMul(rvOG6SpinExp[{name},7],exp{name})==rvOG6SpinExp[{name},4];
            assert exp{name}==add(unit{name},rvOG6CliffordScale[-3](log{name}));
        "
        )
        .unwrap();
        for i in 0..7 {
            writeln!(
                source,
                r"
                let v{name}{i}=rvOG6PolarizedBasis[{cusp},{i}];
                let c{name}{i}=rvOG6CliffordVector[{cusp}](v{name}{i});
                let lhs{name}{i}=rvOG6CliffordMul(exp{name},c{name}{i});
                assert rvOG6CliffordMul(lhs{name}{i},inv{name})
                    ==rvOG6CliffordVector[{cusp}](rvOG6NilpotentExp[{name},-3](v{name}{i}));
            "
            )
            .unwrap();
        }
    }
    assert_eq!(stdout(interpret(&source)), "");
}

#[test]
fn rank_from_exported_spin_images_agrees_with_exact_smith_reports() {
    let mut source = setup().to_owned();
    for name in ["A", "B3"] {
        for mask in (0u32..128).filter(|mask| mask.count_ones().is_multiple_of(2)) {
            writeln!(
                source,
                "show rvOG6CliffordMul(rvOG6SpinLog[{name}],rvOG6CliffordBasis[P,{mask}]);"
            )
            .unwrap();
        }
    }
    let output = stdout(interpret(&source));
    let columns: Vec<_> = output.lines().map(parse_sparse_column).collect();
    assert_eq!(columns.len(), 128);
    assert_eq!(rank_mod_101(columns[..64].to_vec()), 16);
    assert_eq!(rank_mod_101(columns[64..].to_vec()), 32);
}

fn parse_sparse_column(line: &str) -> [i32; 128] {
    let terms = line
        .strip_prefix("RVOG6Clifford(P,0){[")
        .unwrap()
        .strip_suffix("] / 1}")
        .unwrap();
    let mut column = [0; 128];
    if !terms.is_empty() {
        for term in terms.split("), (") {
            let (mask, coefficient) = term.trim_matches(['(', ')']).split_once(", ").unwrap();
            column[mask.parse::<usize>().unwrap()] = coefficient.parse().unwrap();
        }
    }
    column
}

fn rank_mod_101(mut columns: Vec<[i32; 128]>) -> usize {
    let mut rank = 0;
    for row in 0..128 {
        let Some(pivot) = (rank..columns.len()).find(|&i| columns[i][row].rem_euclid(101) != 0)
        else {
            continue;
        };
        columns.swap(rank, pivot);
        let pivot = columns[rank];
        for column in &mut columns[rank + 1..] {
            let c = column[row];
            for (entry, &p) in column.iter_mut().zip(&pivot) {
                *entry = (*entry * pivot[row] - c * p).rem_euclid(101);
            }
        }
        rank += 1;
    }
    rank
}

#[test]
fn scaled_and_sheared_spin_smith_use_the_actual_integral_marking() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{}
            rvog6nilpotentii Huge[P,-170141183460469231731687303715884105728];
            rvog6nilpotentiii Stretch[P,0](add(add(u,u),e));
            let h2=add(e,add(h,h));
            rvog6cusp Q(h2,e,f);
            rvog6nilpotentii C2[Q,1];
            rvog6nilpotentiii C3[Q,0](u);
            rvOG6SpinData[Huge];
            rvOG6SpinData[Stretch];
            rvOG6SpinData[C2];
            rvOG6SpinData[C3];
            ",
            setup()
        ))),
        spin_data("Huge", 16, "170141183460469231731687303715884105728")
            + &spin_data("Stretch", 32, "2")
            + &spin_data("C2", 16, "1")
            + &spin_data("C3", 32, "1")
    );
}

#[test]
fn rejects_unpolarized_vectors_bad_masks_and_different_nominal_marks() {
    for (script, error) in [
        ("show rvOG6PolarizedBasis[P,7];", "between 0 and 6"),
        ("show rvOG6CliffordBasis[P,128];", "between 0 and 127"),
        ("show zero[RVOG6Clifford(P,2)];", "parity must be 0 or 1"),
        (
            "show rvOG6CliffordVector[P](h);",
            "perpendicular to the polarization",
        ),
        ("show rvOG6SpinLog[Missing];", "unknown OG6 nilpotent"),
        (
            "rvOG6CliffordCoeff[128](rvOG6SpinLog[A]);",
            "between 0 and 127",
        ),
        ("rvOG6CliffordCoeff[0](h);", "an RVOG6Clifford element"),
        (
            "rvog6cusp Q(h,e,f); show rvOG6CliffordMul(rvOG6SpinLog[A],rvOG6CliffordBasis[Q,0]);",
            "expected type RVOG6Clifford(P,0)",
        ),
        (
            "show add(rvOG6CliffordBasis[P,0],rvOG6CliffordBasis[P,1]);",
            "expected type RVOG6Clifford(P,0)",
        ),
    ] {
        let actual = diagnostic(interpret(&format!("{}{script}", setup())));
        assert!(actual.contains(error), "{actual}");
    }
}
