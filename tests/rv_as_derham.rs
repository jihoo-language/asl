mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn nonmonomial_cech_pairing_and_binary_field_semilinearity() {
    assert_eq!(
        stdout(example("rank-as-derham")),
        concat!(
            "ASDeRham(C, 11){4*eta2}\n",
            "ASDeRham(C, 11){6*alpha0 + 6*alpha2}\n",
            "BinaryScalar(11){1}\n",
            "BinaryScalar(11){4}\n",
            "BinaryScalar(11){2}\n",
            "rvASEO[C] = [0, 1, 1]\n",
        )
    );
}

fn coordinates(line: &str, genus: usize) -> usize {
    let (_, support) = line.strip_suffix('}').unwrap().split_once('{').unwrap();
    if support == "0" {
        return 0;
    }
    support
        .split(" + ")
        .map(|term| {
            if let Some(index) = term.strip_prefix("alpha") {
                1 << index.parse::<usize>().unwrap()
            } else {
                1 << (genus + term.strip_prefix("eta").unwrap().parse::<usize>().unwrap() - 1)
            }
        })
        .fold(0, |sum, term| sum | term)
}

fn act(columns: &[usize], vector: usize) -> usize {
    columns
        .iter()
        .enumerate()
        .filter(|(i, _)| vector >> i & 1 == 1)
        .fold(0, |image, (_, column)| image ^ column)
}

fn rank(columns: &[usize]) -> usize {
    let mut rows = vec![0; columns.len()];
    for &column in columns {
        let mut image = column;
        while image != 0 {
            let pivot = image.ilog2() as usize;
            if rows[pivot] == 0 {
                rows[pivot] = image;
                break;
            }
            image ^= rows[pivot];
        }
    }
    rows.into_iter().filter(|&row| row != 0).count()
}

#[test]
fn frobenius_verschiebung_pairing_and_canonical_flag_use_the_full_module() {
    for g in 1..=6 {
        let m = 2 * g + 1;
        let mut numerator = format!("rvmon[[0,1]]{{0,{m}}}");
        if g >= 2 {
            numerator = format!("add({numerator},rvmon[[0,1]]{{{},3}})", m - 3);
        }
        let mut input = format!("rvascover C[[0,1]] {{ {numerator},rvmon[[0,1]]{{{m},0}} }};\n");
        for i in 0..2 * g {
            writeln!(input, "let x{i} = rvASDeRhamGen[C,3,{i}];").unwrap();
            writeln!(input, "show rvASF(x{i}); show rvASV(x{i});").unwrap();
        }
        for i in 0..2 * g {
            for j in 0..2 * g {
                writeln!(input, "show rvASDeRhamPair(x{i},x{j});").unwrap();
            }
        }
        writeln!(input, "rvASEO[C];").unwrap();
        let output = stdout(interpret(&input));
        let mut lines = output.lines();
        let mut f = Vec::new();
        let mut v = Vec::new();
        for _ in 0..2 * g {
            f.push(coordinates(lines.next().unwrap(), g));
            v.push(coordinates(lines.next().unwrap(), g));
        }
        let mut j = Vec::new();
        for _ in 0..2 * g {
            let mut column = 0;
            for i in 0..2 * g {
                let entry = lines.next().unwrap();
                if entry == "BinaryScalar(3){1}" {
                    column |= 1 << i;
                } else {
                    assert_eq!(entry, "BinaryScalar(3){0}");
                }
            }
            j.push(column);
        }
        // Rank g for both maps and FV=VF=0 certify both image/kernel equalities.
        assert_eq!(rank(&f), g);
        assert_eq!(rank(&v), g);
        assert_eq!(rank(&j), 2 * g);
        for i in 0..2 * g {
            assert_eq!(act(&f, v[i]), 0);
            assert_eq!(act(&v, f[i]), 0);
            assert_eq!(j[i] >> i & 1, 0);
            for k in 0..2 * g {
                assert_eq!(j[i] >> k & 1, j[k] >> i & 1);
                assert_eq!(
                    act(&j, f[i]) >> k & 1,
                    (j[i] & v[k]).count_ones() as usize & 1
                );
            }
        }
        assert_eq!(
            lines.next().unwrap(),
            format!(
                "rvASEO[C] = {:?}",
                (1..=g).map(|i| i / 2).collect::<Vec<_>>()
            )
        );
    }
}

#[test]
fn zero_genus_and_sparse_large_genus_keep_their_proper_bounds() {
    let program = r"
        rvascover Line[[0,1]] {rvmon[[0,1]]{0,1},rvmon[[0,1]]{1,0}};
        rvASEO[Line];
        let x = zero[ASDeRham(Line,11)];
        assert rvASF(x) == x;
        assert rvASV(x) == x;
        assert rvASDeRhamPair(x,x) == zero[BinaryScalar(11)];
        card (\x : ASDeRham(Line,11) => S.hi);
        card (\y : BinaryScalar(11) => S.hi);
    ";
    assert_eq!(
        stdout(interpret(program)),
        "rvASEO[Line] = []\ncard = 1 / 1\ncard = 8 / 8\n"
    );
    let g = 1usize << (usize::BITS / 2 - 1);
    let degree = 2 * g + 1;
    let setup = format!(
        "rvascover Large[[0,1]] {{rvmon[[0,1]]{{0,{degree}}},rvmon[[0,1]]{{{degree},0}}}};"
    );
    assert_eq!(
        stdout(interpret(&format!(
            "{setup} show rvASF(rvASDeRhamGen[Large,11,{}]);",
            2 * g - 1
        ))),
        "ASDeRham(Large, 11){alpha0}\n"
    );
    let error = diagnostic(interpret(&format!("{setup} rvASEO[Large];")));
    assert!(error.contains("storage limit"), "{error}");
}

#[test]
fn rejects_invalid_coordinates_fields_and_nonpolynomial_charts() {
    let setup = r"
        rvascover C[[0,1]] {rvmon[[0,1]]{0,5},rvmon[[0,1]]{5,0}};
        rvascover Other[[0,1]] {rvmon[[0,1]]{0,5},rvmon[[0,1]]{5,0}};
        rvascover Ordinary[[0,1]] {add(rvmon[[0,1]]{2,0},rvmon[[0,1]]{0,2}),rvmon[[0,1]]{1,1}};
    ";
    for (statement, reason) in [
        ("show rvASDeRhamGen[C,3,4];", "out of range for dimension 4"),
        ("show zero[ASDeRham(C,5)];", "reducible over F2"),
        (
            "show rvFieldElt[11,8];",
            "field coefficient 8 is out of range",
        ),
        (
            "show rvASDeRhamScale[8](zero[ASDeRham(C,11)]);",
            "field coefficient 8 is out of range",
        ),
        ("show rvASV(rvASOmega[C,3,0]);", "an ASDeRham class"),
        (
            "show rvASDeRhamIncl(zero[ASDeRham(C,3)]);",
            "an ASDiff differential",
        ),
        (
            "show rvASDeRhamPair(zero[ASDeRham(C,3)],zero[ASDeRham(Other,3)]);",
            "expected type ASDeRham(C, 3)",
        ),
        (
            "show rvASDeRhamPair(zero[ASDeRham(C,3)],zero[ASDeRham(C,11)]);",
            "expected type ASDeRham(C, 3)",
        ),
        (
            "show zero[ASDeRham(Ordinary,11)];",
            "requires a one-pole polynomial cover",
        ),
        (
            "show zero[ASDeRham(Missing,11)];",
            "unknown Artin–Schreier cover",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}\n{statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
}
