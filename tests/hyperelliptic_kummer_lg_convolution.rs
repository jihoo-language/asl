mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const PREFIX: &str = "rvverbitskytate T{5}; rvlgdefect I[T]{0,1,2}; rvlgdefect J[T]{1,2,3,4}; rvlgdefectconvolution F[I,J];";

#[test]
fn example_checks_shifted_free_blocks_and_an_actual_noninterval_factorization() {
    let output = stdout(example("hyperelliptic-kummer-lg-convolution"));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(lines.len(), 4, "{output}");
    assert!(lines[0].starts_with(
        "rvLGConvolution[F] = W=x^5−z^5 supports=[[0, 1, 2], [1, 2, 3, 4]] rank=12 weights="
    ));
    assert!(lines[0].ends_with("power=verified factorization=verified"));
    assert!(lines[1].contains("(shift,support,free)=[(0, [0, 1, 2, 3, 4], true), (1, [0, 1, 2, 3, 4], true), (2, [3, 4], false)]"));
    assert!(lines[1].ends_with("intertwining=verified determinant=constant-unit inverse=verified"));
    assert!(lines[2].starts_with("rvLGDefect[K] = W=x^5−y^5 support=[3, 4]"));
    assert!(lines[3].starts_with("rvLGConvolution[Q] = W=x^4−z^4 supports=[[0, 2], [1, 3]] rank=4"));
}

#[test]
fn interval_cyclic_splitting_matches_independent_lengths_shifts_and_wraparound() {
    for d in 2usize..=6 {
        let mut source = format!("rvverbitskytate T{{{d}}};");
        let mut expected = Vec::new();
        for r in 1..d {
            for (start, name) in [(d - 1, "L"), (1, "R")] {
                let roots = (0..r)
                    .map(|i| ((start + i) % d).to_string())
                    .collect::<Vec<_>>()
                    .join(",");
                writeln!(source, "rvlgdefect {name}{r}[T]{{{roots}}};").unwrap();
            }
        }
        for r in 1..d {
            for s in 1..d {
                writeln!(
                    source,
                    "rvlgdefectconvolution F{r}_{s}[L{r},R{s}]; rvLGConvolutionSplit[F{r}_{s}];"
                )
                .unwrap();
                let free = (r + s).saturating_sub(d);
                let parts: Vec<_> = (0..r.min(s))
                    .map(|g| {
                        if g < free {
                            (g, (0..d).collect::<Vec<_>>(), true)
                        } else {
                            let l = r + s - 1 - 2 * g;
                            (g, (0..l).map(|i| (g + i) % d).collect(), false)
                        }
                    })
                    .collect();
                expected.push(format!("(shift,support,free)={parts:?}"));
            }
        }
        let output = stdout(interpret(&source));
        assert_eq!(output.lines().count(), expected.len());
        for (actual, expected) in output.lines().zip(expected) {
            assert!(
                actual.contains(&expected),
                "d={d}; {actual} should contain {expected}"
            );
        }
    }
}

#[test]
fn graph_fusion_exposes_all_exact_polynomial_coefficients_and_the_inverse() {
    for d in 2usize..=8 {
        let a = d - 1;
        let b = d - 1;
        let mut source=format!("rvverbitskytate T{{{d}}}; rvlgdefect D[T]{{{a}}}; rvlgdefect E[T]{{{b}}}; rvlgdefectconvolution F[D,E]; rvLGConvolutionSplit[F];
             let one=rvCycloRat[1](rvCyclo[{d},0]); let nothing=zero[RVCycloQ({d})];");
        let power = (a + b) % d;
        writeln!(source, "let z=rvCycloRat[1](rvCyclo[{d},{power}]);").unwrap();
        source.push_str(
            "assert rvLGConvolutionTerm[F,0,0,0,0,1]==z;
             assert rvLGConvolutionTerm[F,1,0,0,0,1]==neg(z);
             assert rvLGConvolutionTerm[F,1,0,0,1,0]==one;
             assert rvLGConvolutionTerm[F,3,0,0,0,0]==one;
             assert rvLGConvolutionTerm[F,4,0,0,0,0]==one;
             assert rvLGConvolutionTerm[F,5,0,0,0,1]==z;",
        );
        for k in 0..d {
            let x = d - 1 - k;
            writeln!(
                source,
                "assert rvLGConvolutionTerm[F,2,0,0,{x},{k}]==rvCycloRat[1](rvCyclo[{d},{}]);",
                power * k
            )
            .unwrap();
            writeln!(
                source,
                "assert rvLGConvolutionTerm[F,2,0,0,{x},{}]==nothing;",
                k + 1
            )
            .unwrap();
        }
        let report = stdout(interpret(&source));
        assert!(report.contains(&format!("(shift,support,free)=[(0, [{power}], false)]")));
    }
}

#[test]
fn arbitrary_supports_have_two_actual_polynomial_factorization_identities() {
    let mut source = String::from(
        "rvverbitskytate T{4}; rvlgdefect I[T]{0,2}; rvlgdefect J[T]{1,3};
         rvlgdefectconvolution F[I,J];
         let one=rvCycloRat[1](rvCyclo[4,0]); let nothing=zero[RVCycloQ(4)];",
    );
    let weights = [0usize, 1, 1, 2];
    for (r, &wr) in weights.iter().enumerate() {
        for (c, &wc) in weights.iter().enumerate() {
            let total = 4 + wc - wr;
            for x in 0..=total {
                let z = total - x;
                let expected = if r == c && x == 4 {
                    "one"
                } else if r == c && x == 0 {
                    "neg(one)"
                } else {
                    "nothing"
                };
                for (first, second) in [(1, 2), (2, 1)] {
                    let mut products = Vec::new();
                    for t in 0..4 {
                        for i in 0..=x {
                            for j in 0..=z {
                                products.push(format!("rvCycloQMul(rvLGConvolutionTerm[F,{first},{r},{t},{i},{j}],rvLGConvolutionTerm[F,{second},{t},{c},{},{}])",x-i,z-j));
                            }
                        }
                    }
                    // Keep both the companion and its polynomial factors in the public language.
                    while products.len() > 1 {
                        products = products
                            .chunks(2)
                            .map(|p| match p {
                                [a, b] => format!("add({a},{b})"),
                                [a] => a.clone(),
                                _ => unreachable!(),
                            })
                            .collect();
                    }
                    writeln!(source, "assert {}=={expected};", products[0]).unwrap();
                }
            }
        }
    }
    assert!(stdout(interpret(&source)).is_empty());
}

#[test]
fn convolution_diagnoses_wrong_marks_degrees_and_unavailable_certificates() {
    for (query, reason) in [
        ("rvlgdefectconvolution F[J,I];", "already defined"),
        (
            "rvverbitskytate U{5}; rvlgdefect K[U]{1}; rvlgdefectconvolution R[I,K];",
            "same Tate marking",
        ),
        ("rvLGConvolution[Z];", "unknown LG convolution"),
        (
            "show rvLGConvolutionTerm[F,6,0,0,0,0];",
            "matrix kind must be between zero and five",
        ),
        (
            "show rvLGConvolutionTerm[F,0,12,0,0,0];",
            "matrix index is out of range",
        ),
        (
            "rvlgconvolutionpart R[F,0];",
            "free block with full root support",
        ),
        (
            "rvlgconvolutionpart R[F,3];",
            "summand index is out of range",
        ),
        (
            "rvlgdefect K[T]{0,2}; rvlgdefectconvolution W[I,K]; rvLGConvolutionSplit[W];",
            "split requires consecutive circular intervals",
        ),
        (
            "rvlgdefect K[T]{0,2}; rvlgdefectconvolution W[I,K]; show rvLGConvolutionTerm[W,3,0,0,1,0];",
            "split requires consecutive circular intervals",
        ),
        (
            "rvverbitskytate U{9}; rvlgdefect K[U]{0,1,2,3,4,5,6,7}; rvlgdefectconvolution W[K,K];",
            "coefficient storage limit",
        ),
    ] {
        let actual = diagnostic(interpret(&format!("{PREFIX} {query}")));
        assert!(actual.contains(reason), "{query}: {actual}");
    }
}
