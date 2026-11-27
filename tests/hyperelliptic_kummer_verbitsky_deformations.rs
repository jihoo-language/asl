mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const MARK: &str = "rvhodgewick B[21]{0,1}{[2,0],[0,-2]};";

#[test]
fn graded_hom_quotients_match_independent_pfaffian_and_complete_intersection() {
    let output = stdout(example("hyperelliptic-kummer-verbitsky-deformations"));
    let lines = output.lines().collect::<Vec<_>>();
    assert_eq!(lines[0], "rvVerbitskyDeformationTable[B, 2] = [(shift, tangent, T1, Ext1), (-5, 0, 0, 4), (-4, 0, 0, 8), (-3, 0, 0, 12), (-2, 6, 5, 6), (-1, 12, 8, 0), (0, 8, 2, 0), (1, 4, 0, 0)] totals = (30, 15, 30)");
    assert_eq!(lines[1], "rvVerbitskyDeformationDims[B, 2, -2] = (6, 5, 6)");
    assert!(lines[3].contains("residual = [(0, {12*D[0]*D[1]}), (1, {-126*D[0]*D[1]^2})]"));
    assert!(lines[3].contains("Ext1 class = []; lift verified; h = [(0, {-12*D[0]})]"));
    assert!(lines[3].contains("C = [(1, [(1, {126})])"));
    // For a regular sequence of weights 3 and 4, the Hom differential
    // vanishes over R. Ext¹(I,R)_w is R_{7+w}, even though every tangent lifts.
    assert_eq!(lines[4], "rvVerbitskyDeformationTable[C, 2] = [(shift, tangent, T1, Ext1), (-7, 0, 0, 1), (-6, 0, 0, 1), (-5, 0, 0, 2), (-4, 1, 1, 1), (-3, 2, 2, 1), (-2, 3, 2, 0), (-1, 3, 1, 0), (0, 2, 0, 0), (1, 1, 0, 0)] totals = (12, 6, 6)");
    assert_eq!(lines[5], "rvVerbitskyObstruction = residual = []; Ext1 class = []; lift verified; h = []; B = [(0, [(1, {-1})])]; C = []");
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Q(i128, i128);

impl Q {
    const ZERO: Self = Self(0, 1);
    fn new(n: i128, d: i128) -> Self {
        let (mut a, mut b) = (n.abs(), d);
        while b != 0 {
            (a, b) = (b, a % b);
        }
        Self(n / a, d / a)
    }
    fn parse(text: &str) -> Self {
        text.split_once('/').map_or_else(
            || Self(text.parse().unwrap(), 1),
            |(a, b)| Self(a.parse().unwrap(), b.parse().unwrap()),
        )
    }
    fn add(self, rhs: Self) -> Self {
        Self::new(self.0 * rhs.1 + rhs.0 * self.1, self.1 * rhs.1)
    }
    fn times(self, rhs: Self) -> Self {
        Self::new(self.0 * rhs.0, self.1 * rhs.1)
    }
}

fn polynomial(text: &str, vars: [i128; 3]) -> Q {
    let text = text.strip_prefix('{').unwrap().strip_suffix('}').unwrap();
    text.split(", ")
        .filter(|term| !term.is_empty())
        .fold(Q::ZERO, |sum, term| {
            let mut factors = term.split('*');
            let mut value = Q::parse(factors.next().unwrap());
            for factor in factors {
                let (variable, exponent) = factor
                    .split_once('^')
                    .map_or((factor, 1), |(v, e)| (v, e.parse().unwrap()));
                let index = match variable {
                    "D[0]" => 0,
                    "D[1]" => 1,
                    "Q[0]" => 2,
                    _ => panic!("unexpected variable {variable}"),
                };
                value = value.times(Q(vars[index].pow(exponent), 1));
            }
            sum.add(value)
        })
}

fn sparse<T: Copy, const N: usize>(
    text: &str,
    absent: T,
    terminator: &str,
    parse: impl Fn(&str) -> T,
) -> [T; N] {
    let mut text = text.strip_prefix('[').unwrap().strip_suffix(']').unwrap();
    let mut values = [absent; N];
    while !text.is_empty() {
        let (index, value) = text.strip_prefix('(').unwrap().split_once(", ").unwrap();
        let end = value.find(terminator).unwrap() + terminator.len();
        values[index.parse::<usize>().unwrap()] = parse(&value[..end - 1]);
        text = value[end..].strip_prefix(", ").unwrap_or(&value[end..]);
    }
    values
}

fn dot(left: &[Q; 5], column: &[Q; 5]) -> Q {
    left.iter()
        .zip(column)
        .fold(Q::ZERO, |sum, (a, b)| sum.add(a.times(*b)))
}

#[test]
fn an_independent_rational_evaluation_checks_both_polynomial_lift_equations() {
    let mut source = String::from(MARK);
    for i in 0..5 {
        writeln!(source, "rvVerbitskyResolutionEntry[B,2,1,0,{i}];").unwrap();
    }
    for j in 0..5 {
        for i in 0..5 {
            writeln!(source, "rvVerbitskyResolutionEntry[B,2,2,{i},{j}];").unwrap();
        }
    }
    source.push_str(
        "let g=add(rvVerbitskyTangentBasis[B,2,-1,0],rvVerbitskyTangentBasis[B,2,-1,3]);
         show g;
         rvVerbitskyObstruction(g);",
    );
    let out = stdout(interpret(&source));
    let lines = out.lines().collect::<Vec<_>>();
    let entries = |index: usize| lines[index].split_once(" = ").unwrap().1;
    let report = lines[31];
    let h = report.split_once("; h = ").unwrap().1;
    let (h, b) = h.split_once("; B = ").unwrap();
    let (b, c) = b.split_once("; C = ").unwrap();
    for vars in [[1, -1, 1], [2, 1, -3], [3, 5, 4]] {
        let vector =
            |text: &str| -> [Q; 5] { sparse(text, Q::ZERO, "})", |value| polynomial(value, vars)) };
        let matrix = |text: &str| -> [[Q; 5]; 5] { sparse(text, [Q::ZERO; 5], "])", vector) };
        let f = std::array::from_fn(|i| polynomial(entries(i), vars));
        let a: [[Q; 5]; 5] = std::array::from_fn(|j| {
            std::array::from_fn(|i| polynomial(entries(5 + 5 * j + i), vars))
        });
        let coordinate_start = lines[30].find(")[").unwrap() + 1;
        let g = vector(&lines[30][coordinate_start..]);
        let h = vector(h);
        let b = matrix(b);
        let c = matrix(c);
        assert!(h.iter().any(|&x| x != Q::ZERO));
        assert!(c.iter().flatten().any(|&x| x != Q::ZERO));
        for j in 0..5 {
            assert_eq!(dot(&f, &a[j]), Q::ZERO);
            assert_eq!(dot(&f, &b[j]).add(dot(&g, &a[j])), Q::ZERO);
            assert_eq!(
                dot(&f, &c[j]).add(dot(&g, &b[j])).add(dot(&h, &a[j])),
                Q::ZERO
            );
        }
    }
}

#[test]
fn nominal_addition_canonical_quotients_and_actual_tangent_coordinates() {
    let output = stdout(interpret(&format!(
        "{MARK}
        let unit=rvVerbitskyBasis[B,2,8,0];
        let g=rvVerbitskyTangent[B,2,1]{{(0,unit)}};
        let x=rvVerbitskyTangent[B,2,1]{{(0,unit),(0,unit)}};
        assert x==add(g,g);
        assert add(x,neg(x))==zero[RVVerbitskyTangent(B,2,1)];
        assert rvVerbitskyDeformationScale[1,2](x)==g;
        assert rvVerbitskyDeformationClass(g)==zero[RVVerbitskyTOne(B,2,1)];
        let t=rvVerbitskyTOneBasis[B,2,-2,4];
        assert rvVerbitskyDeformationClass(rvVerbitskyDeformationRepresentative(t))==t;
        let y=rvVerbitskyTangentBasis[B,2,-2,0];
        assert rvVerbitskyDeformationClass(add(y,y))==rvVerbitskyDeformationScale[2,1](rvVerbitskyDeformationClass(y));
        rvVerbitskyDeformationDims[B,2,-1000000000000];
        rvVerbitskyObstruction(zero[RVVerbitskyTangent(B,2,1000000000000)]);
        "
    )));
    assert_eq!(output.lines().count(), 2);
    assert!(output.starts_with("rvVerbitskyDeformationDims[B, 2, -1000000000000] = (0, 0, 0)"));
    assert!(output.ends_with("h = []; B = []; C = []\n"));
}

#[test]
fn rejects_incompatible_coordinates_nontangents_and_unbounded_resolution() {
    for (query, reason) in [
        ("let g=rvVerbitskyTangent[B,2,-3]{(0,rvVerbitskyBasis[B,2,0,0])};", "fail the syzygy equations"),
        ("let g=rvVerbitskyTangent[B,2,-3]{(0,rvVerbitskyBasis[B,2,2,0])};", "RVVerbitsky(B,2,0) coordinate"),
        ("let g=rvVerbitskyTangent[B,2,-9]{(0,rvVerbitskyBasis[B,2,0,0])};", "no ring component"),
        ("let g=rvVerbitskyTOneBasis[B,2,-2,5];", "basis index out of range"),
        ("rvVerbitskyObstruction(rvVerbitskyTOneBasis[B,2,-2,0]);", "RVVerbitskyTangent vector"),
        ("let x=rvVerbitskyDeformationScale[1,0](rvVerbitskyTangentBasis[B,2,-2,0]);", "denominator must be positive"),
        ("rvVerbitskyDeformationTable[B,1000000000000000000];", "storage limit"),
        ("rvVerbitskyDeformationDims[B,1,0];", "requires total Beauville rank 22"),
        ("rvhodgewick C[21]{0,1}{[2,0],[0,-2]}; let z=add(zero[RVVerbitskyTangent(B,2,-2)],zero[RVVerbitskyTangent(C,2,-2)]);", "RVVerbitskyTangent"),
    ] {
        assert!(diagnostic(interpret(&format!("{MARK} {query}"))).contains(reason), "{query}");
    }
}
