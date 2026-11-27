mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn resolution_detects_the_same_weight_syzygy_and_its_zero_quartic_entry() {
    let report = stdout(example("hyperelliptic-kummer-verbitsky-resolution"));
    let lines = report.lines().collect::<Vec<_>>();
    assert_eq!(
        &lines[..4],
        [
            "rvVerbitskyResolution[B, 2] = ranks=[1, 5, 5, 1] shifts=[[0], [3, 3, 3, 3, 4], [4, 5, 5, 5, 5], [8]]",
            "rvVerbitskyResolutionEntry[B, 2, 1, 0, 0] = {-2/7*D[0]*Q[0], 1*D[0]^3}",
            "rvVerbitskyResolutionEntry[B, 2, 2, 0, 0] = {1*D[1]}",
            "rvVerbitskyResolutionEntry[B, 2, 2, 4, 0] = {}",
        ]
    );
    assert!(lines[4].starts_with("rvVerbitskyPfaffian[B, 2] = scale="));
}

#[derive(Clone, Copy, Debug)]
struct Q(i128, i128);

impl Q {
    fn parse(text: &str) -> Self {
        text.split_once('/').map_or_else(
            || Self(text.parse().unwrap(), 1),
            |(a, b)| Self(a.parse().unwrap(), b.parse().unwrap()),
        )
    }
    fn add(self, rhs: Self) -> Self {
        Self(self.0 * rhs.1 + rhs.0 * self.1, self.1 * rhs.1)
    }
    fn multiply(self, rhs: Self) -> Self {
        Self(self.0 * rhs.0, self.1 * rhs.1)
    }
    fn negate(self) -> Self {
        Self(-self.0, self.1)
    }
    fn equal(self, rhs: Self) -> bool {
        self.0 * rhs.1 == rhs.0 * self.1
    }
}

fn eval(polynomial: &str, xyz: [i128; 3]) -> Q {
    if polynomial.is_empty() {
        return Q(0, 1);
    }
    polynomial.split(", ").fold(Q(0, 1), |sum, term| {
        let mut elements = term.split('*');
        let mut value = Q::parse(elements.next().unwrap());
        for factor in elements {
            let (base, exponent) = factor
                .split_once('^')
                .map_or((factor, 1), |(base, exp)| (base, exp.parse().unwrap()));
            let coordinate = match base {
                "D[0]" => xyz[0],
                "D[1]" => xyz[1],
                "Q[0]" => xyz[2],
                _ => panic!("unknown variable {base}"),
            };
            value = value.multiply(Q(coordinate.pow(exponent), 1));
        }
        sum.add(value)
    })
}

#[test]
fn an_independent_evaluation_sees_every_signed_pfaffian_with_a_common_scale() {
    let output = stdout(interpret(
        "rvhodgewick B[21]{0,1}{[2,0],[0,-2]}; rvVerbitskyPfaffian[B,2];",
    ));
    let certificate = output.split_once(" = scale=").unwrap().1;
    let (scale, entries) = certificate.split_once(" matrix=").unwrap();
    let scale = Q::parse(scale);
    assert_ne!(scale.0, 0);
    let entries = entries
        .split('}')
        .filter_map(|field| field.rsplit_once('{').map(|(_, p)| p))
        .collect::<Vec<_>>();
    assert_eq!(entries.len(), 25);
    let xyz = [2, -1, 3];
    let a = |i: usize, j: usize| eval(entries[5 * i + j], xyz);
    let original = [
        "-2/7*D[0]*Q[0], 1*D[0]^3",
        "1*D[0]*D[1]^2, 2/21*D[0]*Q[0]",
        "1*D[0]^2*D[1], -2/21*D[1]*Q[0]",
        "2/7*D[1]*Q[0], 1*D[1]^3",
        "1*D[1]^2*Q[0], 2/23*Q[0]^2",
    ];
    for (omitted, generator) in original.iter().enumerate() {
        let indices = (0..5).filter(|&i| i != omitted).collect::<Vec<_>>();
        let [i, j, k, l] = <[usize; 4]>::try_from(indices).unwrap();
        let pf = a(i, j)
            .multiply(a(k, l))
            .add(a(i, k).multiply(a(j, l)).negate())
            .add(a(i, l).multiply(a(j, k)));
        let signed = if omitted.is_multiple_of(2) {
            pf
        } else {
            pf.negate()
        };
        assert!(
            signed.equal(scale.multiply(eval(generator, xyz))),
            "signed minor {omitted}"
        );
        if omitted == 4 {
            assert_ne!(signed.0, 0);
        }
    }
}

#[test]
fn complete_intersection_differentials_and_augmentation_are_product_compatible() {
    let out = stdout(interpret(
        r"
        rvhodgewick C[22]{0,1}{[2]};
        rvVerbitskyResolution[C,2];
        let xz=rvVerbitskyResolutionBasis[C,2,0,0]{1,1};
        assert rvVerbitskyResolutionAugment(xz)==rvVerbitskyCup(rvVerbitskyBasis[C,2,2,0],rvVerbitskyClass(rvHodgeQ[C,2]{1}));
        let f=rvVerbitskyResolutionBasis[C,2,1,0];
        let zf=rvVerbitskyFreeTimes{0,1}(f);
        assert zf==rvVerbitskyResolutionBasis[C,2,1,0]{0,1};
        assert rvVerbitskyResolutionDiff(zf)==rvVerbitskyFreeTimes{0,1}(rvVerbitskyResolutionDiff(f));
        assert rvVerbitskyResolutionAugment(rvVerbitskyResolutionDiff(f))==zero[RVVerbitsky(C,2,6)];
        assert rvVerbitskyResolutionAugment(rvVerbitskyResolutionDiff(rvVerbitskyResolutionBasis[C,2,1,1]))==zero[RVVerbitsky(C,2,8)];
        let syz=rvVerbitskyResolutionBasis[C,2,2,0];
        assert rvVerbitskyResolutionDiff(rvVerbitskyResolutionDiff(syz))==zero[RVVerbitskyFree(C,2,0,7)];
        assert add(syz,neg(syz))==zero[RVVerbitskyFree(C,2,2,7)];
        ",
    ));
    assert_eq!(
        out,
        "rvVerbitskyResolution[C, 2] = ranks=[1, 2, 1] shifts=[[0], [3, 4], [7]]\n"
    );
}

#[test]
fn accepts_symbolic_top_overflow_in_free_polynomials_and_checks_grades_and_indices() {
    let text = stdout(interpret(
        r"
      rvhodgewick C[22]{0,1}{[2]};
      let high=rvVerbitskyResolutionBasis[C,2,0,0]{10};
      assert rvVerbitskyResolutionAugment(high)==zero[RVVerbitsky(C,2,20)];
      assert rvVerbitskyResolutionDiff(high)==zero[RVVerbitskyFree(C,2,-1,10)];
      ",
    ));
    assert!(text.is_empty());
    let mark = "rvhodgewick C[22]{0,1}{[2]}; ";
    for (query, reason) in [
        (
            "rvVerbitskyResolutionEntry[C,2,1,1,0];",
            "outside the minimal resolution",
        ),
        ("rvVerbitskyPfaffian[C,2];", "requires codimension three"),
        (
            "let x=rvVerbitskyResolutionAugment(rvVerbitskyResolutionBasis[C,2,1,0]);",
            "homological degree zero",
        ),
        (
            "let x=rvVerbitskyResolutionBasis[C,2,0,0]{0,0,0};",
            "too many variable exponents",
        ),
        (
            "rvVerbitskyResolution[C,1000000000000000000];",
            "storage limit",
        ),
    ] {
        assert!(
            diagnostic(interpret(&format!("{mark}{query}"))).contains(reason),
            "{query}"
        );
    }
}

#[test]
fn further_weight_patterns_admit_exact_alternating_presentations() {
    let text = stdout(interpret(
        r"
        rvhodgewick B[21]{0,1}{[2,0],[0,-2]};
        rvVerbitskyResolution[B,3];
        rvVerbitskyPfaffian[B,3];
        rvhodgewick E[11]{-1,-1,1}{[2]};
        rvVerbitskyResolution[E,2];
        rvVerbitskyPfaffian[E,2];
        ",
    ));
    let rows = text.lines().collect::<Vec<_>>();
    assert_eq!(rows.len(), 4);
    assert_eq!(rows[0], "rvVerbitskyResolution[B, 3] = ranks=[1, 5, 5, 1] shifts=[[0], [4, 4, 4, 4, 4], [6, 6, 6, 6, 6], [10]]");
    assert!(rows[1].starts_with("rvVerbitskyPfaffian[B, 3] = scale="));
    assert_eq!(rows[2], "rvVerbitskyResolution[E, 2] = ranks=[1, 5, 5, 1] shifts=[[0], [3, 3, 4, 4, 4], [5, 5, 5, 6, 6], [9]]");
    assert!(rows[3].starts_with("rvVerbitskyPfaffian[E, 2] = scale="));
}
