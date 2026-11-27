mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

fn marking() -> &'static str {
    include_str!("../examples/hyperelliptic-mukai-moduli.asl")
        .split_once("// Hilbert cube")
        .unwrap()
        .0
}

#[test]
fn example_reports_hilbert_cube_and_a_nonspecial_primitive_vector() {
    assert_eq!(
        stdout(example("hyperelliptic-mukai-moduli")),
        concat!(
            "rvMukaiModuliLattice[Hilb] = n=3, dimension=6, rank=23, signature=(3,20), determinant=4, discriminant=Z/4, q=-1/4 mod 2\n",
            "rvMukaiModuliDiv[Hilb] = 4\n",
            "rvMukaiModuliFujiki[Hilb, 0] = -960\n",
            "rvMukaiModuliRR[Hilb] = 4\n",
            "rvMukaiModuliRR[Hilb] = -20\n",
            "rvMukaiModuliLattice[Mod] = n=27, dimension=54, rank=23, signature=(3,20), determinant=52, discriminant=Z/52, q=-1/52 mod 2\n",
            "rvMukaiPair = 0\n",
            "rvMukaiPair = 0\n",
            "rvMukaiModuliDiv[Mod] = 52\n",
        )
    );
}

#[test]
fn every_kernel_basis_vector_is_perpendicular_and_gram_determinant_proves_saturation() {
    let mut source = format!(
        r"{}
        let h=rvK3FromKummerDiv(rvKummerHyperplane[B]);
        let v=rvMukai[2,-4](add(add(h,h),h));
        rvmukaimoduli M(v);
        let g=rvMukaiModuliGlue[M];
        rvMukaiPair(g,v);
        rvMukaiPair(g,g);
        ",
        marking(),
    );
    for i in 0..23 {
        writeln!(
            source,
            "let x{i}=rvMukaiModuliBasis[M,{i}]; rvMukaiPair(v,x{i});"
        )
        .unwrap();
    }
    for i in 0..23 {
        for j in i..23 {
            writeln!(source, "rvMukaiPair(x{i},x{j});").unwrap();
        }
    }
    let output = stdout(interpret(&source));
    let mut values = output.lines().map(|line| {
        line.strip_prefix("rvMukaiPair = ")
            .unwrap()
            .parse::<i128>()
            .unwrap()
    });
    assert_eq!(values.next(), Some(0));
    let glue_square = values.next().unwrap();
    assert_eq!((glue_square + 52).rem_euclid(2 * 52 * 52), 0);
    for _ in 0..23 {
        assert_eq!(values.next(), Some(0));
    }
    let mut gram = [[0i128; 23]; 23];
    for i in 0..23 {
        for j in i..23 {
            let value = values.next().unwrap();
            gram[i][j] = value;
            gram[j][i] = value;
        }
    }
    assert_eq!(values.next(), None);
    assert_eq!(determinant(gram), 52);
}

fn determinant<const N: usize>(mut matrix: [[i128; N]; N]) -> i128 {
    let (mut sign, mut last) = (1, 1);
    for i in 0..(N - 1) {
        let pivot = (i..N).find(|&j| matrix[j][i] != 0).unwrap();
        if pivot != i {
            matrix.swap(pivot, i);
            sign = -sign;
        }
        let diagonal = matrix[i][i];
        for j in i + 1..N {
            for k in i + 1..N {
                let numerator = diagonal * matrix[j][k] - matrix[j][i] * matrix[i][k];
                assert_eq!(numerator % last, 0);
                matrix[j][k] = numerator / last;
            }
            matrix[j][i] = 0;
        }
        last = diagonal;
    }
    sign * matrix[N - 1][N - 1]
}

#[test]
fn ideal_sheaf_formulas_agree_with_existing_api_and_riemann_roch_covers_both_signs() {
    let mut source = marking().to_owned();
    for n in 2..=4 {
        let s = n - 1;
        writeln!(
            source,
            r"
            rvmukaimoduli M{n}(rvMukaiIdeal[B,{n}]);
            let d{n}=rvMukai[1,{s}](zero[RVK3Lattice(B)]);
            let h{n}=rvMukai[0,0](rvK3FromKummerDiv(rvKummerHyperplane[B]));
            let x{n}=add(h{n},d{n});
            "
        )
        .unwrap();
        for j in 0..=2 * n {
            writeln!(
                source,
                "rvMukaiHilbFujiki[{n},{j}](x{n},d{n}); rvMukaiModuliFujiki[M{n},{j}](x{n},d{n});"
            )
            .unwrap();
        }
    }
    let output = stdout(interpret(&source));
    let values: Vec<_> = output
        .lines()
        .map(|line| line.split_once(" = ").unwrap().1)
        .collect();
    assert_eq!(values.len(), 42);
    assert!(values.chunks_exact(2).all(|pair| pair[0] == pair[1]));

    let output = stdout(interpret(&format!(
        r"{}
        let h=rvMukai[0,0](rvK3FromKummerDiv(rvKummerHyperplane[B]));
        let d2=rvMukai[1,1](zero[RVK3Lattice(B)]);
        let d3=rvMukai[1,2](zero[RVK3Lattice(B)]);
        rvmukaimoduli M2(rvMukaiIdeal[B,2]);
        rvmukaimoduli M3(rvMukaiIdeal[B,3]);
        rvMukaiModuliRR[M2](add(add(d2,d2),d2));
        rvMukaiModuliRR[M3](h);
        rvMukaiModuliRR[M3](d3);
        rvMukaiModuliRR[M3](add(d3,d3));
        rvMukaiModuliDiv[M3](zero[RVMukai(B)]);
        ",
        marking()
    )));
    assert_eq!(
        output,
        concat!(
            "rvMukaiModuliRR[M2] = 21\n",
            "rvMukaiModuliRR[M3] = 20\n",
            "rvMukaiModuliRR[M3] = 0\n",
            "rvMukaiModuliRR[M3] = -20\n",
            "rvMukaiModuliDiv[M3] = 0\n",
        )
    );
}

#[test]
fn line_tensor_and_spherical_reflection_preserve_perpendicular_arithmetic() {
    let output = stdout(interpret(&format!(
        r"{}
        let h=rvKummerHyperplane[B];
        let v=rvMukaiIdeal[B,3];
        let x=rvMukai[1,2](zero[RVK3Lattice(B)]);
        let y=rvMukai[0,0](rvK3FromKummerDiv(h));
        let root=rvMukaiLine(h);
        rvmukaimoduli I(v);
        rvmukaimoduli T(rvMukaiTensor(v,h));
        rvmukaimoduli R(rvMukaiTwist(root,v));
        rvMukaiModuliDiv[I](x);
        rvMukaiModuliDiv[T](rvMukaiTensor(x,h));
        rvMukaiModuliDiv[R](rvMukaiTwist(root,x));
        rvMukaiModuliFujiki[I,2](x,y);
        rvMukaiModuliFujiki[T,2](rvMukaiTensor(x,h),rvMukaiTensor(y,h));
        rvMukaiModuliFujiki[R,2](rvMukaiTwist(root,x),rvMukaiTwist(root,y));
        rvMukaiModuliRR[I](y);
        rvMukaiModuliRR[T](rvMukaiTensor(y,h));
        rvMukaiModuliRR[R](rvMukaiTwist(root,y));
        ",
        marking()
    )));
    let values: Vec<_> = output
        .lines()
        .map(|line| line.split_once(" = ").unwrap().1)
        .collect();
    assert_eq!(
        values,
        ["4", "4", "4", "192", "192", "192", "20", "20", "20"]
    );
}

#[test]
fn dimension_can_exceed_signed_range_and_trivial_rr_needs_one_factor() {
    let n = 1i128 << 126;
    let output = stdout(interpret(&format!(
        "{} rvmukaimoduli M(rvMukaiIdeal[B,{n}]); rvMukaiModuliLattice[M]; rvMukaiModuliRR[M](zero[RVMukai(B)]);",
        marking(),
    )));
    assert!(
        output.contains(&format!("dimension={}", 2 * n as u128)),
        "{output}"
    );
    assert!(
        output.contains(&format!("rvMukaiModuliRR[M] = {}", n + 1)),
        "{output}"
    );
}

#[test]
fn errors_check_primitivity_square_orthogonality_names_grades_and_resource_bounds() {
    let setup = format!("{} let v=rvMukaiIdeal[B,3]; rvmukaimoduli M(v);", marking());
    for (source, expected) in [
        ("rvmukaimoduli Bad(add(v,v));", "requires a primitive Mukai vector"),
        ("rvmukaimoduli Bad(zero[RVMukai(B)]);", "requires a primitive Mukai vector"),
        ("rvmukaimoduli Bad(rvMukaiIdeal[B,0]);", "requires a Mukai vector of positive even square"),
        ("rvmukaimoduli Bad(rvK3LatticeBasis[B,0]);", "an RVMukai vector"),
        ("rvmukaimoduli M(v);", "already defined"),
        ("rvMukaiModuliLattice[Missing];", "unknown Mukai moduli lattice"),
        ("show rvMukaiModuliBasis[M,23];", "between 0 and 22"),
        ("rvMukaiModuliDiv[M](v);", "perpendicular to the moduli vector"),
        ("rvMukaiModuliFujiki[M,0](zero[RVMukai(B)],v);", "perpendicular to the moduli vector"),
        ("rvMukaiModuliRR[M](v);", "perpendicular to the moduli vector"),
        ("rvMukaiModuliFujiki[M,7](zero[RVMukai(B)],zero[RVMukai(B)]);", "between zero and 2n"),
        ("rvjactorsionbasis D[C,2](neg(o)){e1,e2,f1,f2}; rvMukaiModuliDiv[M](zero[RVMukai(D)]);", "expected type RVMukai(B)"),
        ("rvmukaimoduli Large(rvMukaiIdeal[B,85070591730234615865843651857942052864]); rvMukaiModuliFujiki[Large,0](zero[RVMukai(B)],zero[RVMukai(B)]);", "Fujiki dimension exceeds the work limit"),
        ("rvMukaiModuliRR[M](rvMukai[85070591730234615865843651857942052863,170141183460469231731687303715884105726](zero[RVK3Lattice(B)]));", "Riemann–Roch value exceeds the signed 128-bit range"),
    ] {
        let error = diagnostic(interpret(&format!("{setup}{source}")));
        assert!(error.contains(expected), "{error}");
    }
}

const NONTRIVIAL_ISOTROPIC: &str = r"
    let c=rvK3Lattice[B]{[0,2,0,0,4,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]};
    let v=rvMukai[2,2](c);
    rvmukaimoduli Partner(v);
";

#[test]
fn isotropic_example_checks_the_partner_pairing_and_riemann_roch() {
    assert_eq!(
        stdout(example("hyperelliptic-mukai-isotropic")),
        concat!(
            "rvMukaiModuliLattice[Partner] = n=1, dimension=2, rank=22, signature=(3,19), determinant=-1, discriminant=trivial\n",
            "rvMukaiPair = 1\n",
            "rvMukaiPair = 0\n",
            "rvMukaiQuotientPair = 4\n",
            "rvMukaiModuliDiv[Partner] = 1\n",
            "rvMukaiModuliFujiki[Partner, 1] = 4\n",
            "rvMukaiModuliRR[Partner] = 4\n",
            "rvMukaiPair = 0\n",
            "rvMukaiPair = 1\n",
        )
    );
}

#[test]
fn isotropic_partner_basis_is_an_even_unimodular_integral_kernel() {
    let mut source = format!(
        r"{}{NONTRIVIAL_ISOTROPIC}
        let w=rvMukaiModuliDual[Partner];
        let w0=rvMukai[0,0](rvK3LatticeBasis[B,4]);
        rvMukaiPair(v,w0);
        rvMukaiPair(w0,w0);
        rvMukaiPair(v,w);
        rvMukaiPair(w,w);
        assert rvMukaiModuliGlue[Partner]==zero[RVMukai(B)];
        ",
        marking()
    );
    for i in 0..22 {
        writeln!(
            source,
            "let x{i}=rvMukaiModuliBasis[Partner,{i}]; rvMukaiPair(v,x{i}); rvMukaiPair(w,x{i});"
        )
        .unwrap();
    }
    for i in 0..22 {
        for j in i..22 {
            writeln!(source, "rvMukaiPair(x{i},x{j});").unwrap();
        }
    }
    let output = stdout(interpret(&source));
    let mut values = output.lines().map(|line| {
        line.strip_prefix("rvMukaiPair = ")
            .unwrap()
            .parse::<i128>()
            .unwrap()
    });
    assert_eq!(values.next(), Some(1));
    let old_square = values.next().unwrap();
    assert!(old_square < 0);
    assert_eq!(old_square % 2, 0);
    assert_eq!(values.next(), Some(1));
    assert_eq!(values.next(), Some(0));
    for _ in 0..44 {
        assert_eq!(values.next(), Some(0));
    }
    let mut gram = [[0i128; 22]; 22];
    for i in 0..22 {
        for j in i..22 {
            let entry = values.next().unwrap();
            gram[i][j] = entry;
            gram[j][i] = entry;
        }
        assert_eq!(gram[i][i] % 2, 0);
    }
    assert_eq!(values.next(), None);
    assert_eq!(determinant(gram), -1);
}

#[test]
fn classical_isotropic_axes_give_the_marked_k3_lattice() {
    let mut source = format!(
        r"{}
        let zK=zero[RVK3Lattice(B)];
        let point=rvMukai[0,1](zK);
        let rank=rvMukai[1,0](zK);
        rvmukaimoduli Point(point);
        rvmukaimoduli Rank(rank);
        assert rvMukaiModuliDual[Point]==neg(rank);
        assert rvMukaiModuliDual[Rank]==neg(point);
        let big=rvMukai[0,-170141183460469231731687303715884105728](zK);
        assert rvMukaiQuotient[Point](big)==zero[RVMukaiQuotient(Point)];
        let bigRank=rvMukai[-170141183460469231731687303715884105728,0](zK);
        assert rvMukaiQuotient[Rank](bigRank)==zero[RVMukaiQuotient(Rank)];
        rvMukaiModuliRR[Point](big);
        rvMukaiModuliDiv[Point](big);
        rvMukaiModuliFujiki[Point,1](big,point);
        ",
        marking()
    );
    for i in 0..22 {
        writeln!(
            source,
            "let y{i}=rvMukai[0,0](rvK3LatticeBasis[B,{i}]); assert rvMukaiModuliBasis[Point,{i}]==y{i}; assert rvMukaiModuliBasis[Rank,{i}]==y{i};"
        )
        .unwrap();
    }
    source.push_str(
        r"
        let x=rvMukaiQuotient[Point](y0);
        let y=rvMukaiQuotient[Point](add(y1,point));
        assert add(x,neg(x))==zero[RVMukaiQuotient(Point)];
        assert rvMukaiQuotientLift(add(x,y))==add(y0,y1);
        assert rvMukaiQuotient[Point](add(y0,neg(point)))==x;
        rvMukaiQuotientPair(x,y);
        rvMukaiModuliDiv[Point](add(y,y));
        rvMukaiModuliFujiki[Point,1](x,y);
        ",
    );
    assert_eq!(
        stdout(interpret(&source)),
        concat!(
            "rvMukaiModuliRR[Point] = 2\n",
            "rvMukaiModuliDiv[Point] = 0\n",
            "rvMukaiModuliFujiki[Point, 1] = 0\n",
            "rvMukaiQuotientPair = -1\n",
            "rvMukaiModuliDiv[Point] = 2\n",
            "rvMukaiModuliFujiki[Point, 1] = -1\n",
        )
    );
}

#[test]
fn nontrivial_isotropic_projector_is_shift_invariant_and_orthogonal() {
    let output = stdout(interpret(&format!(
        r"{}{NONTRIVIAL_ISOTROPIC}
        let w=rvMukaiModuliDual[Partner];
        let x=rvMukai[1,0](rvK3LatticeBasis[B,1]);
        let y=rvMukaiModuliBasis[Partner,2];
        let q=rvMukaiQuotient[Partner](x);
        let q2=rvMukaiQuotient[Partner](add(x,neg(add(v,v))));
        let qLift=rvMukaiQuotientLift(q);
        assert q==q2;
        assert rvMukaiQuotient[Partner](qLift)==q;
        assert rvMukaiQuotientLift(rvMukaiQuotient[Partner](add(x,y)))==add(qLift,y);
        rvMukaiPair(qLift,v);
        rvMukaiPair(qLift,w);
        rvMukaiModuliDiv[Partner](x);
        rvMukaiModuliDiv[Partner](q);
        rvMukaiModuliRR[Partner](x);
        rvMukaiModuliRR[Partner](q);
        rvMukaiModuliFujiki[Partner,1](add(x,v),add(y,v));
        rvMukaiQuotientPair(q,rvMukaiQuotient[Partner](y));
        ",
        marking()
    )));
    let entries: Vec<_> = output
        .lines()
        .map(|line| line.split_once(" = ").unwrap().1)
        .collect();
    assert_eq!(&entries[..2], ["0", "0"]);
    assert_eq!(entries.len(), 8);
    assert!(entries[2..].chunks_exact(2).all(|pair| pair[0] == pair[1]));
}

#[test]
fn an_isotropic_vector_can_require_bezout_with_no_unit_entry() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{}
            let c=rvK3Lattice[B]{{[0,6,0,0,20,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}};
            let v=rvMukai[6,10](c);
            rvmukaimoduli Multi(v);
            let w=rvMukaiModuliDual[Multi];
            let a=rvMukaiModuliBasis[Multi,0];
            rvMukaiPair(v,w);
            rvMukaiPair(w,w);
            rvMukaiPair(a,v);
            rvMukaiPair(a,w);
            assert rvMukaiQuotientLift(rvMukaiQuotient[Multi](add(a,v)))==a;
            ",
            marking()
        ))),
        "rvMukaiPair = 1\nrvMukaiPair = 0\nrvMukaiPair = 0\nrvMukaiPair = 0\n"
    );
}

#[test]
fn isotropic_rejections_are_nominal_and_checked_at_the_coordinate_boundary() {
    let setup = format!(
        r"{}{NONTRIVIAL_ISOTROPIC}
        rvmukaimoduli Same(v);
        rvmukaimoduli Positive(rvMukaiIdeal[B,2]);
        let ambient=zero[RVMukai(B)];
        ",
        marking()
    );
    for (source, expected) in [
        ("rvmukaimoduli Divisible(add(v,v));", "requires a primitive Mukai vector"),
        ("show rvMukaiModuliDual[Positive];", "requires an isotropic Mukai moduli lattice"),
        ("show zero[RVMukaiQuotient(Positive)];", "requires an isotropic Mukai moduli lattice"),
        ("show rvMukaiQuotient[Positive](ambient);", "requires an isotropic Mukai moduli lattice"),
        ("show rvMukaiQuotientLift(ambient);", "an RVMukaiQuotient class"),
        ("show rvMukaiQuotient[Partner](rvMukaiModuliDual[Partner]);", "perpendicular to the moduli vector"),
        ("show rvMukaiModuliBasis[Partner,22];", "between 0 and 21"),
        ("rvMukaiQuotientPair(zero[RVMukaiQuotient(Partner)],zero[RVMukaiQuotient(Same)]);", "expected type RVMukaiQuotient(Partner)"),
        ("rvMukaiModuliRR[Partner](zero[RVMukaiQuotient(Same)]);", "expected type RVMukaiQuotient(Partner)"),
        ("rvjactorsionbasis D[C,2](neg(o)){e1,e2,f1,f2}; show rvMukaiQuotient[Partner](zero[RVMukai(D)]);", "expected type RVMukai(B)"),
        (
            "let x=rvMukai[85070591730234615865843651857942052866,0](rvK3Lattice[B]{[0,85070591730234615865843651857942052866,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}); show rvMukaiQuotient[Partner](x);",
            "Mukai coordinate exceeds the signed i128 range",
        ),
        (
            "rvmukaimoduli P(rvMukai[0,1](zero[RVK3Lattice(B)])); let y=rvMukai[0,0](rvK3Lattice[B]{[0,-170141183460469231731687303715884105728,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}); show neg(rvMukaiQuotient[P](y));",
            "K3 lattice coordinate exceeds the signed",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}{source}")));
        assert!(error.contains(expected), "{error}");
    }
}
