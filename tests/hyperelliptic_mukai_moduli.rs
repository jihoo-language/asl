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
    let gram = read_symmetric_gram::<23>(&mut values);
    assert_eq!(values.next(), None);
    assert_eq!(determinant(gram), 52);
}

fn read_symmetric_gram<const N: usize>(values: &mut impl Iterator<Item = i128>) -> [[i128; N]; N] {
    let mut gram = [[0; N]; N];
    for i in 0..N {
        let (prefix, later_rows) = gram.split_at_mut(i + 1);
        let row = &mut prefix[i];
        row[i] = values.next().unwrap();
        for (entry, transpose_row) in row[i + 1..].iter_mut().zip(later_rows) {
            *entry = values.next().unwrap();
            transpose_row[i] = *entry;
        }
    }
    gram
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
    let gram = read_symmetric_gram::<22>(&mut values);
    for (i, row) in gram.iter().enumerate() {
        assert_eq!(row[i] % 2, 0);
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

#[test]
fn ogrady_example_checks_exceptional_gluing_fujiki_and_rr() {
    assert_eq!(
        stdout(example("hyperelliptic-ogradys-tenfold")),
        concat!(
            "rvOG10Lattice[M] = n=5, dimension=10, rank=24, signature=(3,21), determinant=-3, discriminant=Z/3, q=-2/3 mod 2, index=2\n",
            "rvOG10Pair = -6\n",
            "rvOG10Pair = -3\n",
            "rvOG10Div = 3\n",
            "rvOG10Fujiki[0] = -7348320\n",
            "rvOG10Fujiki[1] = -3674160\n",
            "rvOG10RR = 6\n",
            "rvOG10RR = 0\n",
            "rvOG10RR = -252\n",
        )
    );
}

#[test]
fn ogrady_gram_determinant_and_pairing_hold_in_a_nontrivial_marking() {
    let mut source = format!(
        r"{}
        let h=rvK3FromKummerDiv(rvKummerHyperplane[B]);
        let e=rvK3FromKummerDiv(rvKummerNodeDiv[B,0]);
        let odd=rvK3LatticeBasis[B,1];
        let w=rvMukai[0,0](add(add(h,odd),e));
        rvmukaimoduli M(w);
        let g=rvMukaiModuliGlue[M];
        let sigma=rvOG10Exc[M];
        let gen=rvOG10Basis[M,23];
        assert add(gen,gen)==add(rvOG10[M,0](g),sigma);
        rvMukaiPair(w,w);
        rvMukaiPair(g,g);
        rvOG10Pair(gen,gen);
        rvOG10Div(sigma);
        ",
        marking(),
    );
    for i in 0..24 {
        writeln!(source, "let x{i}=rvOG10Basis[M,{i}];").unwrap();
        if i < 23 {
            writeln!(source, "let u{i}=rvMukaiModuliBasis[M,{i}]; rvMukaiPair(g,u{i}); rvOG10Pair(gen,x{i}); rvOG10Pair(sigma,x{i});").unwrap();
        }
    }
    for i in 0..24 {
        for j in i..24 {
            writeln!(source, "rvOG10Pair(x{i},x{j});").unwrap();
        }
    }
    let output = stdout(interpret(&source));
    let mut values = output
        .lines()
        .map(|line| line.split_once(" = ").unwrap().1.parse::<i128>().unwrap());
    assert_eq!(values.next(), Some(2));
    let glue_square = values.next().unwrap();
    let gen_square = values.next().unwrap();
    assert_eq!(glue_square.rem_euclid(8), 6);
    assert_eq!(4 * gen_square, glue_square - 6);
    assert_eq!(values.next(), Some(3));

    let mut nontrivial_cross_pair = false;
    for _ in 0..23 {
        let glue = values.next().unwrap();
        let gen = values.next().unwrap();
        assert_eq!(glue, 2 * gen);
        assert_eq!(values.next(), Some(0));
        nontrivial_cross_pair |= gen != 0;
    }
    assert!(nontrivial_cross_pair);

    let mut upper = [[0i128; 24]; 24];
    for (i, row) in upper.iter_mut().enumerate() {
        for entry in row.iter_mut().skip(i) {
            *entry = values.next().unwrap();
        }
        assert_eq!(row[i] % 2, 0);
    }
    assert_eq!(values.next(), None);
    let gram =
        std::array::from_fn(|i: usize| std::array::from_fn(|j: usize| upper[i.min(j)][i.max(j)]));
    assert_eq!(determinant::<24>(gram), -3);
}

#[test]
fn ogrady_pairing_is_bilinear_and_retains_wide_exact_cancellation() {
    let large = 1i128 << 125;
    let doubled = 2 * large;
    assert_eq!(
        stdout(interpret(&format!(
            r"{}
            rvmukaimoduli M(rvMukaiIdeal[B,2]);
            let h=rvK3Lattice[B]{{[0,{doubled},0,0,{doubled},0],[{doubled},0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}};
            let x=rvOG10[M,{large}](rvMukai[0,0](h));
            let exc=rvOG10Exc[M];
            rvOG10Pair(x,x);
            rvOG10Pair(x,neg(x));
            rvOG10Fujiki[5](x,x);
            rvOG10RR(x);
            let a=rvOG10Basis[M,23];
            let b=rvOG10Basis[M,0];
            rvOG10Pair(exc,add(a,b));
            rvOG10Pair(exc,a);
            rvOG10Pair(exc,b);
            rvOG10Div(zero[RVOG10(M)]);
            ",
            marking()
        ))),
        concat!(
            "rvOG10Pair = 0\n",
            "rvOG10Pair = 0\n",
            "rvOG10Fujiki[5] = 0\n",
            "rvOG10RR = 6\n",
            "rvOG10Pair = -3\n",
            "rvOG10Pair = -3\n",
            "rvOG10Pair = 0\n",
            "rvOG10Div = 0\n",
        )
    );
}

#[test]
fn ogrady_requires_square_two_perpendicular_and_nominal_inputs() {
    let setup = format!(
        r"{}
        let w=rvMukaiIdeal[B,2];
        rvmukaimoduli M(w);
        rvmukaimoduli Same(w);
        rvmukaimoduli SquareFour(rvMukaiIdeal[B,3]);
        rvmukaimoduli Isotropic(rvMukaiIdeal[B,1]);
        let zeroM=zero[RVMukai(B)];
        ",
        marking()
    );
    for (source, expected) in [
        (
            "show zero[RVOG10(Isotropic)];",
            "requires a Mukai moduli lattice of square two",
        ),
        (
            "rvOG10Lattice[SquareFour];",
            "requires a Mukai moduli lattice of square two",
        ),
        (
            "show rvOG10Exc[Isotropic];",
            "requires a Mukai moduli lattice of square two",
        ),
        ("show rvOG10Basis[M,24];", "between 0 and 23"),
        (
            "show rvOG10[M,0](w);",
            "requires a Mukai vector perpendicular",
        ),
        (
            "rvOG10Pair(zero[RVOG10(M)],zero[RVOG10(Same)]);",
            "expected type RVOG10(M)",
        ),
        ("rvOG10Div(zeroM);", "an RVOG10 class"),
        (
            "rvOG10Fujiki[11](zero[RVOG10(M)],zero[RVOG10(M)]);",
            "between zero and 2n",
        ),
        (
            "rvjactorsionbasis D[C,2](neg(o)){e1,e2,f1,f2}; show rvOG10[M,0](zero[RVMukai(D)]);",
            "expected type RVMukai(B)",
        ),
        (
            "show add(rvOG10[M,170141183460469231731687303715884105727](zeroM),rvOG10Basis[M,23]);",
            "OG10 gluing coefficient exceeds",
        ),
        (
            "show neg(rvOG10[M,-170141183460469231731687303715884105728](zeroM));",
            "OG10 gluing coefficient exceeds",
        ),
        (
            "let x=rvOG10[M,170141183460469231731687303715884105727](zeroM); rvOG10Pair(x,x);",
            "OG10 pairing exceeds",
        ),
        (
            "rvOG10RR(rvOG10[M,100000000000000000000000000](zeroM));",
            "Riemann–Roch value exceeds",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}{source}")));
        assert!(error.contains(expected), "{error}");
    }
}
