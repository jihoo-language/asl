mod common;

use common::{diagnostic, example, interpret, stdout};

const F23: &str = include_str!("../examples/hyperelliptic-k3-pencil-constant-field.asl");

fn f23() -> &'static str {
    F23.split_once("rvK3PencilLPower[P,1]").unwrap().0
}

fn list<'a>(line: &'a str, name: &str) -> &'a str {
    line.split_once(name).unwrap().1.split_once(']').unwrap().0
}

fn field<'a>(line: &'a str, name: &str) -> &'a str {
    line.split_once(name).unwrap().1.split(',').next().unwrap()
}

#[test]
fn f23_exact_powered_l_factors_leading_values_and_artin_tate_products() {
    let output = stdout(example("hyperelliptic-k3-pencil-constant-field"));
    let powers: Vec<_> = output
        .lines()
        .filter(|line| line.starts_with("rvK3PencilLPower"))
        .collect();
    assert_eq!(powers.len(), 5);
    // Independent powers of the three rational eigenvalues and the two
    // roots with sum -18 and product 529 give these coefficient rows.
    let independent = [
        (1, 23, 2, 18, 16, "[1, -5, -414, -9522, -60835, 6436343", "128/23", 256),
        (2, 529, 3, 20, 17, "[1, -853, -45494, 24066326, 126274613317, -41426511213649", "1792/529", 1792),
        (3, 12167, 2, 18, 16, "[1, -34901, 276604578, 3365447900526, -62862029037720163, 266635235464391245607", "3200/12167", 6400),
        (4, 279841, 3, 20, 17, "[1, -818597, 295676082826, -82742290694110666, 17939245816178538709637, -1716155831334586342923895201", "580608/279841", 580608),
        (6, 148035889, 3, 20, 17, "[1, -664870645, 185701028847633610, -27490416893674087002629290, 2156940707939390609882800355768005, -71094348791151363024389554286420996798449", "75308800/148035889", 75308800),
    ];
    for (line, (power, q, rank, ns_rank, trivial, coefficients, leading, product)) in
        powers.into_iter().zip(independent)
    {
        assert!(line.starts_with(&format!("rvK3PencilLPower[P, {power}] = Q={q}, ")));
        assert_eq!(list(line, "coefficients="), coefficients);
        let root_number = if power % 2 == 1 { 1 } else { -1 };
        assert!(line.contains(&format!("NS rank={ns_rank}, trivial rank={trivial}, Mordell-Weil rank={rank}, root number={root_number}, L*={leading}")), "{line}");
        let (tamagawa, orbits) = if power % 2 == 1 { (512, 2) } else { (1024, 1) };
        assert!(line.contains(&format!("Tamagawa product={tamagawa}, component orbit product={orbits}, arithmetic trivial |det|=1024, signed NS discriminant times Brauer order=-{product}")), "{line}");
    }
    assert!(output.contains("rvK3PencilArtinTateRatio[P, 1, 2] = NS ranks=[18, 20], Artin-Tate magnitudes=[256, 1792], NS rank jump"));
    assert!(output.contains("rvK3PencilArtinTateRatio[P, 1, 3] = NS ranks=[18, 18], Artin-Tate magnitudes=[256, 6400], integral NS lattice fixed, Brauer and Sha order ratio=25"));
    assert!(output.contains("rvK3PencilArtinTateRatio[P, 2, 4] = NS ranks=[20, 20], Artin-Tate magnitudes=[1792, 580608], integral NS lattice fixed, Brauer and Sha order ratio=324"));
    assert!(output.contains("rvK3PencilArtinTateRatio[P, 2, 6] = NS ranks=[20, 20], Artin-Tate magnitudes=[1792, 75308800], integral NS lattice fixed, Brauer and Sha order ratio=42025"));
    assert!(output.contains("rvK3PencilShaPower[P, 3] = base Sha order=4, NS ranks=[18, 18], Artin-Tate magnitudes=[256, 6400], integral NS lattice fixed, Sha order ratio=25, extension Sha order=100, Sha=(Z/10)^2"));
}

#[test]
fn direct_quadratic_base_change_agrees_with_powered_l_and_bsd_data() {
    let source = r"
        rvunram F[5,1]{0,1};
        let z=rvUnramElem[F]{};
        let o=rvUnramElem[F]{1};
        let two=rvUnramElem[F]{2};
        let three=rvUnramElem[F]{3};
        rvhyperelliptic C[F]{z,neg(o),z,z,z,o};
        let r0=rvJacPoint[C](z,z);
        let r1=rvJacPoint[C](o,z);
        let r2=rvJacPoint[C](two,z);
        let r3=rvJacPoint[C](three,z);
        rvjactorsionbasis B[C,2](neg(o)){r0,add(add(r0,r1),r2),r1,add(add(r0,r1),r3)};
        rvk3pencil P[B,0,1];
        rvK3PencilLPower[P,2];
        rvunram Ext[5,1]{2,0,1};
        rvjacbase L[C,Ext];
        let zz=rvUnramElem[Ext]{};
        let oo=rvUnramElem[Ext]{1};
        let tw=rvUnramElem[Ext]{2};
        let th=rvUnramElem[Ext]{3};
        let e0=rvJacPoint[L](zz,zz);
        let e1=rvJacPoint[L](oo,zz);
        let e2=rvJacPoint[L](tw,zz);
        let e3=rvJacPoint[L](th,zz);
        rvjactorsionbasis E[L,2](neg(oo)){e0,add(add(e0,e1),e2),e1,add(add(e0,e1),e3)};
        rvk3pencil Q[E,0,1];
        rvK3PencilLFunction[Q]; rvK3PencilBSD[Q];
        ";
    let output = stdout(interpret(source));
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(lines.len(), 3);
    assert!(lines[0].contains("Q=25"));
    assert_eq!(
        list(lines[0], "coefficients="),
        list(lines[1], "coefficients=")
    );
    for name in ["NS rank=", "trivial rank=", "Mordell-Weil rank="] {
        assert_eq!(field(lines[0], name), field(lines[1], name), "{name}");
    }
    for name in [
        "root number=",
        "L*=",
        "Tamagawa product=",
        "component orbit product=",
        "signed NS discriminant times Brauer order=",
        "BSD quotient=",
    ] {
        assert_eq!(field(lines[0], name), field(lines[2], name), "{name}");
    }
}

#[test]
fn large_extension_uses_base_data_and_rank_jumps_keep_only_certified_orders() {
    let output = stdout(interpret(&format!("{} rvK3PencilLPower[P,12]; rvK3PencilShaPower[P,2]{{p,q}}; rvK3PencilArtinTateRatio[P,3,3]; rvK3PencilArtinTateRatio[P,4,12]; rvK3PencilArtinTateRatio[P,2,12];",f23())));
    assert!(output.contains(
        "rvK3PencilLPower[P, 12] = Q=21914624432020321, coefficients=[1, -70650916887448805"
    ));
    assert!(output.contains("L*=38922205272652800/21914624432020321"));
    assert!(output.contains("Brauer order=-38922205272652800"));
    assert!(output.contains("rvK3PencilShaPower[P, 2] = base Sha order=4, NS ranks=[18, 20], Artin-Tate magnitudes=[256, 1792], NS rank jump, extension Sha order unresolved"));
    assert!(output.contains("rvK3PencilArtinTateRatio[P, 3, 3] = NS ranks=[18, 18], Artin-Tate magnitudes=[6400, 6400], integral NS lattice fixed, Brauer and Sha order ratio=1"));
    assert!(output.contains("rvK3PencilArtinTateRatio[P, 4, 12] = NS ranks=[20, 20], Artin-Tate magnitudes=[580608, 38922205272652800], integral NS lattice fixed, Brauer and Sha order ratio=67036977225"));
    assert!(output.contains("rvK3PencilArtinTateRatio[P, 2, 12] = NS ranks=[20, 20], Artin-Tate magnitudes=[1792, 38922205272652800], integral NS lattice fixed, Brauer and Sha order ratio=21719980620900"));
}

#[test]
fn invalid_or_uncertified_extension_inputs_have_clear_diagnostics() {
    for (operation, expected) in [
        (
            "rvK3PencilLPower[P,0];",
            "constant-field degree must be positive",
        ),
        (
            "rvK3PencilArtinTateRatio[P,0,2];",
            "constant-field degree must be positive",
        ),
        (
            "rvK3PencilArtinTateRatio[P,2,0];",
            "constant-field degree must be positive",
        ),
        (
            "rvK3PencilArtinTateRatio[P,2,3];",
            "must satisfy m dividing n",
        ),
        (
            "rvK3PencilLPower[P,1000];",
            "checked 128-bit arithmetic range",
        ),
        (
            "rvK3PencilShaPower[P,3]{add(p,p),q};",
            "requires an exact base Sha order",
        ),
        ("rvK3PencilShaPower[P,3]{p};", "regulator needs 2 sections"),
        (
            "rvk3pencil Other[B,14,0]; rvK3PencilShaPower[P,3]{p,zero[RVK3PencilSection(Other)]};",
            "expected type RVK3PencilSection(P)",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{}\n{operation}", f23())));
        assert!(error.contains(expected), "{operation}: {error}");
    }
}
