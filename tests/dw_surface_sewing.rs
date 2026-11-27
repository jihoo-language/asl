mod common;

use std::fmt::Write;

use common::{diagnostic, example, interpret, stdout};

const EXAMPLE: &str = include_str!("../examples/twisted-double-surface-sewing.asl");

fn type_iii() -> &'static str {
    EXAMPLE.split_once("// Genus two").unwrap().0
}

const CYCLIC: &str = "rvgroup G[3];rvdwcoeff P[G,3];
    let k=rvPostnikovCochain[P,3]{0,1,1,1,0,2,2,2};rvdouble D[P](k);";

#[test]
fn runnable_example_counts_handles_and_sews_nontrivial_cut_channels() {
    let output = stdout(example("twisted-double-surface-sewing"));
    assert!(output.contains("charges=[22, 6, 6, 6, 6, 6, 6, 6], verlinde_verified=true"));
    assert!(
        output.contains("charges=[736, 480, 480, 480, 480, 480, 480, 480], verlinde_verified=true")
    );
    assert!(output.contains(
        "charges=[36352, 32256, 32256, 32256, 32256, 32256, 32256, 32256], verlinde_verified=true"
    ));
    assert!(output.contains("dimension=736, channels=[(0, 0, 484), (1, 1, 36), (2, 2, 36), (3, 3, 36), (4, 4, 36), (5, 5, 36), (6, 6, 36), (7, 7, 36)]"));
    assert!(output.contains("output=[5, 2, 5, 7], charge=0, conductor=4, phase=2"));
    assert!(output.contains("RVDoubleSurface(D,2,0){([5, 2, 5, 7],[-1]/1)}"));
    assert!(output
        .contains("neutral_states=22, conductor=4, fourier=true, S4=true, ST3_equals_S2=true"));
    assert!(output.contains("charges=[81, 0, 0], verlinde_verified=true"));
    assert!(
        output.contains("neutral_states=9, conductor=9, fourier=true, S4=true, ST3_equals_S2=true")
    );
    assert!(output.ends_with("RVCycloQ(9){[0, 0, 0, 1] / 1}\n"));
}

fn shuffle(f: impl Fn(usize, usize, usize) -> usize, [a, b, c]: [usize; 3], p: usize) -> usize {
    (p + f(a, b, c) + f(c, a, b) - f(a, c, b)) % p
}

#[test]
fn every_type_iii_neutral_slide_matches_the_independent_charge_and_phase_model() {
    let cocycle = |a: usize, b: usize, c: usize| (a & 1) * ((b >> 1) & 1) * ((c >> 2) & 1);
    let gauge = |a, b| {
        (0..8).fold(0, |n, t| {
            n | ((shuffle(cocycle, [a, b, t], 2) ^ shuffle(cocycle, [b, a, t], 2)) << t)
        })
    };
    let mut script = type_iii().to_owned();
    let mut expected = Vec::new();
    let mut transitions = [[0; 2]; 2];
    let mut counts = [0; 2];
    for a in 0..8 {
        for b in 0..8 {
            for c in 0..8 {
                for d in 0..8 {
                    let q = gauge(a, b);
                    if q != gauge(c, d) {
                        continue;
                    }
                    let m = [a ^ c, b, c, d ^ b];
                    let phase = shuffle(cocycle, [a, c, b], 2) ^ shuffle(cocycle, [d ^ b, b, c], 2);
                    transitions[usize::from(q != 0)][usize::from(gauge(m[0], m[1]) != 0)] += 1;
                    counts[phase] += 1;
                    writeln!(
                        script,
                        "rvDoubleSurfaceChain[D,Slide,0,1]{{{a},{b},{c},{d}}};"
                    )
                    .unwrap();
                    expected.push(format!(
                        "output={m:?}, charge=0, conductor=4, phase={}",
                        2 * phase
                    ));
                }
            }
        }
    }
    let output = stdout(interpret(&script));
    assert_eq!(output.lines().count(), 736);
    for (line, want) in output.lines().zip(expected) {
        assert!(
            line.contains(&want) && line.ends_with("boundary_verified=true"),
            "{want}: {line}"
        );
    }
    assert_eq!(counts, [624, 112]);
    assert_eq!(transitions, [[316, 168], [168, 84]]);
}

#[test]
fn every_twisted_cyclic_three_slide_is_checked_with_its_nonreal_phase() {
    let k = |a: usize, b: usize, c: usize| a * ((b + c) / 3);
    let mut source = CYCLIC.to_owned();
    let mut phases = [0; 3];
    let mut rows = Vec::new();
    for a in 0..3 {
        for b in 0..3 {
            for c in 0..3 {
                for d in 0..3 {
                    let e = (d + 3 - b) % 3;
                    let w = (6 - shuffle(k, [a, c, b], 3) - shuffle(k, [e, b, c], 3)) % 3;
                    phases[w] += 1;
                    let m = [(a + c) % 3, b, c, e];
                    rows.push(format!(
                        "output={m:?}, charge=0, conductor=9, phase={}",
                        3 * w
                    ));
                    writeln!(
                        source,
                        "rvDoubleSurfaceChain[D,Slide,0,1]{{{a},{b},{c},{d}}};"
                    )
                    .unwrap();
                }
            }
        }
    }
    let output = stdout(interpret(&source));
    assert_eq!(output.lines().count(), 81);
    assert_eq!(phases, [58, 7, 16]);
    for (line, want) in output.lines().zip(rows) {
        assert!(
            line.contains(&want) && line.ends_with("boundary_verified=true"),
            "{want}: {line}"
        );
    }
}

#[test]
fn normalized_coboundary_gauge_intertwines_both_torus_moves_and_cross_handle_slides() {
    let source = format!(
        "{CYCLIC}
        let b=rvPostnikovCochain[P,2]{{0,1,0,0}};rvdouble E[P](rvPostnikovGauge(k,b));
        let t=rvDoubleSurfaceBasis[D,0]{{1,2}};
        let v=rvDoubleSurfaceBasis[D,0]{{1,2,2,0}};
        let u=rvDoubleSurfaceGauge[E]{{0,3,0,0}}(v);
        assert rvDoubleSurfaceCoefficient{{1,2,2,0}}(u)==rvCycloRat[1](rvCyclo[9,3]);
        assert rvDoubleSurfaceGauge[D]{{0,-3,0,0}}(u)==v;
        assert rvDoubleSurfaceGauge[E]{{0,3,0,0}}(rvDoubleSurfaceMove[S,0](t))
            ==rvDoubleSurfaceMove[S,0](rvDoubleSurfaceGauge[E]{{0,3,0,0}}(t));
        assert rvDoubleSurfaceGauge[E]{{0,3,0,0}}(rvDoubleSurfaceMove[T,0](t))
            ==rvDoubleSurfaceMove[T,0](rvDoubleSurfaceGauge[E]{{0,3,0,0}}(t));
        assert rvDoubleSurfaceGauge[E]{{0,3,0,0}}(rvDoubleSurfaceMove[Slide,0,1](v))
            ==rvDoubleSurfaceMove[Slide,0,1](u);
        assert rvDoubleSurfaceMove[Swap,0,1](rvDoubleSurfaceMove[Swap,0,1](v))==v;
        assert rvDoubleSurfaceMove[TInv,0](rvDoubleSurfaceMove[T,0](t))==t;
        assert rvDoubleSurfaceMove[SInv,0](rvDoubleSurfaceMove[S,0](t))==t;"
    );
    assert_eq!(stdout(interpret(&source)), "");
    let error = diagnostic(interpret(&format!(
        "{source}show rvDoubleSurfaceGauge[E]{{0,0,0,0}}(v);"
    )));
    assert!(
        error.contains("does not relate the two bulk cocycles"),
        "{error}"
    );
}

#[test]
fn complex_torus_fourier_maps_the_actual_marking_to_the_published_tube_actions() {
    let source=format!("{CYCLIC}
        rvdoubleboundary H[D]{{0}};rvdoubledefects C[H,H];rvdoubletube T[C];
        let color=rvDoubleTubeCentral[T,3];let v=rvDoubleSurfaceFromTube(color);
        assert rvDoubleSurfaceToTube[T](v)==color;
        assert rvDoubleSurfaceMove[T,0](v)==rvDoubleSurfaceScale(rvCycloRat[1](rvCyclo[9,1]),v);
        assert rvDoubleSurfaceToTube[T](rvDoubleSurfaceMove[S,0](v))==rvDoubleTubeWord{{S}}(color);
        assert rvDoubleSurfaceToTube[T](rvDoubleSurfaceMove[T,0](v))==rvDoubleTubeWord{{T[1]}}(color);
        assert rvDoubleSurfaceToTube[T](rvDoubleSurfaceMove[S,0](rvDoubleSurfaceMove[S,0](v)))==rvDoubleTubeCentral[T,7];
        rvDoubleSurfaceModular[D];");
    assert_eq!(stdout(interpret(&source)),
        "rvDoubleSurfaceModular = neutral_states=9, conductor=9, fourier=true, S4=true, ST3_equals_S2=true\n");
}

#[test]
fn tube_fourier_exactly_embeds_and_descends_the_bulk_cyclotomic_field() {
    let source = "rvgroup G[2];rvdwcoeff P[G,2];rvdouble D[P](rvPostnikovK[P]);
        rvdoubleboundary H[D]{0,1}{(2,3)};rvdoubledefects C[H,H];rvdoubletube T[C];
        let v=rvDoubleSurfaceScale(rvCycloRat[1](rvCyclo[4,1]),rvDoubleSurfaceBasis[D,0]{1,1});
        let lift=rvDoubleSurfaceToTube[T](v);
        assert rvDoubleSurfaceFromTube(lift)==v;";
    assert_eq!(stdout(interpret(source)), "");
    let error=diagnostic(interpret(&format!("{source}
        show rvDoubleSurfaceFromTube(rvDoubleTubeScale(rvCycloRat[1](rvCyclo[12,1]),rvDoubleTubeCentral[T,0]));")));
    assert!(error.contains("does not descend"), "{error}");
}

#[test]
fn surface_terms_enforce_charge_nominal_bulk_and_integral_marking_domains() {
    let setup = format!(
        "{}let u=rvDoubleSurfaceBasis[D,0]{{0,2,5,5}};let v=rvDoubleSurfaceBasis[D,5]{{5,2}};
        rvdoubletube T[PP];rvdouble X[P](k);",
        type_iii()
    );
    for (op, want) in [
        ("show rvDoubleSurfaceBasis[D,0]{5,2};", "gauge charge 5"),
        (
            "show rvDoubleSurfaceBasis[D,0]{2,8};",
            "valid pairs of handle holonomies",
        ),
        (
            "show rvDoubleSurfaceSew(v,u);",
            "requires inverse electric charges",
        ),
        (
            "show rvDoubleSurfaceMove[Slide,0,0](u);",
            "valid distinct handle indices",
        ),
        ("show rvDoubleSurfaceCut[3,0](u);", "valid separating cut"),
        (
            "show rvDoubleSurfaceToTube[T](rvDoubleSurfaceBasis[X,0]{0,0});",
            "identical nominal bulk",
        ),
        (
            "show rvDoubleSurfaceToTube[T](u);",
            "requires the neutral torus",
        ),
        (
            "show rvDoubleSurfaceFromTube(rvDoubleTubeMatrixUnit[T,8,0,1]);",
            "tube center",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}{op}")));
        assert!(error.contains(want), "{want}: {error}");
    }
}

#[test]
fn odd_electric_charges_sew_with_the_distinct_inverse_character() {
    let mut source = String::from("rvgroup G{");
    for a in 0..27 {
        source.push('{');
        for b in 0..27 {
            let c =
                (a % 3 + b % 3) % 3 + 3 * ((a / 3 % 3 + b / 3 % 3) % 3) + 9 * ((a / 9 + b / 9) % 3);
            write!(source, "{c},").unwrap();
        }
        source.push_str("},");
    }
    source.push_str("};rvdwcoeff P[G,3];let k=rvPostnikovCochain[P,3]{");
    for a in 1..27 {
        for b in 1..27 {
            for c in 1..27 {
                write!(source, "{},", (a % 3 * (b / 3 % 3) * (c / 9)) % 3).unwrap();
            }
        }
    }
    source.push_str(
        "};rvdouble D[P](k);
        let v=rvDoubleSurfaceBasis[D,1]{1,3};let w=rvDoubleSurfaceBasis[D,2]{3,1};
        let x=rvDoubleSurfaceSew(v,w);let y=rvDoubleSurfaceBasis[D,0]{0,0,0,0};
        let mixed=rvDoubleSurfaceAdd(x,y);
        assert x==rvDoubleSurfaceBasis[D,0]{1,3,3,1};
        assert rvDoubleSurfaceCut[1,0](mixed)==y;
        assert rvDoubleSurfaceCut[1,1](mixed)==x;
        rvDoubleSurfaceSewing[D,1,1];",
    );
    let output = stdout(interpret(&source));
    assert!(output.contains("dimension=26001, channels=[(0, 0, 11025), (1, 2, 576), (2, 1, 576)"));
    assert!(output.ends_with("verified=true\n"));
    let error = diagnostic(interpret(&format!("{source}show rvDoubleSurfaceSew(v,v);")));
    assert!(
        error.contains("requires inverse electric charges"),
        "{error}"
    );
}
