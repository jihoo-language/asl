use std::io::Write;
use std::process::{Command, Output, Stdio};

const ASL: &str = env!("CARGO_BIN_EXE_asl");

fn interpret(source: &str) -> Output {
    let mut child = Command::new(ASL)
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(source.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}

fn example(name: &str) -> Output {
    Command::new(ASL)
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .arg(format!("examples/{name}.asl"))
        .output()
        .unwrap()
}

fn stdout(output: Output) -> String {
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

fn diagnostic(output: Output) -> String {
    assert!(!output.status.success(), "{output:?}");
    let error = String::from_utf8(output.stderr).unwrap();
    assert!(!error.is_empty());
    assert!(!error.contains("panicked"), "{error}");
    error
}

#[test]
fn rotation_descends_and_its_finite_homotopy_class_is_rigid() {
    assert_eq!(
        stdout(example("descent")),
        "Circle.bc\nbetti(Circle) = [1, 1]\nbetti(Circle * Circle) = [1, 2, 1]\nhomotopy = no\n"
    );
}

#[test]
fn a_point_evaluates_opens_and_meets_compute_intersections() {
    assert_eq!(
        stdout(example("observables")),
        "S.hi\nS.lo\nhomotopy = yes (2 steps)\n"
    );
}

#[test]
fn projective_plane_has_a_mod_two_fundamental_class() {
    assert_eq!(
        stdout(example("projective-plane")),
        "betti(RP2) = [1, 1, 1]\n"
    );
}

#[test]
fn double_cover_sections_disagree_on_the_last_overlap() {
    let error = diagnostic(example("rejected/double-cover"));
    for witness in ["double-cover.asl:", "ca0", "ca1"] {
        assert!(error.contains(witness), "{error}");
    }
}

#[test]
fn currying_and_uncurrying_are_inverse_on_the_whole_exponential() {
    let program = r"
        let curry = \k : (S * S) -> S => \x : S => \y : S => k((x, y));
        let uncurry = \h : S -> S -> S => \p : S * S => h(fst(p))(snd(p));
        assert (\k : (S * S) -> S => uncurry(curry(k))) == (\k : (S * S) -> S => k);
        assert (\h : S -> S -> S => curry(uncurry(h))) == (\h : S -> S -> S => h);
        let eval = \p : (S -> S) * S => fst(p)(snd(p));
        show eval(((\x : S => x), S.hi));
    ";
    assert_eq!(stdout(interpret(program)), "S.hi\n");
}

#[test]
fn restriction_and_gluing_compose_on_nested_stars() {
    let program = r"
        space Chain { a < b < c };
        let id = \x : Chain => x;
        let local = res[b](id);
        let germ = res[c](local);
        assert germ(Star(Star(Chain, b), c).c) == Chain.c;
        assert (sew Star(Chain, b) -> Chain { b => res[b](local) }) == local;
        show local(Star(Chain, b).c);
    ";
    assert_eq!(stdout(interpret(program)), "Chain.c\n");
}

#[test]
fn overlap_agreement_is_extensional_for_function_valued_charts() {
    let program = r"
        space Joint { w < j, e < j };
        let id = \x : S => x;
        let eta = \y : S => id(y);
        let bottom = \x : S => S.lo;
        let west = tab Star(Joint, w) -> (S -> S) { w => bottom, j => eta };
        let east = \x : Star(Joint, e) => id;
        let glued = sew Joint -> (S -> S) { e => east, w => west };
        assert res[w](glued) == west;
        show glued(Joint.j)(S.hi);
    ";
    assert_eq!(stdout(interpret(program)), "S.hi\n");
}

#[test]
fn empty_domains_and_discrete_targets_have_their_categorical_counts() {
    let program = r"
        space Void {};
        space Two { a, b };
        let empty = tab Void -> Void {};
        assert empty == (sew Void -> Void {});
        assert empty == (\x : Void => x);
        betti Void;
        betti Void * S;
        betti Void -> Void;
        betti S -> Void;
        betti Two -> Two;
        betti S -> Two;
    ";
    let output = stdout(interpret(program));
    let ranks: Vec<_> = output
        .lines()
        .map(|line| line.split_once(" = ").unwrap().1)
        .collect();
    assert_eq!(ranks, ["[]", "[]", "[1]", "[]", "[4]", "[2]"]);
}

#[test]
fn homotopy_measures_shortest_comparability_zigzags() {
    let program = r"
        space Fence { a < b, c < b, c < d, e < d, isolated };
        let at = \x : Fence => \u : One => x;
        homotopy at(Fence.a) ~ at(Fence.a);
        homotopy at(Fence.a) ~ at(Fence.b);
        homotopy at(Fence.a) ~ at(Fence.e);
        homotopy at(Fence.b) ~ at(Fence.isolated);
        space D { x, y };
        let left = tab D -> S { x => S.hi, y => S.lo };
        let right = tab D -> S { x => S.lo, y => S.hi };
        homotopy left ~ right;
    ";
    assert_eq!(
        stdout(interpret(program)),
        "homotopy = yes (0 steps)\nhomotopy = yes (1 step)\nhomotopy = yes (4 steps)\nhomotopy = no\nhomotopy = yes (2 steps)\n"
    );
}

#[test]
fn individual_maps_work_beyond_the_function_space_enumeration_limit() {
    let domain = "space D { a, b, c, d, e, f, g, h, i, j, k, l, m };";
    let program = format!(
        r"{domain}
        let low = \x : D => S.lo;
        let high = \x : D => S.hi;
        show high(D.m);
        homotopy low ~ high;"
    );
    assert_eq!(
        stdout(interpret(&program)),
        "S.hi\nhomotopy = yes (1 step)\n"
    );

    let error = diagnostic(interpret(&format!("{domain}\nbetti D -> S;")));
    assert!(error.contains(":2:") && error.contains("limit"), "{error}");
}

#[test]
fn discontinuity_points_to_the_reversed_comparison() {
    let error = diagnostic(interpret(
        "\nlet swap = tab S -> S { lo => S.hi, hi => S.lo };",
    ));
    for witness in [":2:", "S.lo", "S.hi"] {
        assert!(error.contains(witness), "{error}");
    }
}

#[test]
fn space_and_star_names_survive_type_checking() {
    for program in [
        r"space Copy { lo < hi }; show (\x : S => x)(Copy.lo);",
        r"let f = res[lo](\x : S => x); show f(S.hi);",
        r"let f = \x : S => x; let bad = sew S -> S { lo => f };",
    ] {
        diagnostic(interpret(program));
    }
}

#[test]
fn empty_binders_still_require_well_typed_bodies() {
    for body in [
        "ghost",
        "S.absent",
        "fst(S.lo)",
        "tab S -> S { hi => S.hi }",
    ] {
        let program = format!(r"space Void {{}}; let f = \x : Void => {body};");
        diagnostic(interpret(&program));
    }
}

#[test]
fn invalid_orders_incomplete_atlases_and_false_assertions_fail() {
    for program in [
        "space Cycle { a < b < c < a };",
        "space D { a, b }; let bad = tab D -> S { a => S.hi, a => S.lo };",
        r"space D { a, b }; let f = \x : D => x; let bad = sew D -> D { a => res[a](f) };",
        "assert S.lo == S.hi;",
        r"homotopy (\x : S => x) ~ (\x : One => S.lo);",
    ] {
        diagnostic(interpret(program));
    }
}
