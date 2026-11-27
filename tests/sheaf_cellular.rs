#[allow(dead_code)]
mod common;

use std::fmt::Write;
use std::process::Command;

use common::{diagnostic, interpret, stdout};

fn generated(order: usize, prime: usize, split: bool, cellular: bool) -> String {
    stdout(
        Command::new("awk")
            .current_dir(env!("CARGO_MANIFEST_DIR"))
            .args(["-v", &format!("r={order}"), "-v", &format!("prime={prime}")])
            .args(["-v", if split { "untwisted=1" } else { "untwisted=0" }])
            .args(["-v", if cellular { "cellular=1" } else { "cellular=0" }])
            .args(["-f", "examples/sheaf-transgression.awk"])
            .output()
            .unwrap(),
    )
}

#[test]
fn signed_incidence_uses_restrictions_and_identity_at_equal_anchors() {
    let source = r"
        space X {a<b};
        sheaf F[X,3] {a:1,b:1} {a<b:[2]};
        shresolution R[X,3] {u:(0,a),v:(0,b),e:(1,b)} {e:[u:1,v:-1]} {u:1,v:1};
        shmap id[F,F] {a:[1],b:[1]};
        shcomplex K {F,id,F};
        shcellular C[R,K];
        let u = shCellTerm[C,0,u,0];
        let v = shCellTerm[C,0,v,0];
        let e = shCellTerm[C,0,e,0];
        let upperU = shCellTerm[C,1,u,0];
        let upperV = shCellTerm[C,1,v,0];
        let upperE = shCellTerm[C,1,e,0,-1];
        assert shD(u) == add(shScale[-1](e),upperU);
        assert shD(v) == add(shScale[-1](e),upperV);
        assert shD(e) == upperE;
        assert shD(shD(u)) == zero[ShTot(C,2)];
        assert shD(shD(v)) == zero[ShTot(C,2)];
        assert shED(shEClass[0,1](e)) == shEClass[0,1](upperE);
        show shD(e);
        shhyperdim[C,0];
        shhyperdim[C,1];
        shhyperdim[C,2];
        card (\x : ShTot(C,0) => S.hi);
        shcomplex Ordinary {F};
        shcellular Single[R,Ordinary];
        let c = shhypergen[Single,0,0];
        assert shClass(shRep(c)) == c;
        shedim[Single,2,0,0];
        shedim[Single,2,1,0];
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "ShTot(C, 2){2*q1:e#0}\n",
            "shhyperdim[C, 0] = 0\n",
            "shhyperdim[C, 1] = 0\n",
            "shhyperdim[C, 2] = 0\n",
            "card = 9 / 9\n",
            "shedim[Single, 2, 0, 0] = 1\n",
            "shedim[Single, 2, 1, 0] = 0\n",
        ),
    );
}

#[test]
fn roos_and_cellular_pages_agree_in_all_bidegrees_from_e2() {
    for (r, p) in [(2, 3), (3, 3), (3, 5)] {
        for split in [false, true] {
            let mut dimensions = String::new();
            for page in 2..=r + 1 {
                for horizontal in 0..=r {
                    for vertical in 0..r {
                        writeln!(dimensions, "shedim[K,{page},{horizontal},{vertical}];").unwrap();
                    }
                }
            }
            for degree in 0..=2 * r {
                writeln!(dimensions, "shhyperdim[K,{degree}];").unwrap();
            }
            let roos = stdout(interpret(&(generated(r, p, split, false) + &dimensions)));
            let cellular = stdout(interpret(&(generated(r, p, split, true) + &dimensions)));
            let verify = |output: &str| {
                let observed: Vec<_> = output.lines().collect();
                let hit = observed[1];
                assert_eq!(hit.ends_with("{0}"), split, "r={r}, p={p}: {hit}");
                assert!(hit.starts_with(&format!("ShE(K, {r}, {r}, 0)")), "{hit}");
                let mut expected = vec![format!("shedim[K, 2, 0, {}] = 1", r - 1)];
                expected.push(format!(
                    "shedim[K, {}, 0, {}] = {}",
                    r + 1,
                    r - 1,
                    usize::from(split)
                ));
                expected.push(format!(
                    "shedim[K, {}, {r}, 0] = {}",
                    r + 1,
                    usize::from(split)
                ));
                for page in 2..=r + 1 {
                    for horizontal in 0..=r {
                        for vertical in 0..r {
                            let permanent = (horizontal, vertical) == (0, 0)
                                || (horizontal, vertical) == (r, r - 1);
                            let paired = (horizontal, vertical) == (0, r - 1)
                                || (horizontal, vertical) == (r, 0);
                            let rank = usize::from(permanent || (paired && (split || page <= r)));
                            expected.push(format!(
                                "shedim[K, {page}, {horizontal}, {vertical}] = {rank}"
                            ));
                        }
                    }
                }
                for degree in 0..=2 * r {
                    let rank = usize::from(
                        degree == 0
                            || degree == 2 * r - 1
                            || (split && (degree == r || degree == r - 1)),
                    );
                    expected.push(format!("shhyperdim[K, {degree}] = {rank}"));
                }
                assert!(observed
                    .into_iter()
                    .enumerate()
                    .filter_map(|(i, line)| (i != 1).then_some(line))
                    .eq(expected));
            };
            verify(&roos);
            verify(&cellular);
        }
    }
}

#[test]
fn cellular_resolution_computes_the_twisted_d20_with_split_control() {
    for split in [false, true] {
        let mut source = generated(20, 3, split, true);
        source.push_str(
            "let lift = shELiftTo[20,0](shERep(v2));
             assert shEClass[20,0](lift) == v20;
             assert shELiftTo[20,0](lift) == lift;
             assert shENext(shEPrev(v20)) == v20;",
        );
        let output = stdout(interpret(&source));
        assert_eq!(
            output,
            format!(
                "shedim[K, 2, 0, 19] = 1\nShE(K, 20, 20, 0){{{}}}\n\
                 shedim[K, 21, 0, 19] = {}\nshedim[K, 21, 20, 0] = {}\n",
                if split { "0" } else { "2*g0" },
                usize::from(split),
                usize::from(split),
            ),
        );
    }
}

#[test]
fn resolution_validation_checks_each_principal_stalk_and_boundary() {
    for (declaration, reason) in [
        (
            "{u:(0,b)} {} {u:1}",
            "augmentation is not surjective at 'a'",
        ),
        (
            "{u:(0,a),v:(0,b)} {} {u:1,v:1}",
            "not exact in degree 0 at 'b'",
        ),
        (
            "{u:(0,a),v:(0,b),e:(1,b)} {e:[u:1]} {u:1,v:1}",
            "augmentation of the boundary of 'e' is nonzero",
        ),
        (
            "{u:(0,a),v:(0,b),e:(1,b),f:(2,b)} {e:[u:1,v:-1],f:[e:1]} {u:1,v:1}",
            "boundary squared is nonzero at 'f'",
        ),
        (
            "{u:(0,b),e:(1,a)} {e:[u:1]} {u:1}",
            "must have degree 0 and lie at or below",
        ),
        ("{u:(0,a),e:(2,a)} {e:[u:1]} {u:1}", "must have degree 1"),
        (
            "{u:(0,a),u:(0,b)} {} {u:1}",
            "duplicate sheaf resolution generator 'u'",
        ),
        (
            "{u:(0,a)} {} {u:1,u:2}",
            "duplicate sheaf resolution augmentation entry 'u'",
        ),
        (
            "{u:(0,a),e:(1,b)} {e:[],e:[]} {u:1}",
            "duplicate sheaf resolution boundary 'e'",
        ),
        (
            "{u:(0,a)} {u:[]} {u:1}",
            "a degree-zero generator has no boundary",
        ),
        (
            "{u:(0,a),e:(1,b)} {} {u:1,e:0}",
            "only accepts degree-zero generators",
        ),
        (
            "{u:(0,a)} {} {v:1}",
            "unknown sheaf resolution generator 'v'",
        ),
    ] {
        let source = format!("space X {{a<b}}; shresolution R[X,3] {declaration};");
        let error = diagnostic(interpret(&source));
        assert!(error.contains(reason), "{declaration}\n{error}");
    }
    let source = r"
        space X {a};
        sheaf F[X,2] {a:1} {};
        shresolution R[X,3] {u:(0,a)} {} {u:1};
        shcomplex K {F};
        shcellular C[R,K];
    ";
    let error = diagnostic(interpret(source));
    assert!(error.contains("same named base and prime field"), "{error}");
}

#[test]
fn empty_site_and_missing_coordinates_remain_checked() {
    let empty = r"
        space Empty {};
        sheaf F[Empty,5] {} {};
        shcomplex K {F};
        shresolution R[Empty,5] {} {} {};
        shcellular C[R,K];
        assert shD(zero[ShTot(C,0)]) == zero[ShTot(C,1)];
        shhyperdim[C,0];
        shedim[C,12345,7,2];
    ";
    assert_eq!(
        stdout(interpret(empty)),
        "shhyperdim[C, 0] = 0\nshedim[C, 12345, 7, 2] = 0\n"
    );
    let source = r"
        space X {a};
        sheaf F[X,5] {a:1} {};
        shcomplex K {F};
        shresolution R[X,5] {u:(0,a)} {} {u:1};
        shcellular C[R,K];
    ";
    for (statement, reason) in [
        (
            "show shCellTerm[C,0,u,1];",
            "stalk coordinate 1 is out of range",
        ),
        (
            "show shCellTerm[C,1,u,0];",
            "complex stalk degree 1 is out of range",
        ),
        (
            "show shCellTerm[K,0,u,0];",
            "requires a declared cellular total complex",
        ),
        (
            "show shTotPart[0](shCellTerm[C,0,u,0]);",
            "requires ordinary Roos cochains",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{source}{statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
}

#[test]
fn shared_site_and_transport_keep_high_pages_within_the_default_storage() {
    for (r, p) in [(32, 3), (32, 5), (64, 3), (64, 5)] {
        for split in [false, true] {
            let output = stdout(interpret(&generated(r, p, split, true)));
            let coefficient = if split {
                "0".into()
            } else {
                format!("{}*g0", p - 1)
            };
            assert_eq!(
                output,
                format!(
                    "shedim[K, 2, 0, {q}] = 1\nShE(K, {r}, {r}, 0){{{coefficient}}}\n\
                     shedim[K, {next}, 0, {q}] = {rank}\nshedim[K, {next}, {r}, 0] = {rank}\n",
                    q = r - 1,
                    next = r + 1,
                    rank = usize::from(split),
                ),
            );
        }
    }
}
