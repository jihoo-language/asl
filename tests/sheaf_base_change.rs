mod common;

use std::process::Command;

use common::{diagnostic, example, interpret, stdout};

fn generated(kind: &str, dimension: usize, prime: u32) -> String {
    stdout(
        Command::new("awk")
            .current_dir(env!("CARGO_MANIFEST_DIR"))
            .args([
                "-v",
                &format!("kind={kind}"),
                "-v",
                &format!("dim={dimension}"),
                "-v",
                &format!("prime={prime}"),
                "-v",
                "witness=1",
                "-f",
                "examples/sheaf-base-change.awk",
            ])
            .output()
            .unwrap(),
    )
}

#[test]
fn duplicate_images_are_degenerate_and_the_actual_roos_map_is_natural() {
    assert_eq!(
        stdout(example("sheaf-base-change")),
        concat!(
            "shmaprank[Degree0, lo] = 1\n",
            "shmaprank[Degree0, hi] = 1\n",
            "shbcuniversal = exact\n",
        )
    );
}

#[test]
fn connected_strict_squares_first_fail_at_degree_six() {
    for (dimension, prime) in [(2, 3), (5, 3), (6, 5)] {
        let input = generated("sphere", dimension, prime);
        let result = stdout(interpret(&input));
        let lines: Vec<_> = result.lines().collect();
        assert_eq!(lines[0], "shbcrank = [1, 1, 1, 0, 0]");
        assert_eq!(lines[1], "shbcrank = [0, 1, 0, 0, 1]");
        assert_eq!(lines[2], format!("shbccomma = [({dimension}, 1)]"));
        assert!(lines[3].starts_with(&format!(
            "shbcuniversal = fails at (pt, top) in degree {dimension} with cycle"
        )));
        assert_eq!(
            lines[3].split_once("with cycle ").unwrap().1,
            lines[4].strip_prefix("shbcwitness = ").unwrap()
        );
        assert_eq!(lines[4].matches(", [\"").count(), 1 << (dimension + 1));
    }
}

#[test]
fn projective_plane_comma_is_acyclic_at_odd_primes_only() {
    for prime in [2, 3, 5, 7, 11, 13] {
        let input = generated("rp2", 2, prime);
        let result = stdout(interpret(&input));
        let lines: Vec<_> = result.lines().collect();
        assert_eq!(lines[0], "shbcrank = [1, 1, 1, 0, 0]");
        if prime == 2 {
            assert_eq!(
                &lines[1..4],
                [
                    "shbcrank = [0, 1, 0, 0, 1]",
                    "shbcrank = [0, 1, 0, 0, 1]",
                    "shbccomma = [(1, 1), (2, 1)]",
                ]
            );
            assert!(
                lines[4].starts_with("shbcuniversal = fails at (pt, top) in degree 1 with cycle ")
            );
        } else {
            assert_eq!(
                &lines[1..],
                [
                    "shbcrank = [0, 0, 0, 0, 0]",
                    "shbcrank = [0, 0, 0, 0, 0]",
                    "shbccomma = []",
                    "shbcuniversal = exact",
                ]
            );
        }
    }
}

#[test]
fn empty_comma_is_detected_as_an_actual_kernel_in_degree_zero() {
    let output = stdout(interpret(
        r"
        space Empty {};
        let f=tab One->S {pt=>S.hi};
        let u=tab Empty->One {};
        let v=tab One->S {pt=>S.lo};
        let fp=tab Empty->One {};
        shbcsquare Q[f,u,v,fp];
        sheaf F[One,3] {pt:1} {};
        shbc C[Q,F];
        shbcrank[C,pt,0];
        shbckernel[C,pt,0];
        shbccokernel[C,pt,0];
        shbccomma[Q,3,pt,pt];
        shbcwitness[Q,3,pt,pt,-1,0];
        shbcuniversal[Q,3];
        assert shBCApply[C](shhypergen[C__from,0,0])==zero[ShHyper(C__to,0)];
    ",
    ));
    assert_eq!(
        output,
        concat!(
            "shbcrank = [1, 0, 0, 1, 0]\nshbckernel = [(0, 1)]\nshbccokernel = []\n",
            "shbccomma = [(-1, 1)]\nshbcwitness = [(1, [])]\n",
            "shbcuniversal = fails at (pt, pt) in degree -1 with cycle [(1, [])]\n",
        )
    );
}

#[test]
fn disconnected_nonempty_comma_has_a_cokernel_witness_in_degree_zero() {
    let mut input = generated("sphere", 0, 3);
    input.push_str("shbccokernel[B,pt,0];");
    let result = stdout(interpret(&input));
    let lines: Vec<_> = result.lines().collect();
    assert_eq!(lines[0], "shbcrank = [1, 2, 1, 0, 1]");
    assert_eq!(lines[1], "shbccomma = [(0, 1)]");
    assert!(lines[2].starts_with("shbcuniversal = fails at (pt, top) in degree 0 with cycle "));
    assert_eq!(lines[4], "shbccokernel = [(0, 1)]");
}

#[test]
fn signed_coefficient_complexes_are_pulled_before_totalization() {
    let output = stdout(interpret(
        r"
        space X {a<b}; space XP {x<y<z};
        let f=tab X->S {a=>S.lo,b=>S.hi};
        let u=tab XP->X {x=>X.a,y=>X.a,z=>X.b};
        let v=\s:S=>s;
        let fp=tab XP->S {x=>S.lo,y=>S.lo,z=>S.hi};
        shbcsquare Q[f,u,v,fp];
        sheaf F[X,5] {a:2,b:2} {a<b:[2,0,0,1]};
        shmap nonzero[F,F] {a:[1,0,0,0],b:[1,0,0,0]};
        shcomplex K {F,nonzero,F};
        shbccomplex C[Q,K];
        let x=shTotEmbed[C__from,0](shTerm[C__from__0,[lo,hi],0]);
        assert shBCApply[C](shD(x))==shD(shBCApply[C](x));
        shbchom H1[C,1];
        shmaprank[H1,lo];
        shmaprank[H1,hi];
        shbcrank[C,lo,0];
        shbcrank[C,lo,1];
        shbcrank[C,lo,1000000000];
    ",
    ));
    assert_eq!(
        output,
        concat!(
            "shmaprank[H1, lo] = 1\nshmaprank[H1, hi] = 1\n",
            "shbcrank = [1, 1, 1, 0, 0]\nshbcrank = [1, 1, 1, 0, 0]\n",
            "shbcrank = [0, 0, 0, 0, 0]\n",
        )
    );
}

#[test]
fn raw_base_change_maps_compose_with_the_canonical_layer_identifications() {
    let output = stdout(interpret(
        r"
        space X {a<b}; space Y {a0<a1<b0}; space Z {aa0<aa1<bb0<bb1};
        let f=\x:X=>One.pt; let g=\y:Y=>One.pt; let h=\z:Z=>One.pt;
        let id=\y:One=>y;
        let u=tab Y->X {a0=>X.a,a1=>X.a,b0=>X.b};
        let v=tab Z->Y {aa0=>Y.a0,aa1=>Y.a1,bb0=>Y.b0,bb1=>Y.b0};
        let uv=\z:Z=>u(v(z));
        shbcsquare Q1[f,u,id,g]; shbcsquare Q2[g,v,id,h]; shbcsquare Q[f,uv,id,h];
        sheaf F[X,3] {a:1,b:1} {a<b:[2]}; shpull P[u,F];
        shbc A[Q1,F]; shbc B[Q2,P]; shbc C[Q,F];
        shbcmap M[A,1]; shbcmap N[B,1]; shbcmap O[C,1];
        shmap middle[A__to__1,B__from__1] {pt:[1,0,0,0,1,0,0,0,1]};
        shmap start[A__from__1,C__from__1] {pt:[1]};
        shmap end[C__to__1,B__to__1] {
           pt:[1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,
               0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1]};
        assert shExtCompose(shMapC[N],shExtCompose(shMapC[middle],shMapC[M]))==
               shExtCompose(shMapC[end],shExtCompose(shMapC[O],shMapC[start]));
    ",
    ));
    assert!(output.is_empty(), "{output}");
}

#[test]
fn invalid_squares_and_witness_grades_are_diagnosed() {
    let setup = include_str!("../examples/sheaf-base-change.asl");
    for (tail, reason) in [
        ("shbcsquare Bad[f,u,fp,v];", "compatible endpoints"),
        (
            "let w=\\s:S=>S.hi; shbcsquare Bad[f,u,w,fp];",
            "square does not commute",
        ),
        ("shbc Bad[Q,B__to__0];", "source of f"),
        (
            "shbccomma[Q,5,hi,a];",
            "outside the pulled-back principal open",
        ),
        ("shbcwitness[Q,5,lo,a,-2,0];", "starts in degree -1"),
        ("shbcwitness[Q,5,lo,a,0,0];", "no comma homology generator"),
        (
            "show shBCApply[B](zero[ShTot(B__to,0)]);",
            "requires ShTot or ShHyper from",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{setup}{tail}")));
        assert!(error.contains(reason), "{tail}\n{error}");
    }
}
