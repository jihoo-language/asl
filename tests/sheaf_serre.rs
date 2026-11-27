mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn the_a3_serre_orbit_returns_to_its_simple_after_four_actual_steps() {
    let output = stdout(example("sheaf-serre"));
    assert_eq!(
        output,
        concat!(
            "shserreshift[K1] = 2\nshdim[F1a, 0] = 0\nshdim[F1b, 0] = 1\nshdim[F1c, 0] = 0\n",
            "shserreshift[K2] = 2\nshdim[F2a, 0] = 0\nshdim[F2b, 0] = 0\nshdim[F2c, 0] = 1\n",
            "shserreshift[K3] = 2\nshdim[F3a, 0] = 1\nshdim[F3b, 0] = 1\nshdim[F3c, 0] = 1\n",
            "shserreshift[K4] = 2\nshdim[F4a, 0] = 1\nshdim[F4b, 0] = 0\nshdim[F4c, 0] = 0\n",
        )
    );
}

#[test]
fn signed_bar_faces_and_last_face_projection_are_actual_stalk_maps() {
    let program = r"
        space X {a<b};
        sheaf P[X,5] {a:1,b:1} {a<b:[2]};
        shserre K[P];
        let at_a = shTotEmbed[K,0](shTerm[K__0,[a],0]);
        let at_b = shTotEmbed[K,0](shTerm[K__0,[b],0]);
        assert shTotPart[1](shD(at_a)) ==
            add(shTerm[K__1,[a],0,-1],shTerm[K__1,[a],1,2]);
        assert shTotPart[1](shD(at_b)) == shTerm[K__1,[b],0,2];
        let e = shTerm[K__1,[a],0];
        assert shD(e) == zero[ShC(K__1,1)];
        assert shD(shD(at_a)) == zero[ShTot(K,2)];
        shcohomology I[K,1];
        let p = tab One -> X {pt=>X.a};
        let q = tab One -> X {pt=>X.b};
        shpull A[p,I]; shdim[A,0];
        shpull B[q,I]; shdim[B,0];
        shserreshift[K];
    ";
    assert_eq!(
        stdout(interpret(program)),
        "shdim[A, 0] = 1\nshdim[B, 0] = 0\nshserreshift[K] = 1\n"
    );
}

#[test]
fn a_nonprojective_a4_interval_cycles_back_with_accumulated_shift_three() {
    let mut source = String::from(
        "space X {a<b<c<d};\n\
         sheaf F0[X,5] {a:1,b:1,c:0,d:0} {a<b:[1],b<c:[],c<d:[]};\n",
    );
    let labels = ['a', 'b', 'c', 'd'];
    for point in labels {
        source.push_str(&format!(
            "let at{point} = tab One -> X {{pt=>X.{point}}};\n"
        ));
    }
    let mut expected = Vec::new();
    let mut restrictions = Vec::new();
    for (index, (lo, hi, degree)) in [(1, 2, 2), (2, 3, 2), (0, 2, 3), (1, 3, 2), (0, 1, 3)]
        .into_iter()
        .enumerate()
    {
        let j = index + 1;
        source.push_str(&format!(
            "shserre K{j}[F{index}]; shcohomology F{j}[K{j},{degree}];\n"
        ));
        for (x, label) in labels.into_iter().enumerate() {
            source.push_str(&format!(
                "shpull P{j}{label}[at{label},F{j}]; shdim[P{j}{label},0];\n"
            ));
            expected.push(usize::from((lo..=hi).contains(&x)));
        }
        for x in lo..hi {
            source.push_str(&format!("show shD(shTerm[F{j},[{}],0]);\n", labels[x]));
            restrictions.push(format!("[{}<{}]#0", labels[x], labels[x + 1]));
        }
    }
    let output = stdout(interpret(&source));
    let dims: Vec<_> = output
        .lines()
        .filter(|line| line.starts_with("shdim"))
        .map(|line| line.rsplit_once(" = ").unwrap().1.parse::<usize>().unwrap())
        .collect();
    assert_eq!(dims, expected);
    let columns: Vec<_> = output
        .lines()
        .filter(|line| line.starts_with("ShC(F"))
        .collect();
    assert_eq!(columns.len(), restrictions.len());
    for (column, target) in columns.into_iter().zip(restrictions) {
        assert!(column.contains(&target), "{column} is missing {target}");
    }
}

#[test]
fn cohomology_sheaves_take_the_stalk_quotient_and_handle_empty_sites() {
    let setup = r"
        sheaf X[One,7] {pt:1} {};
        sheaf Y[One,7] {pt:2} {};
        shmap d[X,Y] {pt:[1,-1]};
        shcomplex C {X,d,Y};
        shcohomology H0[C,0];
        shcohomology H1[C,1];
        shcohomology Hbig[C,1000000];
        shdim[H0,0]; shdim[H1,0]; shdim[Hbig,0];
        space Empty {};
        sheaf Z[Empty,2] {} {};
        shserre K[Z];
        shserreshift[K];
        shcohomology E[K,0]; shdim[E,0];
    ";
    assert_eq!(
        stdout(interpret(setup)),
        "shdim[H0, 0] = 0\nshdim[H1, 0] = 1\nshdim[Hbig, 0] = 0\nshserreshift[K] = 0\nshdim[E, 0] = 0\n"
    );
    for (statement, reason) in [
        ("shserreshift[C];", "is not a declared Serre complex"),
        ("shserre K[X];", "sheaf complex 'K' is already defined"),
        ("shserre Q[Missing];", "unknown sheaf 'Missing'"),
        ("shcohomology H1[C,2];", "sheaf 'H1' is already defined"),
    ] {
        let error = diagnostic(interpret(&format!("{setup}\n{statement}")));
        assert!(error.contains(reason), "{statement}\n{error}");
    }
}
