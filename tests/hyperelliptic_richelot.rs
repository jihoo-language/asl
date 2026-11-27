mod common;

use common::{diagnostic, example, interpret, stdout};

const F7: &str = r"
    rvunram F[7,1]{0,1};
    let z=rvUnramElem[F]{};
    let o=rvUnramElem[F]{1};
    let t=rvUnramElem[F]{2};
    let h=rvUnramElem[F]{4};
    let e=rvUnramElem[F]{5};
    let s=rvUnramElem[F]{6};
    rvhyperelliptic C[F]{z,h,t,e,t,o};
";

const SPLITTING: &str = "{[],[1]}{[1],[],[1]}{[4],[2],[1]}";

#[test]
fn quintic_source_and_dual_have_exact_equations_and_matching_weil_polynomials() {
    assert_eq!(
        stdout(example("hyperelliptic-richelot")),
        concat!(
            "rvRichelotData[D] = source=C, determinant=[3], f=[[5], [6], [1], [3], [4], [5], [4]], factors=[[[4], [2], [3]], [[3], [0], [1]], [[1], [0], [6]]]\n",
            "rvHWWPoints[C, 1] = 10, trace=-2\n",
            "rvHWWPoints[D, 1] = 10, trace=-2\n",
            "rvHWWWeil[C] = [1, 2, -2, 14, 49]\n",
            "rvHWWWeil[D] = [1, 2, -2, 14, 49]\n",
            "rvRichelotData[E] = source=D, determinant=[1], f=[[0], [2], [1], [6], [1], [4]], factors=[[[0], [1]], [[5], [0], [5]], [[6], [3], [5]]]\n",
            "rvHWWWeil[E] = [1, 2, -2, 14, 49]\n",
        )
    );
}

#[test]
fn sextic_source_and_permuted_or_iterated_splittings() {
    assert_eq!(
        stdout(interpret(&format!(
            r"{F7}
            rvrichelot D[C]{SPLITTING};
            rvrichelot P[C]{{[4],[2],[1]}}{{[1],[],[1]}}{{[],[1]}};
            rvrichelotdual E[D];
            rvrichelotdual T[E];
            rvhyperelliptic H[F]{{z,o,h,e,e,h,o}};
            rvrichelot Q[H]{{[],[1],[1]}}{{[1],[],[1]}}{{[1],[3],[1]}};
            rvRichelotData[P];
            rvRichelotData[T];
            rvRichelotData[Q];
            rvHWWWeil[H];
            rvHWWWeil[Q];
            "
        ))),
        concat!(
            "rvRichelotData[P] = source=C, determinant=[4], f=[[5], [6], [1], [3], [4], [5], [4]], factors=[[[5], [0], [2]], [[4], [0], [6]], [[2], [1], [5]]]\n",
            "rvRichelotData[T] = source=E, determinant=[5], f=[[6], [3], [4], [5], [2], [6], [2]], factors=[[[4], [2], [3]], [[1], [0], [5]], [[5], [0], [2]]]\n",
            "rvRichelotData[Q] = source=H, determinant=[3], f=[[1], [4], [4], [5], [0], [5], [2]], factors=[[[6], [0], [1]], [[6], [5], [5]], [[1], [2], [6]]]\n",
            "rvHWWWeil[H] = [1, 2, 6, 14, 49]\n",
            "rvHWWWeil[Q] = [1, 2, 6, 14, 49]\n",
        )
    );
}

#[test]
fn extension_field_determinant_and_target_preserve_coordinates() {
    assert_eq!(
        stdout(interpret(
            r"
            rvunram U[3,1]{1,0,1};
            let z=rvUnramElem[U]{};
            let o=rvUnramElem[U]{1};
            let a=rvUnramElem[U]{0,1};
            let b=rvUnramElem[U]{1,1};
            rvhyperelliptic C[U]{z,a,o,b,o,o};
            rvrichelot D[C]{[],[1]}{[1],[],[1]}{[0,-2],[1],[1]};
            rvrichelotdual E[D];
            rvRichelotData[D];
            rvRichelotData[E];
            rvHWWPoints[C,1];
            rvHWWPoints[D,1];
            rvHWWWeil[C];
            rvHWWWeil[D];
            "
        )),
        concat!(
            "rvRichelotData[D] = source=C, determinant=[2, 1], f=[[2, 1], [0, 1], [1, 0], [2, 2], [1, 0], [1, 0], [2, 2]], factors=[[[2, 2], [2, 0], [1, 1]], [[0, 2], [0, 0], [1, 0]], [[1, 0], [0, 0], [2, 0]]]\n",
            "rvRichelotData[E] = source=D, determinant=[2, 1], f=[[0, 0], [0, 1], [1, 0], [1, 1], [1, 0], [1, 0]], factors=[[[0, 0], [1, 0]], [[1, 0], [0, 0], [1, 0]], [[0, 1], [1, 0], [1, 0]]]\n",
            "rvHWWPoints[C, 1] = 8, trace=2\n",
            "rvHWWPoints[D, 1] = 8, trace=2\n",
            "rvHWWWeil[C] = [1, -2, 2, -18, 81]\n",
            "rvHWWWeil[D] = [1, -2, 2, -18, 81]\n",
        )
    );
}

#[test]
fn degenerate_splittings_and_invalid_or_conflicting_inputs_are_diagnosed() {
    for (input, expected) in [
        (
            "rvhyperelliptic H[F]{s,z,h,z,s,z,o}; rvrichelot D[H]{[1],[],[1]}{[2],[],[1]}{[3],[],[1]};".to_string(),
            "product of elliptic curves over the algebraic closure",
        ),
        (
            "rvrichelot D[C]{[],[1]}{[1],[],[1]}{[5],[2],[1]};".to_string(),
            "must multiply to the source equation",
        ),
        (
            "rvrichelot D[C]{[1]}{[1],[],[1]}{[4],[2],[1]};".to_string(),
            "must have degree one or two",
        ),
        (
            "rvrichelot D[C]{[],[1],[1],[1]}{[1],[],[1]}{[4],[2],[1]};".to_string(),
            "must have degree one or two",
        ),
        (
            "rvrichelot D[C]{[],[1,1]}{[1],[],[1]}{[4],[2],[1]};".to_string(),
            "has at most 1 coordinates",
        ),
        (
            format!("rvrichelot C[C]{SPLITTING};"),
            "hyperelliptic curve is already defined",
        ),
        (
            format!("rvhyperelliptic I[7]{{0,4,2,5,2,1}}; rvrichelot D[I]{SPLITTING};"),
            "unknown unramified hyperelliptic curve 'I'",
        ),
        (
            format!("rvhyperelliptic I[F]{{o,o,z,o}}; rvrichelot D[I]{SPLITTING};"),
            "smooth genus-two curve",
        ),
        (
            format!("rvunram W[7,2]{{0,1}}; let w=rvUnramElem[W]{{1}}; rvhyperelliptic I[W]{{w,w,w,w,w,w}}; rvrichelot D[I]{SPLITTING};"),
            "named precision-one field",
        ),
        ("rvrichelotdual D[C];".to_string(), "has no stored Richelot splitting"),
    ] {
        let problem = diagnostic(interpret(&format!("{F7}{input}")));
        assert!(problem.contains(expected), "{problem}");
    }
}
