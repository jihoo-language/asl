mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn a_homologically_zero_map_has_a_four_cell_nonsplit_cone() {
    assert_eq!(
        stdout(example("mixed-triangle")),
        concat!(
            "MixedH(F, 0){g0}\n",
            "mixHdim[K, 34] = 2\n",
            "mixHdim[K, 0] = 1\n",
            "mixbars[K] = free [], torsion [(34, 17), (200000000, 199999983)]\n",
            "mixbars[Split] = free [], torsion [(34, 100000000), (200000000, 100000000)]\n",
        ),
    );
}

const RANK_ONE: &str = r"
    mixedcomplex X[5] {x0:0,x1:0} {} {};
    mixedcomplex Y[5] {y0:0,y1:0} {} {};
    mixhom F[X,Y];
    let f = pscale[2](mixU[0](mixunit[F,x0,y0]));
    mixcone K(f);
    let x0 = mixHgen[X,0,0];
    let x1 = mixHgen[X,0,1];
    let y0 = mixHgen[Y,0,0];
    let y1 = mixHgen[Y,0,1];
";

#[test]
fn the_actual_triangle_maps_have_nonzero_images_and_exact_middle() {
    let source = format!(
        r"{RANK_ONE}
        let z0 = mixHgen[K,0,0];
        let z1 = mixHgen[K,1,0];
        assert mixConeConnect[K](x0) == pscale[2](y0);
        assert mixConeConnect[K](x1) == zero[MixedH(Y,0)];
        assert mixConeConnect[K](x0) == mixInduce(f,x0);
        assert mixConeP[K](z1) == x1;
        assert mixConeI[K](y0) == zero[MixedH(K,0)];
        assert mixConeI[K](y1) == z0;
        assert mixConeP[K](z0) == zero[MixedH(X,-1)];
        assert mixConeP[K](mixU[3](z1)) == mixU[3](x1);
        assert mixConeConnect[K](mixU[3](x0)) == pscale[2](mixU[3](y0));
        assert mixConeI[K](mixU[3](y1)) == mixU[3](z0);
        let exact = mixConeI[K](mixHRep(y0));
        let witness = mixHPrimitive(exact);
        assert witness == pscale[3](mixConeLift[K](mixHRep(x0)));
        assert mixD(witness) == exact;
        let id = (\z : MixedH(X,0) => z);
        let composite = (\z : MixedH(X,0) => mixHClass(mixHRep(z)));
        assert id == composite;
        mixbars[K];
        "
    );
    assert_eq!(
        stdout(interpret(&source)),
        "mixbars[K] = free [0, 1], torsion []\n"
    );
}

#[test]
fn homogeneous_induction_respects_composition_and_polynomial_degree() {
    let source = format!(
        r"{RANK_ONE}
        mixedcomplex E[5] {{e:3}} {{}} {{}};
        mixhom G[Y,E];
        mixhom L[X,E];
        let g = mixU[1](mixunit[G,y0,e]);
        let gf = mixComp[L](g,f);
        let gfy = mixInduce(gf,x0);
        assert gfy == mixInduce(g,mixInduce(f,x0));
        assert gfy == pscale[2](mixHgen[E,1,0]);
        assert mixInduce(gf,mixU[10](x0)) == mixU[10](gfy);
        assert mixInduce(gf,x1) == zero[MixedH(E,1)];
        "
    );
    assert!(stdout(interpret(&source)).is_empty());
}

#[test]
fn a_boundary_can_mix_multiple_incoming_pivots_and_omit_killed_bars() {
    let source = r"
        mixedhigher C[7] {x:4,y:2,z:0,m:1,n:-1}
            {m=>[(2,x,2),(1,y,-1)],n=>[(3,x,1),(1,z,3)]};
        let mx = mixU[1](mixchain[C,1][(m,1)]);
        let n0 = mixU[0](mixchain[C,-1][(n,1)]);
        let composite = add(pscale[4](mx),pscale[2](n0));
        let a = mixD(composite);
        let w = mixHPrimitive(a);
        assert w == composite;
        assert mixD(w) == a;
        assert mixHClass(a) == zero[MixedH(C,-2)];
        let z = mixU[0](mixchain[C,0][(z,1)]);
        let h = mixHClass(z);
        let zAtMinusTwo = mixU[1](z);
        assert mixHClass(add(zAtMinusTwo,a)) == mixHClass(zAtMinusTwo);
        assert mixHClass(mixHRep(h)) == h;
        let u = mixU[8](h);
        assert u == mixHClass(mixU[8](z));
        mixHdim[C,4];
        mixHdim[C,2];
        mixHdim[C,-10];
    ";
    assert_eq!(
        stdout(interpret(source)),
        "mixHdim[C, 4] = 1\nmixHdim[C, 2] = 2\nmixHdim[C, -10] = 1\n"
    );
}

#[test]
fn invalid_homology_operations_fail_at_the_operation() {
    for (statement, expected) in [
        ("show mixHClass(mixConeLift[K](mixHRep(x0)));", "requires a total cycle"),
        ("show mixHPrimitive(mixHRep(x0));", "requires an exact total chain"),
        ("show mixHPrimitive(mixConeLift[K](mixHRep(x0)));", "requires an exact total chain"),
        ("show mixConeConnect[K](mixHRep(x0));", "a MixedH class"),
        ("show mixConeLift[K](x0);", "a MixedChain or MixedU chain"),
        ("show mixInduce(f,y0);", "wrong source complex"),
        ("show mixHgen[X,0,2];", "generator index"),
        ("show mixHPrimitive(zero[MixedU(X,9223372036854775807)]);", "degree is out of range"),
        (
            "mixedcomplex Z[5] {z:0,w:-1} {z=>[(w,1)]} {}; mixhom G[X,Z]; show mixInduce(mixU[0](mixunit[G,x0,z]),zero[MixedH(X,0)]);",
            "requires a closed total Hom map",
        ),
    ] {
        let source = format!("{RANK_ONE} {statement}");
        let error = diagnostic(interpret(&source));
        assert!(error.contains(expected), "{statement}: {error}");
    }
}

#[test]
fn a_free_bar_shifts_across_the_full_signed_degree_range() {
    let source = r"
        mixedcomplex C[3] {x:9223372036854775806} {} {};
        let h = mixHgen[C,9223372036854775806,0];
        let shifted = mixU[9223372036854775807](h);
        assert shifted == mixHgen[C,-9223372036854775808,0];
        assert mixHClass(mixHRep(shifted)) == shifted;
        show mixHRep(shifted);
    ";
    assert_eq!(
        stdout(interpret(source)),
        "MixedU(C, -9223372036854775808){u^9223372036854775807*x}\n"
    );
}
