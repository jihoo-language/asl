mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn actual_billionth_differential_and_its_torsion_witness_fit_two_cells() {
    assert_eq!(
        stdout(example("mixed-cone")),
        concat!(
            "MixedE(Cofiber, 1000000000, 0, 1000000000){g0}\n",
            "mixdim[Cofiber, 1000000001, 1, 0] = 0\n",
            "mixdim[Cofiber, 1000000001, 0, 1000000000] = 0\n",
            "mixbars[Cofiber] = free [], torsion [(2000000000, 1000000000)]\n",
            "mixbars[Split] = free [1, 2000000000], torsion []\n",
        )
    );
}

const STAIRCASE: &str = r"
    mixedcomplex C[3] {x:0} {} {};
    mixedcomplex E[3] {x0:0,x1:2,x2:4,y1:1,y2:3,y3:5,t:1,z:0}
        {x1 => [(y1,-1)],x2 => [(y2,-1)],t => [(z,1)]}
        {x0 => [(y1,1)],x1 => [(y2,1)]};
    mixhom F[C,E];
    let x = mixchain[C,0][(x,1)];
    let sx0 = mixU[0](mixunit[F,x,x0]);
    let f = add(sx0,add(mixU[1](mixunit[F,x,x1]),mixU[2](mixunit[F,x,x2])));
    let h = mixU[0](mixunit[F,x,t]);
    let g = add(f,mixD(h));
    mixcone K(f);
    mixcone L(g);
";

#[test]
fn higher_relation_shear_and_hom_are_actual_signed_chain_operations() {
    let program = format!(
        r"{STAIRCASE}
        let sx = mixConeLift[K](x);
        let sxu = mixU[0](sx);
        let t = mixchain[E,1][(t,1)];
        let k = add(sxu,mixU[0](mixConeI[K](t)));
        let y2 = mixchain[E,3][(y2,1)];
        let d11 = mixDelta[1](mixDelta[1](sx));
        assert d11 == mixConeI[K](y2);
        assert add(d11,mixDelta[0](mixDelta[2](sx))) == zero[MixedChain(K,3)];
        assert mixD(mixD(k)) == zero[MixedU(K,-1)];
        assert mixD(sxu) == mixConeI[K](mixEval(f,mixU[0](x)));
        let image = mixConeShear[L](h,k);
        assert image == mixU[0](mixConeLift[L](x));
        assert mixConeShear[K](pscale[2](h),image) == k;
        assert mixD(image) == mixConeShear[L](h,mixD(k));
        assert mixConeP[L](image) == mixU[0](x);
        mixhom End[K,K];
        let m = mixunit[End,w8,w8];
        let id = mixId[End];
        assert mixD(mixU[0](id)) == zero[MixedU(End,-1)];
        assert mixDelta[2](m) == mixunit[End,w8,w2];
        mixhom Alias[K,K];
        let alias = mixComp[Alias](mixU[0](id),mixU[0](m));
        assert alias == mixU[0](mixunit[Alias,w8,w8]);
        assert mixDelta[2](mixunit[Alias,w8,w8]) == mixunit[Alias,w8,w2];
        assert mixEval(alias,sxu) == sxu;
        assert mixD(mixEval(mixU[0](m),sxu)) ==
            add(mixEval(mixD(mixU[0](m)),sxu),mixEval(mixU[0](m),mixD(sxu)));
        mixbars[K];
        mixbars[L];
        "
    );
    assert_eq!(
        stdout(interpret(&program)),
        "mixbars[K] = free [5], torsion []\nmixbars[L] = free [5], torsion []\n"
    );
}

#[test]
fn declared_higher_terms_check_homogeneity_and_square_zero() {
    let input = r"
        mixedhigher M[5] {x:-2,y:3,z:3} {x => [(3,y,-1),(3,y,2),(3,z,5)]};
        let x = mixchain[M,-2][(x,1)];
        let y = mixchain[M,3][(y,1)];
        assert mixDelta[3](x) == y;
        assert mixDelta[4](x) == zero[MixedChain(M,5)];
        assert mixD(mixU[7](x)) == mixU[10](y);
        assert mixD(mixD(mixU[7](x))) == zero[MixedU(M,-18)];
        assert mixD(mixgen[M,3,-2,0,0]) == mixgen[M,3,-3,3,0];
        mixedhigher N[5] {a:0,b:1,c:3}
            {a => [(1,b,1),(2,c,1),(2,c,-1)]};
        let a = mixchain[N,0][(a,1)];
        assert mixB(a) == mixDelta[1](a);
        assert mixDelta[2](a) == zero[MixedChain(N,3)];
        mixbars[M];
    ";
    assert_eq!(
        stdout(interpret(input)),
        "mixbars[M] = free [3], torsion [(3, 3)]\n"
    );
    for (program, expected) in [
        (
            "mixedhigher K[3] {x:0,y:2} {x => [(2,y,1)]};",
            "wrong homological degree",
        ),
        (
            "mixedhigher K[3] {x:0,y:3,z:6} {x => [(2,y,1)],y => [(2,z,1)]};",
            "higher mixed D squared is nonzero",
        ),
        (
            "mixedhigher K[3] {x:0} {x => [],x => []};",
            "duplicate higher mixed differential column",
        ),
        (
            "mixedhigher K[3] {x:0,y:0} {x => [(0,y,0)]};",
            "wrong homological degree",
        ),
        ("mixedhigher K[4] {x:0} {};", "requires a prime field"),
    ] {
        let error = diagnostic(interpret(program));
        assert!(error.contains(expected), "{program}: {error}");
    }
}

#[test]
fn suspension_uses_the_negative_of_the_source_differential() {
    let input = r"
        mixedcomplex C[7] {x:0,t:1} {t => [(x,2)]} {};
        mixedcomplex E[7] {y:4} {} {};
        mixhom F[C,E];
        mixcone K(zero[MixedU(F,0)]);
        let x = mixU[0](mixchain[C,0][(x,1)]);
        let t = mixU[0](mixchain[C,1][(t,1)]);
        let y = mixU[0](mixchain[E,4][(y,1)]);
        let st = mixConeLift[K](t);
        assert mixD(st) == pscale[5](mixConeLift[K](x));
        assert mixConeP[K](mixD(st)) == pscale[6](mixD(mixConeP[K](st)));
        assert mixConeP[K](mixConeI[K](y)) == zero[MixedU(C,3)];
        assert mixD(mixConeI[K](y)) == mixConeI[K](mixD(y));
        mixbars[K];
    ";
    assert_eq!(
        stdout(interpret(input)),
        "mixbars[K] = free [4], torsion []\n"
    );
}

#[test]
fn ill_typed_operations_and_nonclosed_cone_maps_have_precise_diagnostics() {
    for (statement, expected) in [
        ("mixcone Z(sx0);", "requires a closed degree-zero map"),
        ("mixcone Z(h);", "requires a degree-zero MixedU map"),
        ("mixcone Z(zero[MixedU(C,0)]);", "defined by mixhom"),
        ("show mixConeI[K](x);", "wrong endpoint"),
        ("show mixConeLift[E](x);", "requires a mixed mapping cone"),
        (
            "show mixConeShear[L](zero[MixedU(F,1)],mixU[0](mixConeLift[K](x)));",
            "must differ by the differential",
        ),
        (
            "mixhom G[E,C]; show mixConeShear[L](zero[MixedU(G,1)],mixU[0](mixConeLift[K](x)));",
            "incompatible endpoints",
        ),
        ("show mixb(mixConeLift[K](x));", "require a strict mixed complex"),
        ("show mixB(mixConeLift[K](x));", "require a strict mixed complex"),
        ("mixminimal N[K];", "requires a strict mixed complex"),
        ("mixtensor N[K,L];", "requires two strict mixed complexes"),
        (
            "mixedcomplex W[3] {z:9223372036854775807} {} {}; mixedcomplex R[3] {} {} {}; mixhom G[W,R]; mixcone Z(zero[MixedU(G,0)]);",
            "suspension degree in mixed mapping cone is out of range",
        ),
    ] {
        let result = diagnostic(interpret(&format!("{STAIRCASE} {statement}")));
        assert!(result.contains(expected), "{statement}: {result}");
    }
}
