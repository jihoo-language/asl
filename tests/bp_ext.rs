mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn corrected_height_two_class_realizes_order_sixteen() {
    assert_eq!(
        stdout(example("bp-ext")),
        concat!(
            "bpext[2, 1, 1, 8] = free 0, torsion [8]\n",
            "bporder[8] = 8\n",
            "bpext[2, 2, 1, 8] = free 0, torsion [16]\n",
            "bporder[8] = 16\n",
            "bpext[3, 2, 1, 12] = free 0, torsion [9]\n",
            "bpext[3, 2, 2, 12] = free 0, torsion [3]\n",
            "bporder[12] = 3\n",
        )
    );
}

#[test]
fn normalized_powers_and_exact_cocycles_have_the_expected_orders() {
    let source = r"
        let t = bpt[2,2,2,1,1];
        let u = bpt[2,2,2,2,1];
        let eta_square = bpMul(t,u);
        let twice = bpMul(bpint[2,2,2,2],eta_square);
        bpext[2,2,2,4];
        bporder[4](eta_square);
        bporder[4](twice);
        bporder[4](bpMul(bpint[2,2,2,3],eta_square));
        assert bpD(bpNull[4](eta_square)) == twice;
        assert bpD(bpNull[4](twice)) == twice;
        assert bpD(bpNull[4](bpMul(bpint[2,2,2,3],eta_square))) == bpMul(bpint[2,2,2,3],twice);
        let alpha = bpt[5,2,1,1,1];
        bporder[8](alpha);
        assert bpD(bpNull[8](alpha)) == bpMul(bpint[5,2,1,5],alpha);
    ";
    assert_eq!(
        stdout(interpret(source)),
        concat!(
            "bpext[2, 2, 2, 4] = free 0, torsion [2]\n",
            "bporder[4] = 2\n",
            "bporder[4] = 1\n",
            "bporder[4] = 2\n",
            "bporder[8] = 5\n",
        )
    );
}

#[test]
fn zero_and_empty_grades_bypass_unnecessary_bases() {
    let max = usize::MAX;
    let source = format!(
        "bpext[2,0,0,0]; bpext[2,0,3,0]; bpext[2,0,0,{max}];\n\
         bpext[2,1,0,0]; bpext[3,2,3,8]; bpext[2,1,16,32];\n\
         bporder[0](bpint[2,0,0,-5]); bporder[{max}](zero[BP(3,2,2)]);\n\
         assert bpNull[{max}](zero[BP(3,2,2)]) == zero[BP(3,2,1)];"
    );
    assert_eq!(
        stdout(interpret(&source)),
        format!(
            "bpext[2, 0, 0, 0] = free 1, torsion []\n\
             bpext[2, 0, 3, 0] = free 0, torsion []\n\
             bpext[2, 0, 0, {max}] = free 0, torsion []\n\
             bpext[2, 1, 0, 0] = free 1, torsion []\n\
             bpext[3, 2, 3, 8] = free 0, torsion []\n\
             bpext[2, 1, 16, 32] = free 0, torsion [2]\n\
             bporder[0] = infinite\n\
             bporder[{max}] = 1\n"
        )
    );
}

#[test]
fn integral_cobar_diagnostics_identify_which_condition_failed() {
    for (source, expected) in [
        ("bpext[4,2,1,8];", "BP requires a prime"),
        ("bpext[2,9,1,8];", "BP height exceeds the limit"),
        ("bpext[2,1,17,34];", "BP arrow arity exceeds the limit"),
        (
            "bporder[0](bpReduce[1](bpint[2,2,0,1]));",
            "bporder requires an integral BP cochain",
        ),
        (
            "show bpNull[2](bpReduce[1](bpt[2,2,1,1,1]));",
            "bpNull requires an integral BP cochain",
        ),
        (
            "bporder[2](bpv[2,1,1,1]);",
            "normalized cochain using every arrow slot",
        ),
        (
            "bporder[4](bpt[2,2,1,1,1]);",
            "does not have the requested homogeneous internal degree",
        ),
        (
            "show bpNull[6](bpt[2,2,1,1,2]);",
            "requires a closed cochain",
        ),
        (
            "show bpNull[0](bpint[2,2,0,1]);",
            "bpNull requires positive BP cobar arity",
        ),
    ] {
        let error = diagnostic(interpret(source));
        assert!(error.contains(expected), "{source}\n{error}");
    }
    let error = diagnostic(interpret(&format!("bpext[2,1,1,{}];", usize::MAX - 1)));
    assert!(error.contains("basis exceeds the storage limit"), "{error}");
}

#[test]
fn weight_sixteen_uses_certified_modular_ranks_and_witnesses() {
    let source = r"
        bpext[2,2,1,16];
        let q = bpv[2,2,2,1];
        let x = bpt[2,2,2,1,1];
        let y = bpt[2,2,2,2,2];
        let boundary = bpD(bpMul(bpPow[4](q),bpMul(x,y)));
        bporder[16](boundary);
        assert bpD(bpNull[16](boundary)) == boundary;
    ";
    assert_eq!(
        stdout(interpret(source)),
        "bpext[2, 2, 1, 16] = free 0, torsion [32]\nbporder[16] = 1\n"
    );
}

#[test]
fn connecting_cocycles_feed_integral_orders_and_lifts() {
    let source = r"
        let y = bpv[3,2,0,2];
        let first = bpDelta(bpReduce[2](y));
        let beta = bpDelta(first);
        bporder[12](beta);
        assert bpD(bpNull[12](beta)) == bpMul(bpint[3,2,2,3],beta);
        let fourth = bpPow[4](bpv[2,2,0,1]);
        let v1v2 = bpMul(bpv[2,2,0,1],bpv[2,2,0,2]);
        let primitive = add(fourth,bpMul(bpint[2,2,0,8],v1v2));
        let corrected = bpBockP[4](primitive);
        bporder[8](corrected);
        assert bpNull[8](corrected) == primitive;
    ";
    assert_eq!(
        stdout(interpret(source)),
        "bporder[12] = 3\nbporder[8] = 16\n"
    );
}
