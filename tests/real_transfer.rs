mod common;

use common::{diagnostic, example, interpret, stdout};

const GENERATORS: &str = r"
    let rho = extClass(rmcword[1,0]{});
    let h0 = extClass(rmcword[0,0]{([],[1])});
    let h1 = extClass(rmcword[0,0]{([1],[])});
    let h02 = yoneda(h0,h0);
    let rhoh1 = yoneda(rho,h1);
    let unit = extClass(rmcword[0,0]{});
";

#[test]
fn real_transfer_example_has_nonzero_higher_operations() {
    assert_eq!(
        stdout(example("real-motivic-transfer")),
        concat!(
            "rmTripleContains = yes\n",
            "rhoorder = 4\n",
            "RMExt(1, 1, -2){g0}\n",
            "RMExt(2, 5, 0){g0}\n",
        )
    );
}

#[test]
fn minimality_strict_units_and_multilinearity_over_f2() {
    let source = r"
        assert ainfty(rho) == zero[RMExt(1,-1,-1)];
        assert alift(rho) == rmcword[1,0]{};
        assert ainfty(rho,rho) == extClass(rmcword[2,0]{});
        assert ainfty(rho,h0) == zero[RMExt(1,0,-1)];
        assert ainfty(unit,h0) == h0;
        assert ainfty(h0,unit) == h0;
        assert alift(unit,h0) == zero[RMCobar(0,1,0)];
        assert alift(h0,unit) == zero[RMCobar(0,1,0)];
        assert ainfty(unit,rho,h0) == zero[RMExt(0,0,-1)];
        assert ainfty(h0,unit,h1) == zero[RMExt(1,3,1)];
        assert alift(h0,unit,h1) == zero[RMCobar(0,3,1)];
        assert ainfty(rho,rho,h0) == zero[RMExt(0,-1,-2)];
        let mix = add(h0,rhoh1);
        assert ainfty(rho,h0,mix,h1,h02) == add(ainfty(rho,h0,h0,h1,h02),ainfty(rho,h0,rhoh1,h1,h02));
        assert alift(rho,h0,mix,h1,h02) == add(alift(rho,h0,h0,h1,h02),alift(rho,h0,rhoh1,h1,h02));
    ";
    assert!(stdout(interpret(&format!("{GENERATORS}{source}"))).is_empty());
}

type Input = (&'static str, i32, i32, i32);

fn cohomological_degree(inputs: &[Input]) -> i32 {
    inputs.iter().map(|x| x.1).sum::<i32>() + 2 - inputs.len() as i32
}

fn transfer_call<'a>(operation: &str, inputs: impl IntoIterator<Item = &'a str>) -> String {
    let names = inputs.into_iter().collect::<Vec<_>>().join(",");
    format!("{operation}({names})")
}

fn sum(terms: Vec<String>, zero: String) -> String {
    terms.into_iter().fold(zero, |a, b| format!("add({a},{b})"))
}

fn identities(inputs: &[Input]) -> String {
    let n = inputs.len();
    let s = cohomological_degree(inputs);
    let t = inputs.iter().map(|x| x.2).sum::<i32>();
    let w = inputs.iter().map(|x| x.3).sum::<i32>();
    let zero_class = format!("zero[RMExt({},{t},{w})]", s + 1);
    let zero_cochain = format!("zero[RMCobar({s},{t},{w})]");
    let mut class_terms = Vec::new();
    let mut chain_terms = Vec::new();

    for length in 2..=n {
        for i in 0..=n - length {
            let end = i + length;
            if cohomological_degree(&inputs[i..end]) < 0 {
                continue;
            }
            let inside = transfer_call("ainfty", inputs[i..end].iter().map(|x| x.0));
            let args = inputs[..i]
                .iter()
                .map(|x| x.0)
                .chain(std::iter::once(inside.as_str()))
                .chain(inputs[end..].iter().map(|x| x.0));
            if length < n {
                class_terms.push(transfer_call("ainfty", args.clone()));
            }
            chain_terms.push(transfer_call("alift", args));
        }
    }
    for i in 1..n {
        let (a, b) = inputs.split_at(i);
        if cohomological_degree(a) < 1 || cohomological_degree(b) < 1 {
            continue;
        }
        let a = transfer_call("alift", a.iter().map(|x| x.0));
        let b = transfer_call("alift", b.iter().map(|x| x.0));
        chain_terms.push(format!("concat({a},{b})"));
    }
    let class = sum(class_terms, zero_class.clone());
    let chain = sum(chain_terms, zero_cochain.clone());
    let differential = if s == 0 {
        zero_cochain
    } else {
        let lift = transfer_call("alift", inputs.iter().map(|x| x.0));
        format!("cobarD({lift})")
    };
    format!("assert {class} == {zero_class}; assert {chain} == {differential};")
}

#[test]
fn stasheff_and_morphism_identities_with_degree_zero_inputs() {
    let rho = ("rho", 0, -1, -1);
    let h0 = ("h0", 1, 1, 0);
    let h1 = ("h1", 1, 2, 1);
    let h02 = ("h02", 2, 2, 0);
    let rhoh1 = ("rhoh1", 1, 1, 0);
    for inputs in [
        [rho, rho, h0].as_slice(),
        &[rho, h0, rho, h02],
        &[rho, h0, h0, h1, h02],
        &[rho, h0, h02, rho, rhoh1, h02],
    ] {
        assert!(stdout(interpret(&format!("{GENERATORS}{}", identities(inputs)))).is_empty());
    }
}

#[test]
fn real_transfer_checks_final_degrees_operand_families_and_work_limits() {
    for (expression, expected) in [
        ("alift(rho,rho)", "negative cohomological degree"),
        ("ainfty(h0,extgen[1,1,0])", "requires RMExt classes"),
        ("ainfty(h0,extRep(h0))", "requires RMExt classes"),
        (
            "ainfty(rho,zero[RMExt(1,0,-9223372036854775808)])",
            "tridegree exceeds",
        ),
    ] {
        let error = diagnostic(interpret(&format!("{GENERATORS}show {expression};")));
        assert!(error.contains(expected), "{error}");
    }

    let inputs = std::iter::repeat_n("h0", 800).collect::<Vec<_>>().join(",");
    let error = diagnostic(interpret(&format!("{GENERATORS}show ainfty({inputs});")));
    assert!(
        error.contains("transfer lift table exceeds the storage limit"),
        "{error}"
    );
}
