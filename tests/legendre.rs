mod common;

use common::{diagnostic, example, interpret, stdout};

#[test]
fn actual_dwork_sums_agree_with_geometric_trace_and_hensel_at_high_precision() {
    assert_eq!(
        stdout(example("legendre-dwork")),
        concat!(
            "legendreTrace[7, 3] = 4, points=4\n",
            "legendreDwork[7, 3, 2] mod 49 = z=19, B_2=16, B_1=31, ratio=10, unit=39, trace=4\n",
            "dworkPeriod[7, 19, 1, 2] mod 49 = 31\n",
            "dworkPeriod[7, 19, 2, 2] mod 49 = 16\n",
            "legendreHensel[7, 3, 2] mod 49 = trace=4, unit=39\n",
            "legendreDwork[7, 3, 6] mod 117649 = z=82682, B_6=109726, B_5=14530, ratio=69933, unit=47716, trace=4\n",
            "legendreDwork[7, 5, 2] mod 49 = z=48, B_2=44, B_1=25, ratio=39, unit=10, trace=-4\n",
            "legendreDwork[5, 2, 2] mod 25 = z=7, B_2=9, B_1=18, ratio=13, unit=13, trace=-2\n",
            "legendreDwork[5, 2, 8] mod 390625 = z=280182, B_8=135706, B_7=5837, ratio=169738, unit=169738, trace=-2\n",
            "legendreDwork[7, 2, 2] mod 49 = z=1, B_2=42, B_1=0, supersingular, trace=0\n",
            "dworkCongruence[3, 3] mod 27 = yes, coefficients=105\n",
            "dworkCongruence[5, 2, 2] mod 25 = left=1, right=1\n",
        )
    );
}

#[test]
fn precision_one_and_supersingular_inputs_follow_the_unit_condition() {
    assert_eq!(
        stdout(interpret(
            "legendreDwork[7,3,1]; legendreDwork[7,2,1]; legendreTrace[7,2]; legendreHensel[7,2,4];"
        )),
        concat!(
            "legendreDwork[7, 3, 1] mod 7 = z=5, B_1=3, B_0=1, ratio=3, unit=4, trace=4\n",
            "legendreDwork[7, 2, 1] mod 7 = z=1, B_1=0, B_0=1, supersingular, trace=0\n",
            "legendreTrace[7, 2] = 0, points=8\n",
            "legendreHensel[7, 2, 4] mod 2401 = trace=0, supersingular\n",
        )
    );
}

#[test]
fn raw_periods_keep_units_after_vanishing_coefficients_and_support_wide_moduli() {
    assert_eq!(
        stdout(interpret(
            "dworkPeriod[3,2,2,2]; dworkPeriod[3,-170141183460469231731687303715884105728,2,40]; dworkPeriod[3,18446744073709551615,2,40]; dworkPeriod[11,10,0,4]; legendreHensel[5,2,27];"
        )),
        concat!(
            "dworkPeriod[3, 2, 2, 2] mod 9 = 6\n",
            "dworkPeriod[3, -170141183460469231731687303715884105728, 2, 40] mod 12157665459056928801 = 10172903787935937439\n",
            "dworkPeriod[3, 18446744073709551615, 2, 40] mod 12157665459056928801 = 1551867429027496588\n",
            "dworkPeriod[11, 10, 0, 4] mod 14641 = 1\n",
            "legendreHensel[5, 2, 27] mod 7450580596923828125 = trace=-2, unit=5881564383791575988\n",
        )
    );
}

#[test]
fn congruences_compare_full_coefficients_and_non_teichmueller_evaluations() {
    assert_eq!(
        stdout(interpret(
            "dworkCongruence[5,2]; dworkCongruence[11,2]; dworkCongruence[11,2,23]; dworkCongruence[5,2,-3];"
        )),
        concat!(
            "dworkCongruence[5, 2] mod 25 = yes, coefficients=145\n",
            "dworkCongruence[11, 2] mod 121 = yes, coefficients=1441\n",
            "dworkCongruence[11, 2, 23] mod 121 = left=36, right=36\n",
            "dworkCongruence[5, 2, -3] mod 25 = left=21, right=21\n",
        )
    );
}

#[test]
fn fields_precisions_and_resource_limits_report_precise_errors() {
    for (source, expected) in [
        ("legendreTrace[2,1];", "require an odd prime"),
        ("legendreTrace[9,2];", "9 is not prime"),
        ("legendreTrace[5,0];", "0 and 1 are singular"),
        ("legendreDwork[5,1,3];", "0 and 1 are singular"),
        ("legendreTrace[5,6];", "represented by an integer in 0..p"),
        ("legendreHensel[3,2,0];", "precision must be positive"),
        ("dworkPeriod[3,1,2,41];", "modulus is out of range"),
        (
            "dworkPeriod[3,1,41,2];",
            "truncation term count is out of range",
        ),
        (
            "legendreDwork[5,2,16];",
            "statement exceeds the computation limit",
        ),
        ("dworkCongruence[7,6];", "polynomial storage limit"),
    ] {
        let result = diagnostic(interpret(source));
        assert!(result.contains(expected), "{source} returned {result}");
    }
}
