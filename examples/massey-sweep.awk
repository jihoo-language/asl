# Emit an ASL program for a range of arities. Use first=last for a single run.
# awk -v m=64 -v first=64 -v last=65 -f examples/massey-sweep.awk
BEGIN {
    if (m == "") m = 8
    if (first == "") first = 2
    if (last == "") last = m
    if (m !~ /^[0-9]+$/ || first !~ /^[0-9]+$/ || last !~ /^[0-9]+$/ ||
            m + 0 < 2 || first + 0 < 2 || last + 0 < first + 0) {
        print "use integer m >= 2 and 2 <= first <= last" > "/dev/stderr"
        exit 2
    }
    for (k = 1; k < last; k++) {
        printf "let b%d = binom[%s, %d];\n", k, m, k
        args = args (k == 1 ? "" : ", ") "b" k
        if (k + 1 >= first) printf "show powerdef(%s);\n", args
    }
}
