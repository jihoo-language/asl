# Generate two disjoint arms with unequal lengths. Run with -v left=3 -v right=5 -v prime=5.
# Set -v explicit=1 to verify every prescribed interval witness against the solver.
function name(arm, i, len) {
    return i == 0 ? "bottom" : (i == len ? "top" : arm i)
}
function path(arm, len,    i) {
    printf "%s", name(arm, 0, len)
    for (i = 1; i <= len; ++i)
        printf "<%s", name(arm, i, len)
}
function inputs(arm, len,    i) {
    for (i = 1; i <= len; ++i)
        printf "%s%s%d", (i > 1 ? "," : ""), arm, i
}
function branches(arm, len,    i, j, a, b, sep, factor) {
    for (i = 1; i <= len; ++i) {
        a = name(arm, i - 1, len)
        b = name(arm, i, len)
        printf "let %s%d = shExtGen[F%s,F%s,1,0];\n", arm, i, a, b
    }
    printf "let result%s = shMasseySolve(", arm
    inputs(arm, len)
    print ");"
    if (!explicit)
        return
    printf "assert result%s == shMassey(", arm
    inputs(arm, len)
    printf "){"
    sep = ""
    for (i = 0; i < len - 1; ++i)
        for (j = i + 2; j <= len; ++j) {
            if (i == 0 && j == len)
                continue
            a = name(arm, i, len)
            b = name(arm, j, len)
            factor = (j - i) % 2 ? 1 : -1
            printf "%s%d<%d=>shHomTerm[F%s,F%s,[%s,%s],0,0,%d]", sep, i, j, a, b, a, b, factor
            sep = ","
        }
    print "};"
}
BEGIN {
    if (left == "")
        left = 3
    if (right == "")
        right = 5
    if (prime == "")
        prime = 5
    if (left != int(left) || right != int(right) || left < 2 || right < 2 || left > 128 || right > 128 ||
            left == right || prime != int(prime) || prime < 2 || prime > 4294967295 ||
            (explicit != "" && explicit != 0 && explicit != 1)) {
        print "Use distinct integers left,right in 2..128, a u32 prime, and optionally explicit=1." > "/dev/stderr"
        exit 1
    }
    printf "space Arms {"
    path("a", left)
    printf ","
    path("b", right)
    print "};"
    points[0] = "bottom"
    for (i = 1; i < left; ++i)
        points[i] = "a" i
    for (i = 1; i < right; ++i)
        points[left+i-1] = "b" i
    count = left + right
    points[count-1] = "top"
    for (i = 0; i < count; ++i) {
        printf "sheaf F%s[Arms,%.0f] {", points[i], prime
        for (j = 0; j < count; ++j)
            printf "%s%s:%d", (j ? "," : ""), points[j], (i == j)
        printf "} {"
        for (j = 1; j <= left; ++j)
            printf "%s%s<%s:[]", (j > 1 ? "," : ""), name("a", j-1, left), name("a", j, left)
        for (j = 1; j <= right; ++j)
            printf ",%s<%s:[]", name("b", j-1, right), name("b", j, right)
        print "};"
    }
    branches("a", left)
    branches("b", right)
    print "shExtDim[Fbottom,Ftop,2];"
    print "show resulta;"
    print "show resultb;"
    printf "assert resultb == shScale[%d](resulta);\n", ((left + right) % 2 ? 1 : -1)
}
