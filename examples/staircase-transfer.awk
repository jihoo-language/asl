# Emit the staircase quotient with a nonzero transferred m_n.
BEGIN {
    if (n == "") n = 5
    if (n !~ /^[0-9]+$/ || n < 3 || n > 31) {
        print "n must be an integer between 3 and 31 (primitive labels are 1 through 2n)" > "/dev/stderr"
        exit 2
    }
    if (bar == "") bar = 0
    if (bar !~ /^[01]$/) {
        print "bar must be 0 or 1" > "/dev/stderr"
        exit 2
    }
    if (local == "") local = 0
    if (local !~ /^[01]$/) {
        print "local must be 0 (total grading) or 1 (multigrading)" > "/dev/stderr"
        exit 2
    }
    if (bar) local = 1
    constructor = local ? "rvmkterm" : "rvkterm"
    for (k = 1; k <= 2*n; k++) selection = append(selection, k)

    emit("rvideal Stair[[" selection "]] {")
    for (j = 1; j <= n; j++)
        for (k = 1; k <= j; k++)
            if (j != n || k != 1)
                emit("    " monomial(j, n+k) (j == n && k == n ? "" : ","))
    emit("};")

    for (j = 1; j <= n; j++)
        emit("let a" j " = rvTorClass(" constructor "[Stair, [" (n+j) "]](" monomial(j, 0) "));")

    # Test every proper contiguous interval on the chosen adjacent classes.
    for (span = 2; span < n; span++)
        for (first = 1; first + span - 1 <= n; first++) {
            inputs = arguments(first, span)
            grade = local ? ("RVMTor(Stair, " (2*span-2) ", [" multidegree(first, span) "])")                           : ("RVTor(Stair, " (2*span-2) ", " (2*span) ")")
            emit("assert rvTHigher(" inputs ") == zero[" grade "];")
        }

    for (k = 1; k <= 2*n; k++)
        if (k != n && k != n+1)
            wedge = append(wedge, k)
    emit("let expected = rvTorClass(" constructor "[Stair, [" wedge "]](" monomial(n, n+1) "));")
    emit("let result = rvTHigher(" arguments(1, n) ");")
    emit("assert result == expected;")
    if (local)
        emit(sprintf("assert result == rvmtorgen[Stair, %d, [%s], 0];", 2*n-2, multidegree(1, n)))
    else
        emit(sprintf("assert result == rvtorgen[Stair, %d, %d, 0];", 2*n-2, 2*n))
    emit("show result;")

    if (bar) {
        emit("rvbarpath StairBar[Stair] {")
        for (j = 1; j <= n; j++)
            emit("    [" multidegree(j, 1) "]" (j == n ? "" : ","))
        emit("};")
        for (j = 0; j <= n; j++) cuts = append(cuts, j)
        emit("let word = rvBarWord[StairBar, [" cuts "]](" arguments(1, n) ");")
        emit("let top = rvBarWord[StairBar, [0, " n "]](result);")
        emit("assert rvBarD(word) == top;")
        emit(sprintf("assert rvBarD(top) == zero[RVBar(StairBar, 0, %d, %d)];", n, 2*n-2))
        emit(sprintf("rvbardim[StairBar, 0, %d, %d, %d];", n, 2*n, n))
        emit(sprintf("rvbardim[StairBar, 0, %d, %d, 1];", n, 2*n-1))
        emit("show top;")
    }

    for (i = 1; i <= line_count; i++) print lines[i]
}

function append(list, value) {
    return list (list == "" ? "" : ", ") value
}

function monomial(i, j,    k, powers) {
    for (k = 1; k <= 2*n; k++) powers = append(powers, (k == i || k == j) ? 1 : 0)
    return "rvmon[[" selection "]]{" powers "}"
}

function multidegree(first, span,    k, j, mu) {
    for (k = 1; k <= 2*n; k++) {
        j = k > n ? k-n : k
        mu = append(mu, (first <= j && j < first+span) ? 1 : 0)
    }
    return mu
}

function arguments(first, span,    k, args) {
    for (k = first; k < first + span; k++) args = append(args, "a" k)
    return args
}

function emit(line) {
    source_bytes += length(line) + 1
    if (source_bytes > 512 * 1024) {
        print "generated program exceeds the ASL source limit of 512 KiB" > "/dev/stderr"
        exit 2
    }
    lines[++line_count] = line
}
