# Emit the staircase quotient with a nonzero transferred m_n.
BEGIN {
    if (n == "") n = 5
    if (n !~ /^[0-9]+$/ || n < 3 || n > 31) {
        print "n must be an integer between 3 and 31 (primitive labels are 1 through 2n)" > "/dev/stderr"
        exit 2
    }
    for (k = 1; k <= 2*n; k++) selection = append(selection, k)

    emit("rvideal Stair[[" selection "]] {")
    for (j = 1; j <= n; j++)
        for (k = 1; k <= j; k++)
            if (j != n || k != 1)
                emit("    " monomial(j, n+k) (j == n && k == n ? "" : ","))
    emit("};")

    for (j = 1; j <= n; j++)
        emit("let a" j " = rvTorClass(rvkterm[Stair, [" (n+j) "]](" monomial(j, 0) "));")

    # Test every proper contiguous interval on the chosen adjacent classes.
    for (span = 2; span < n; span++)
        for (first = 1; first + span - 1 <= n; first++) {
            inputs = arguments(first, span)
            emit(sprintf("assert rvTHigher(%s) == zero[RVTor(Stair, %d, %d)];", inputs, 2*span-2, 2*span))
        }

    for (k = 1; k <= 2*n; k++)
        if (k != n && k != n+1)
            wedge = append(wedge, k)
    emit("let expected = rvTorClass(rvkterm[Stair, [" wedge "]](" monomial(n, n+1) "));")
    emit("let result = rvTHigher(" arguments(1, n) ");")
    emit("assert result == expected;")
    emit(sprintf("assert result == rvtorgen[Stair, %d, %d, 0];", 2*n-2, 2*n))
    emit("show result;")

    for (i = 1; i <= line_count; i++) print lines[i]
}

function append(list, value) {
    return list (list == "" ? "" : ", ") value
}

function monomial(i, j,    k, powers) {
    for (k = 1; k <= 2*n; k++) powers = append(powers, (k == i || k == j) ? 1 : 0)
    return "rvmon[[" selection "]]{" powers "}"
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
