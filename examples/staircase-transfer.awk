# Emit the staircase quotient with a nonzero transferred m_n.
BEGIN {
    if (n == "") n = 5
    if (n !~ /^[0-9]+$/ || n < 3 || n > 31) {
        print "n must be an integer between 3 and 31 (primitive labels are 1 through 2n)" > "/dev/stderr"
        exit 2
    }
    if (hh == "") hh = 0
    if (hh !~ /^[01]$/) {
        print "hh must be 0 or 1" > "/dev/stderr"
        exit 2
    }
    if (omega == "") omega = 0
    if (omega !~ /^[01]$/) {
        print "omega must be 0 or 1" > "/dev/stderr"
        exit 2
    }
    if (dual == "") dual = 0
    if (dual !~ /^[01]$/) {
        print "dual must be 0 or 1" > "/dev/stderr"
        exit 2
    }
    if (pages == "") pages = 0
    if (pages !~ /^[01]$/) {
        print "pages must be 0 or 1" > "/dev/stderr"
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
    if (pages || dual || omega) bar = 1
    if (bar || hh) local = 1
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

    if (pages) {
        for (r = 0; r <= n-1; r++) {
            emit(sprintf("let e%d = rvBarClass[%d,%d](word);", r,r,n))
            emit(sprintf("assert e%d == rvbaregen[StairBar,0,%d,%d,%d,%d,0];", r,n,r,2*n,n))
            if (r < n-1)
                emit(sprintf("assert rvBarD(e%d) == zero[RVBarE(StairBar,0,%d,%d,%d,%d)];", r,n,r,2*n-1,n-r))
            if (r > 0) {
                emit(sprintf("assert rvBarPrev(e%d) == e%d;", r,r-1))
                emit(sprintf("assert rvBarNext(e%d) == e%d;", r-1,r))
            }
        }
        emit(sprintf("let targetPage = rvBarClass[%d,1](top);", n-1))
        emit(sprintf("assert rvBarD(e%d) == targetPage;", n-1))
        emit(sprintf("assert targetPage == rvbaregen[StairBar,0,%d,%d,%d,1,0];", n,n-1,2*n-1))
        emit(sprintf("assert rvBarNext(targetPage) == zero[RVBarE(StairBar,0,%d,%d,%d,1)];", n,n,2*n-1))
        emit(sprintf("rvbaredim[StairBar,0,%d,%d,%d,%d];", n,n,2*n,n))
        emit(sprintf("rvbaredim[StairBar,0,%d,%d,%d,1];", n,n,2*n-1))
        emit(sprintf("show rvBarD(e%d);", n-1))
    }

    if (dual) {
        for (j = 1; j <= n; j++) {
            emit(sprintf("let dual%d = rvBarDualize(rvBarWord[StairBar,[%d,%d]](a%d));", j,j-1,j,j))
            emit(sprintf("let ext%d_%d = rvPathClass(dual%d);", j,j,j))
            emit(sprintf("assert ext%d_%d == rvpathextgen[StairBar,%d,%d,2,0];",j,j,j-1,j))
        }
        emit("let cup1 = dual1;")
        for (j = 2; j <= n; j++)
            emit(sprintf("let cup%d = rvPathCup(cup%d,dual%d);", j,j-1,j))
        for (span = 2; span <= n; span++) {
            for (first = 1; first+span-1 <= n; first++) {
                last = first+span-1
                emit(sprintf("let ext%d_%d = rvPathCup(ext%d_%d,ext%d_%d);",first,last,first,last-1,last,last))
                if (span < n)
                    emit(sprintf("assert ext%d_%d == rvpathextgen[StairBar,%d,%d,%d,0];",first,last,first-1,last,2*span))
            }
        }
        emit(sprintf("assert cup%d == rvBarDualize(word);",n))
        emit(sprintf("assert rvBarDualD(rvBarDualize(top)) == cup%d;",n))
        emit(sprintf("assert rvPathClass(cup%d) == ext1_%d;",n,n))
        emit(sprintf("assert ext1_%d == zero[RVPathExt(StairBar,0,%d,%d)];",n,n,2*n))
        emit(sprintf("show ext1_%d;",n))
    }

    if (omega) {
        emit("let rectified = rvBCDesuspend(word);")
        emit("let inclusion = rvBCDesuspend(top);")
        emit("let boundary0 = inclusion;")
        for (j = 1; j < n; j++) {
            left = "rvBCDesuspend(rvBarWord[StairBar,[" cuts_between(0,j) "]](" arguments(1,j) "))"
            right = "rvBCDesuspend(rvBarWord[StairBar,[" cuts_between(j,n) "]](" arguments(j+1,n-j) "))"
            emit(sprintf("let split%d = rvBCMul(%s,%s);",j,left,right))
            emit(sprintf("let boundary%d = add(boundary%d,split%d);",j,j-1,j))
        }
        emit(sprintf("assert rvBCD(rectified) == boundary%d;",n-1))
        emit(sprintf("assert rvBCD(rvBCD(rectified)) == zero[RVBarCobar(StairBar,0,%d,%d)];",n,2*n-3))
        emit("assert rvBCHomotopy(rvBCD(rectified)) == rectified;")
        emit(sprintf("assert rvBCHomotopy(rectified) == zero[RVBarCobar(StairBar,0,%d,%d)];",n,2*n))
        emit("assert rvBCProject(split1) == result;")
        emit("assert rvBCProject(inclusion) == result;")
        emit("assert add(rvBCD(rvBCHomotopy(split1)),rvBCHomotopy(rvBCD(split1))) == add(split1,inclusion);")
        emit(sprintf("assert rvBCHomotopy(rvBCHomotopy(split1)) == zero[RVBarCobar(StairBar,0,%d,%d)];",n,2*n))
        emit("show rvBCProject(split1);")
    }

    if (hh) {
        weight = multidegree(1,n)
        emit("let unit = rvmtorgen[Stair,0,[" multidegree(1,0) "],0];")
        emit("let circle = rvHHWord(unit," arguments(1,n) ");")
        emit("let cyclic_top = rvHHWord(unit,result);")
        emit("let two_ends = add(rvHHWord(" arguments(1,n) "),rvHHWord(a" n "," arguments(1,n-1) "));")
        emit("let hh_boundary = rvHHd(circle);")
        emit("assert hh_boundary == add(cyclic_top,two_ends);")
        emit("assert rvHHLength[1](hh_boundary) == cyclic_top;")
        if (n <= 9)
            emit(sprintf("assert rvHHd(hh_boundary) == zero[RVHHChain(Stair,%d,[%s])];",2*n-2,weight))
        emit(sprintf("assert rvConnesB(hh_boundary) == zero[RVHHChain(Stair,%d,[%s])];",2*n,weight))
        emit("show rvHHLength[1](hh_boundary);")
    }

    for (i = 1; i <= line_count; i++) print lines[i]
}

function cuts_between(first, last,    k, cuts) {
    for (k = first; k <= last; k++) cuts = append(cuts, k)
    return cuts
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
