# Use -v kind=sphere -v dim=6 -v prime=3 or -v kind=rp2 -v prime=2.
function addmark(name) { mark[n++] = name }
function edge(a,b,key) {
    key = a SUBSEP b
    if (!(key in seen)) {
        seen[key] = 1
        left[m] = a
        right[m++] = b
    }
}
function sphere(    i,a,b) {
    for (i = 0; i <= dim; ++i) {
        addmark("a" i)
        addmark("b" i)
    }
    for (i = 0; i < dim; ++i)
        for (a = 0; a < 2; ++a)
            for (b = 0; b < 2; ++b)
                edge(mark[2*i+a], mark[2*i+2+b])
    edge("a" dim,"top")
    edge("b" dim,"top")
}
function rp2(    i,j,t,h,a,b,e,tri,k) {
    for (i = 1; i <= 6; ++i)
        addmark("f" i)
    for (i = 1; i <= 6; ++i)
        for (j = i+1; j <= 6; ++j) {
            e = "f" i j
            addmark(e)
            edge("f" i,e)
            edge("f" j,e)
        }
    k = split("123 124 135 146 156 236 245 256 345 346",tri," ")
    for (t = 1; t <= k; ++t) {
        h = "f" tri[t]
        addmark(h)
        for (i = 1; i <= 3; ++i)
            for (j = i+1; j <= 3; ++j)
                edge("f" substr(tri[t],i,1) substr(tri[t],j,1),h)
        edge(h,"top")
    }
}
function site(name,cone,    i) {
    printf "space %s {",name
    for (i = 0; i < n; ++i)
        printf "%s%s", (i ? "," : ""), mark[i]
    if (cone)
        printf ",top"
    for (i = 0; i < m; ++i)
        if (cone || right[i] != "top")
            printf ",%s<%s",left[i],right[i]
    print "};"
}
function coefficient(as_set,    i) {
    if (as_set)
        printf "setfamily A[X] {"
    else
        printf "sheaf F[X,%.0f] {",prime
    for (i = 0; i < n; ++i)
        printf "%s%s:%d", (i ? "," : ""), mark[i], (as_set ? 2 : 1)
    printf ",top:%d} {",(as_set ? 2 : 1)
    for (i = 0; i < m; ++i)
        printf "%s%s<%s:[%s]",(i ? "," : ""),left[i],right[i],(as_set ? "0,1" : "1")
    print "};"
}
BEGIN {
    n = 0
    m = 0
    if (kind == "")
        kind = "sphere"
    if (dim == "")
        dim = 6
    if (prime == "")
        prime = 3
    if ((kind != "rp2" && kind != "sphere") || dim != int(dim) || dim < 0 || dim > 8 ||
            prime != int(prime) || prime < 2 || prime > 4294967295 ||
            (witness != "" && witness != 0 && witness != 1)) {
        print "Use kind=sphere or kind=rp2, integer dim in 0..8, and a u32 prime." > "/dev/stderr"
        exit 1
    }
    if (kind == "rp2")
        rp2()
    else
        sphere()
    site("X",1)
    site("Link",0)
    printf "let f=tab X -> S {"
    for (i = 0; i < n; ++i)
        printf "%s%s=>S.lo",(i ? "," : ""),mark[i]
    print ",top=>S.hi};"
    printf "let u=tab Link -> X {"
    for (i = 0; i < n; ++i)
        printf "%s%s=>X.%s",(i ? "," : ""),mark[i],mark[i]
    print "};"
    print "let v=tab One -> S {pt=>S.lo};"
    print "let fp=\\x : Link => One.pt;"
    print "shbcsquare Q[f,u,v,fp];"
    coefficient(0)
    print "shbc B[Q,F];"
    print "shbcrank[B,pt,0];"
    if (kind == "rp2") {
        print "shbcrank[B,pt,1];"
        print "shbcrank[B,pt,2];"
    } else if (dim != 0)
        printf "shbcrank[B,pt,%d];\n",dim
    printf "shbccomma[Q,%.0f,pt,top];\n",prime
    printf "shbcuniversal[Q,%.0f];\n",prime
    if (witness && (kind == "sphere" || prime == 2))
        printf "shbcwitness[Q,%.0f,pt,top,%d,0];\n",prime,(kind=="rp2" ? 1 : dim)
    if (kind == "rp2" || dim > 0) {
        coefficient(1)
        print "setpullback APrime[A](u);"
        print "setran R[A](f);"
        print "setran RPrime[APrime](fp);"
        print "setpullback RPull[R](v);"
        print "setnat inverse[RPrime,RPull] {pt:[0,1]};"
        print "assert setComp(setRanBC[RPrime,RPull],inverse) == setId[RPrime];"
        print "assert setComp(inverse,setRanBC[RPrime,RPull]) == setId[RPull];"
    }
}
