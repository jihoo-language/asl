# Emit a finite mixed complex and check the high differential it produces.
BEGIN {
    if (r == "") r = 20
    if (p == "") p = 3
    if (control == "") control = 0
    if (r !~ /^[0-9]+$/ || r < 1 ||
        p !~ /^[0-9]+$/ || p < 2 || p > 4294967295 ||
        control !~ /^[01]$/) {
        print "r must be positive, p a prime fitting u32, and control 0 or 1" > "/dev/stderr"
        exit 2
    }
    emit("mixedcomplex Stair[" p "] {")
    for (i=0; i<r; i++) emit(sprintf("    x%d:%d, y%d:%d%s", i,2*i,i+1,2*i+1,(i+1==r?"":",")))
    emit("} {")
    for (i=1; i<r; i++) emit(sprintf("    x%d => [(y%d,-1)]%s",i,i,(i+1==r?"":",")))
    emit("} {")
    for (i=0; i<r-control; i++) emit(sprintf("    x%d => [(y%d,1)]%s",i,i+1,(i+1==r-control?"":",")))
    emit("};")
    emit(sprintf("let x = mixgen[Stair,%d,0,0,0];", r))
    emit(sprintf("let y = mixU[%d](mixchain[Stair,%d][(y%d,1)]);", r,2*r-1,r))
    for (i=0; i<r; i++) {
        emit(sprintf("assert mixD(mixgen[Stair,%d,0,0,0]) == zero[MixedE(Stair,%d,-1,%d)];", i,i,i))
    }
    emit(sprintf("assert mixPrev(x) == mixgen[Stair,%d,0,0,0];", r-1))
    if (control) {
        emit("assert mixD(mixPageRep(x)) == zero[MixedU(Stair,-1)];")
        emit(sprintf("assert mixNext(x) == mixgen[Stair,%d,0,0,0];", r+1))
        emit("show mixNext(x);")
    } else {
        emit("assert mixD(mixPageRep(x)) == y;")
        emit(sprintf("assert mixD(x) == mixgen[Stair,%d,-1,%d,0];",r,r))
        emit(sprintf("assert mixNext(mixD(x)) == zero[MixedE(Stair,%d,-1,%d)];",r+1,r))
        emit(sprintf("mixdim[Stair,%d,0,0];",r+1))
        emit("show mixD(x);")
    }

    for (i=1; i<=line_count; i++) print lines[i]
}

function emit(line) {
    source_bytes += length(line) + 1
    if (source_bytes > 512 * 1024) {
        print "generated program exceeds the ASL source limit of 512 KiB" > "/dev/stderr"
        exit 2
    }
    lines[++line_count] = line
}
