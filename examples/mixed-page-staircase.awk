# Emit a finite mixed complex and check the high differential it produces.
BEGIN {
    if (r == "") r = 20
    if (p == "") p = 3
    if (control == "") control = 0
    if (bars == "") bars = 0
    if (direct == "") direct = 0
    if (hom == "") hom = 0
    if (r !~ /^[0-9]+$/ || r < 1 ||
        p !~ /^[0-9]+$/ || p < 2 || p > 4294967295 ||
        control !~ /^[01]$/ || bars !~ /^[01]$/ || direct !~ /^[01]$/ || hom !~ /^[01]$/) {
        print "r must be positive and p must be a prime fitting u32. Set control, bars, hom, and direct to 0 or 1." > "/dev/stderr"
        exit 2
    }
    emit("mixedcomplex Stair[" p "] {")
    for (i=0; i<r; i++) emit(sprintf("    x%d:%d, y%d:%d%s", i,2*i,i+1,2*i+1,(i+1==r?"":",")))
    emit("} {")
    for (i=1; i<r; i++) emit(sprintf("    x%d => [(y%d,-1)]%s",i,i,(i+1==r?"":",")))
    emit("} {")
    for (i=0; i<r-control; i++) emit(sprintf("    x%d => [(y%d,1)]%s",i,i+1,(i+1==r-control?"":",")))
    emit("};")
    C = hom ? "Maps" : "Stair"
    y = sprintf("mixU[%d](mixchain[Stair,%d][(y%d,1)])",r,2*r-1,r)
    if (hom) {
        emit("mixedcomplex Unit[" p "] { x:0 } {} {};")
        emit("mixhom Maps[Unit,Stair];")
        emit("mixhom Ids[Unit,Unit];")
        emit("let actual = " y ";")
        emit("let input = mixU[0](mixchain[Unit,0][(x,1)]);")
        y = sprintf("mixU[%d](mixunit[Maps,x,y%d])",r,r)
    }

    emit(sprintf("let x = mixgen[%s,%d,0,0,0];", C,r))
    emit("let y = " y ";")
    for (i=(direct ? r-1 : 0); i<r; i++) {
        emit(sprintf("assert mixD(mixgen[%s,%d,0,0,0]) == zero[MixedE(%s,%d,-1,%d)];", C,i,C,i,i))
    }
    emit(sprintf("assert mixPrev(x) == mixgen[%s,%d,0,0,0];", C,r-1))
    if (control) {
        emit(sprintf("assert mixD(mixPageRep(x)) == zero[MixedU(%s,-1)];",C))
        emit(sprintf("assert mixNext(x) == mixgen[%s,%d,0,0,0];", C,r+1))
        if (hom) {
            emit("let candidate = mixPageRep(x);")
            emit("let evaluated = mixEval(candidate,input);")
            emit("assert mixD(evaluated) == zero[MixedU(Stair,-1)];")
            emit("assert mixEval(mixComp[Maps](candidate,mixU[0](mixId[Ids])),input) == evaluated;")
            emit("assert mixPageClass[0,0](evaluated) == mixgen[Stair,0,0,0,0];")
        }
        emit("show mixNext(x);")
    } else {
        emit("assert mixD(mixPageRep(x)) == y;")
        emit(sprintf("assert mixD(x) == mixgen[%s,%d,-1,%d,0];",C,r,r))
        emit(sprintf("assert mixNext(mixD(x)) == zero[MixedE(%s,%d,-1,%d)];",C,r+1,r))
        if (hom) {
            emit("let candidate = mixPageRep(x);")
            emit("assert mixEval(mixD(candidate),input) == actual;")
            emit("assert mixD(mixEval(candidate,input)) == actual;")
            emit("assert mixEval(mixComp[Maps](candidate,mixU[0](mixId[Ids])),input) == mixEval(candidate,input);")
        }
        emit(sprintf("mixdim[%s,%d,0,0];",C,r+1))
        emit("show mixD(x);")
    }

    if (bars) {
        if (control) {
            emit(sprintf("assert mixBarRep[%s,0] == mixPageRep(x);",C))
            emit(sprintf("assert mixD(mixBarRep[%s,1]) == zero[MixedU(%s,%d)];",C,C,2*r-2))
        } else {
            emit(sprintf("assert mixBarNull[%s,0] == mixPageRep(x);",C))
            emit(sprintf("assert mixD(mixBarNull[%s,0]) == mixU[%d](mixBarRep[%s,0]);",C,r,C))
        }
        emit(sprintf("mixbars[%s];",C))
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
