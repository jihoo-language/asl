# Emit a declared module with a nonzero d_m in the localized Bockstein sequence.
BEGIN {
    if (m == "") m = 8
    if (m !~ /^[0-9]+$/ || m + 0 < 1) {
        print "m must be a positive integer" > "/dev/stderr"
        exit 2
    }

    emit("qmodule Flash[0, 1] {")
    for (j = 1; j <= m; j++) {
        emit(sprintf("  a%d:%d, b%d:%d%s", j, 2*j, j, 2*j-1, j == m ? "" : ","))
    }
    emit("} {")
    for (j = 1; j <= m; j++) {
        rule = sprintf("  0(a%d)=>[b%d]", j, j)
        if (j > 1) rule = rule sprintf(", 1(a%d)=>[b%d]", j, j-1)
        emit(rule (j == m ? "" : ","))
    }
    emit("};")

    emit("let lift1 = vLocal(qterm[One,0,0,Flash,2,0,[0,1]]{0,0});")
    for (j = 2; j <= m; j++) {
        emit(sprintf("let a%d = vLocal(qterm[One,0,0,Flash,%d,0,[0,1]]{%d,0});", j, 2*j, j-1))
        emit(sprintf("let lift%d = add(lift%d, vShift[%d](a%d));", j, j-1, 1-j, j))
    }
    emit(sprintf("let corrected = lift%d;", m))
    emit(sprintf("let target = vShift[%d](vLocal(qterm[One,0,0,Flash,%d,0,[0,1]]{%d,0}));", 1-m, 2*m-1, m))
    emit("assert qD(corrected) == target;")
    emit(sprintf("let source = vbClass[%d,0](corrected);", m))
    emit(sprintf("assert source == vbockgen[One,Flash,[0,1],%d,0,2,0,0];", m))
    emit("let previous = vbPrev(source);")
    emit(sprintf("assert vbD(previous) == zero[VBock(One,Flash,[0,1],%d,1,2,%d)];", m-1, m-1))
    emit("assert vbNext(previous) == source;")
    emit("let result = vbD(source);")
    emit(sprintf("assert result == vbClass[%d,%d](target);", m, m))
    emit(sprintf("assert result == vbockgen[One,Flash,[0,1],%d,1,2,%d,0];", m, m))
    emit(sprintf("assert vbNext(result) == zero[VBock(One,Flash,[0,1],%d,1,2,%d)];", m+1, m))
    emit("show result;")

    for (j = 1; j <= line_count; j++) print lines[j]
}

function emit(line) {
    # ASL accepts at most 512 KiB of source. Check before writing any output.
    source_bytes += length(line) + 1
    if (source_bytes > 512 * 1024) {
        print "generated program exceeds the ASL source limit of 512 KiB" > "/dev/stderr"
        exit 2
    }
    lines[++line_count] = line
}
