# Emit a parallel-arm restriction deformation with its first obstruction in order r.
BEGIN {
    if (r == "") r = 5
    if (p == "") p = 5
    if (r !~ /^[0-9]+$/ || r < 2 || p !~ /^[0-9]+$/ || p < 2 || p > 4294967295) {
        print "r must be an integer at least 2 and p must be a prime fitting u32." > "/dev/stderr"
        exit 2
    }
    printf "space P {x0"
    for (i=1; i<=r; i++) printf "<x%d", i
    printf ",x0<z<x%d};\n", r
    printf "sheaf F[P,%.0f] {", p
    for (i=0; i<=r; i++) printf "x%d:1,", i
    printf "z:0} {"
    for (i=0; i<r; i++) printf "x%d<x%d:[0],", i, i+1
    printf "x0<z:[],z<x%d:[]};\n", r
    print "let a0 = zero[ShHomC(F,F,1)];"
    for (i=0; i<r; i++) {
        printf "let a%d = add(a%d,shHomTerm[F,F,[x%d,x%d],0,0]);\n", i+1, i, i, i+1
    }
    print "let j0 = shDefStart[F];"
    printf "let j1 = shDefAppend(j0,a%d);\n", r
    for (i=2; i<r; i++) {
        printf "let j%d = shDefLift(j%d);\n", i, i-1
        printf "assert shD(shDefTerm[%d](j%d)) == shScale[-1](shDefCurvature(j%d));\n", i, i, i-1
    }
    printf "let c = shDefCurvature(j%d);\n", r-1
    printf "let other = shHomTerm[F,F,[x0,z,x%d],0,0];\n", r
    printf "let obstruction = shDefObstruction(j%d);\n", r-1
    print "assert shD(c) == zero[ShHomC(F,F,3)];"
    print "assert obstruction == shClass(c);"
    print "assert obstruction == shScale[-1](shClass(other));"
    print "show obstruction;"
}
