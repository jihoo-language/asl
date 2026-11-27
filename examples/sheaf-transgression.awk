# Generate a length-r complex on the 2r+2-point model of the r-sphere.
# Use -v r=4 -v prime=3 and optionally -v untwisted=1 or -v cellular=1.
function dim(q, point,    level, edge) {
    level = int(point / 2)
    if (q == 0)
        return 1 + (level >= r - 1)
    edge = r - 1 - q
    return (level >= edge) + (level > edge)
}
function restriction(q, a, b, row, col,    level, edge) {
    level = int(a / 2)
    if (q == 0) {
        if (!untwisted && level == r - 1 && a % 2 == 0 && b % 2 == 0 && col == 1)
            return 1
        return row == col
    }
    edge = r - 1 - q
    return row == (level == edge ? a % 2 : col)
}
function vertical(q, row, col) {
    if (q == 0 && col == 0)
        return 0
    return row == 0 ? 1 : -1
}
BEGIN {
    if (r == "")
        r = 4
    if (prime == "")
        prime = 3
    if (r != int(r) || r < 2 || r > 255 || prime != int(prime) ||
            prime < 2 || prime > 4294967295 || (untwisted != "" && untwisted != 0 && untwisted != 1) ||
            (cellular != "" && cellular != 0 && cellular != 1)) {
        print "Use an integer r in 2..255, a u32 prime, and optional 0 or 1 for untwisted and cellular." > "/dev/stderr"
        exit 1
    }

    printf "space Sphere {"
    separator = ""
    for (a = 0; a < 2 * r; a++) {
        next_level = 2 * (int(a / 2) + 1)
        for (b = next_level; b < next_level + 2; b++) {
            printf "%sa%d<a%d", separator, a, b
            separator = ","
        }
    }
    print "};"

    for (q = 0; q < r; q++) {
        printf "sheaf F%d[Sphere,%.0f] {", q, prime
        for (a = 0; a < 2 * r + 2; a++)
            printf "%sa%d:%d", a ? "," : "", a, dim(q, a)
        printf "} {"
        separator = ""
        for (a = 0; a < 2 * r; a++) {
            next_level = 2 * (int(a / 2) + 1)
            for (b = next_level; b < next_level + 2; b++) {
                printf "%sa%d<a%d:[", separator, a, b
                separator = ","
                for (row = 0; row < dim(q, b); row++)
                    for (col = 0; col < dim(q, a); col++)
                        printf "%s%d", row || col ? "," : "", restriction(q, a, b, row, col)
                printf "]"
            }
        }
        print "};"
    }

    for (q = 0; q < r - 1; q++) {
        printf "shmap d%d[F%d,F%d] {", q, q, q + 1
        for (a = 0; a < 2 * r + 2; a++) {
            printf "%sa%d:[", a ? "," : "", a
            for (row = 0; row < dim(q + 1, a); row++)
                for (col = 0; col < dim(q, a); col++)
                    printf "%s%d", row || col ? "," : "", vertical(q, row, col)
            printf "]"
        }
        print "};"
    }
    printf "shcomplex %s {", cellular ? "Ordinary" : "K"
    for (q = 0; q < r - 1; q++)
        printf "F%d,d%d,", q, q
    printf "F%d};\n", r - 1
    if (cellular) {
        printf "shresolution R[Sphere,%.0f] {", prime
        for (a = 0; a < 2 * r + 2; a++)
            printf "%se%d:(%d,a%d)", a ? "," : "", a, int(a / 2), a
        printf "} {"
        for (a = 2; a < 2 * r + 2; a++) {
            previous = 2 * (int(a / 2) - 1)
            printf "%se%d:[e%d:1,e%d:-1]", (a > 2 ? "," : ""), a, previous, previous + 1
        }
        print "} {e0:1,e1:1};"
        print "shcellular K[R,Ordinary];"
    }

    printf "shedim[K,2,0,%d];\n", r - 1
    printf "let v2 = shegen[K,2,0,%d,0];\n", r - 1
    for (j = 2; j < r; j++) {
        printf "assert shED(v%d) == zero[ShE(K,%d,%d,%d)];\n", j, j, j, r - j
        printf "let v%d = shENext(v%d);\n", j + 1, j
    }
    printf "let hit = shED(v%d);\n", r
    print "show hit;"
    printf "assert shED(hit) == zero[ShE(K,%d,%d,%d)];\n", r, 2 * r, 1 - r
    printf "shedim[K,%d,0,%d];\n", r + 1, r - 1
    printf "shedim[K,%d,%d,0];\n", r + 1, r
    if (!untwisted)
        printf "assert shENext(hit) == zero[ShE(K,%d,%d,0)];\n", r + 1, r
}
