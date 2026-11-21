# Emit two 127-fold defining systems over C64. The second graph omits
# one edge. Both graphs induce exactly the same 127 adjacent inputs.
function edge(i, j, omit) {
    if (omit && i == 61 && j == 64) return 0
    # Between two internal layers the adjacency matrix is [1 1; 1 0].
    return i == 0 || j == n || i % 2 || j % 2
}

function cochain(i, j) {
    return paths[i,j] ? "b" (layer[j] - layer[i]) : "z"
}

BEGIN {
    m = 64
    n = 2 * m - 1
    for (j = 0; j <= n; j++) layer[j] = int((j + 1) / 2)
    printf "let z = bar[%d]{};\n", m
    for (k = 1; k < m; k++) printf "let b%d = binom[%d, %d];\n", k, m, k

    for (omit = 0; omit <= 1; omit++) {
        for (i = 0; i < n; i++) {
            paths[i,i] = 1
            for (j = i + 1; j <= n; j++) {
                # The only possible predecessors are in the previous layer.
                first = layer[j] == 1 ? 0 : 2 * layer[j] - 3
                last = layer[j] == 1 ? 0 : first + 1
                count = 0
                for (p = first; p <= last; p++) {
                    if (p >= i && edge(p, j, omit)) count += paths[i,p]
                }
                paths[i,j] = count % 2
            }
        }

        printf "show massey("
        for (i = 0; i < n; i++) printf "%s%s", (i ? ", " : ""), cochain(i, i + 1)
        print "){"
        for (gap = 2; gap < n; gap++) {
            for (i = 0; i + gap <= n; i++) {
                printf "%d < %d => %s,\n", i, i + gap, cochain(i, i + gap)
            }
        }
        print "};"
    }
}
