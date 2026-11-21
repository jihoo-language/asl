# Attach a disk to a bouquet of circles along a word, then emit the requested
# Massey products. Uppercase letters traverse their generator backwards.
# awk -v word=aaaa -v queries='aa aaa aaaa aaaaa' -f examples/bracket-relator.awk
function relation(x, y) { print "  " x " < " y "," }
function fail(message) { print message > "/dev/stderr"; exit 2 }

# The character marks o < f_x in the bouquet. Pull it back to the mapping
# cylinder, then integrate around the attaching circle to extend over the cone.
function character(x,    i, j, parity, mark) {
    mark = "f_" x
    printf "let %s = cocycle[K, 1]{\n", x
    relation("o", mark)
    parity = 0
    for (i = 0; i < n; i++) {
        j = (i + 1) % n
        if (vertex[i] == "o") relation("v" i, mark)
        if (parity) relation("v" i, "apex")
        if (vertex[i] == "o" && edge[i] == mark) {
            relation("v" i, "e" i)
            parity = 1 - parity
        }
        if (parity) relation("e" i, "apex")
        if (vertex[j] == "o" && edge[i] == mark) {
            relation("v" j, "e" i)
            parity = 1 - parity
        }
    }
    print "};"
}

BEGIN {
    alphabet = "abcdefghijklmnopqrstuvwxyz"
    if (word == "") word = "asASabABcbaBAC"
    if (queries == "" && word == "asASabABcbaBAC") queries = "abc bca abca"
    if (queries == "" && word == "abAbaBAB") queries = "abb bab"
    if (word !~ /^[a-zA-Z]+$/ || queries !~ /^[a-z]+( [a-z]+)*$/)
        fail("supply an ASCII relator word and space-separated bracket words")
    n = 2 * length(word)
    for (i = 0; i < n; i += 2) {
        letter = substr(word, i / 2 + 1, 1)
        x = tolower(letter)
        occurrences[x]++
        vertex[i] = "o"
        vertex[i + 1] = "p_" x
        edge[i] = (x == letter ? "f_" : "r_") x
        edge[i + 1] = (x == letter ? "r_" : "f_") x
    }
    count = split(queries, query, " ")
    for (q = 1; q <= count; q++) {
        if (length(query[q]) < 2) fail("each bracket needs at least two inputs")
        for (i = 1; i <= length(query[q]); i++) {
            x = substr(query[q], i, 1)
            if (!occurrences[x] || occurrences[x] % 2)
                fail("the single-letter character " x " does not exist for this presentation")
            needed[x] = 1
        }
    }
    print "space K {"
    for (i = 1; i <= length(alphabet); i++) {
        x = substr(alphabet, i, 1)
        if (occurrences[x]) {
            relation("o", "f_" x)
            relation("o", "r_" x)
            relation("p_" x, "f_" x)
            relation("p_" x, "r_" x)
        }
    }
    for (i = 0; i < n; i++) {
        relation("v" i, "e" i)
        relation("v" ((i + 1) % n), "e" i)
        relation("v" i, vertex[i])
        relation("e" i, edge[i])
        relation("e" i, "apex")
    }
    print "};\nbetti K;"
    for (i = 1; i <= length(alphabet); i++) {
        x = substr(alphabet, i, 1)
        if (needed[x]) character(x)
    }
    for (q = 1; q <= count; q++) {
        printf "show bracket("
        for (i = 1; i <= length(query[q]); i++)
            printf "%s%s", (i == 1 ? "" : ", "), substr(query[q], i, 1)
        print ");"
    }
}
