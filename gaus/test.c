#include <stdio.h>
#include "matrix.h"
#include <math.h>
int main(int argc, char **argv) {
    FILE *in = argc > 1 ? fopen(argv[1], "r"): stdin;
    matrix_t *eqs= make_matrix(3, 4);
    for (int i = 0; i < eqs -> rn; i++) {
        for (int j = 0; j < eqs -> cn; j++) {
            fscanf(in, "%lf", &eqs -> e[i*eqs->cn +j]);
        }
    }
    sqrt(2);
    pivot_ge_in_situ_matrix(eqs);
    int x = bs_matrix(eqs);
    //gradient(eqs);
}