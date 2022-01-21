#include <stdio.h>
#include "matrix.h"
#include <stdlib.h>

int main(){
    macierz A;
    macierz b;
    macierz *c;
    A.e = malloc(3 * sizeof *A.e);
    A.rn = 3;
    A.cn = 2;
    for( int i =0 ;i<3;i++){
        A.e[i] = malloc(2* sizeof (double));
    }
    A.e[0][0] = 1;
    A.e[0][1] = 2;
    A.e[1][0] = 3;
    A.e[1][1] = 4;
    A.e[2][0] = 0;
    A.e[2][1] = 5;
    for(int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            printf("%lf ", A.e[i][j]);
        }
        printf("\n");
    }
    c = transponuj(A);
    
    for(int i = 0; i < 2; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%lf ", c->e[i][j]);
        }
        printf("\n");
    }
}