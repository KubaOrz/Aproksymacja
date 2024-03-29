#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdio.h>

typedef struct {
	int rn;
    int cn;
    double *e;
} matrix_t;

//Nasza nowa struktura
typedef struct {
    int rn;
    int cn;
    double **e;
}   macierz;

matrix_t * make_matrix( int rn, int cn );

matrix_t * read_matrix( FILE *in );

void write_matrix( matrix_t *, FILE *out );

void put_entry_matrix( matrix_t *, int i, int j, double val );

void add_to_entry_matrix( matrix_t *, int i, int j, double val );

double get_entry_matrix( matrix_t *, int i, int j );

matrix_t * copy_matrix( matrix_t *s );

matrix_t * transpose_matrix( matrix_t * s );

void xchg_rows( matrix_t *m, int i, int j );

void xchg_cols( matrix_t *m, int i, int j );

matrix_t * mull_matrix( matrix_t *, matrix_t * );

matrix_t * ge_matrix( matrix_t * );

int bs_matrix( matrix_t * );

matrix_t * pivot_ge_matrix( matrix_t *, int *row_per );

void pivot_ge_in_situ_matrix( matrix_t * );

void gradient(matrix_t *); /* metoda gradientow sprzezonych */

matrix_t * symm_pivot_ge_matrix( matrix_t *, int *per );

int *pivot_get_inv_per( matrix_t *, int *row_per );

macierz *utworz(int wiersze, int kolumny);

macierz *pomnoz(macierz A, macierz x);

macierz *transponuj(macierz A);

macierz *odejmij(macierz A,macierz x);

macierz *dodaj(macierz A,macierz x);

macierz *pomnoz_przez_liczbe(double n, macierz A);

double wartosc(macierz A);

#endif
