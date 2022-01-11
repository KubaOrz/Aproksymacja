#include "matrix.h"
#include <stdlib.h>
#include <math.h>

matrix_t *
pivot_ge_matrix(matrix_t *a, int *row_per)
{
  matrix_t *c = copy_matrix(a);
  if (c != NULL)
  {
    int i, j, k;
    int cn = c->cn;
    int rn = c->rn;
    double *e = c->e;
    for (i = 0; i < rn; i++)
      row_per[i] = i;
    for (k = 0; k < rn - 1; k++)
    {              /* eliminujemy (zerujemy) kolumnę nr k */
      int piv = k; /* wybór eleemntu dominującego - maks. z k-tej kol., poniżej diag */
      for (i = k + 1; i < rn; i++)
        if (fabs(*(e + i * cn + k)) > fabs(*(e + piv * cn + k)))
          piv = i;
      if (piv != k)
      { /* jeśli diag. nie jest pivtem - wymień wiersze */
        int tmp;
        xchg_rows(c, piv, k);
        tmp = row_per[k];
        row_per[k] = row_per[piv];
        row_per[piv] = tmp;
      }
      for (i = k + 1; i < rn; i++)
      { /* pętla po kolejnych
                                           wierszach poniżej diagonalii k,k */
        double d = *(e + i * cn + k) / *(e + k * cn + k);
        for (j = k; j < cn; j++)
          *(e + i * cn + j) -= d * *(e + k * cn + j);
      }
    }
  }
  return c;
}

void pivot_ge_in_situ_matrix(matrix_t *c)
{
  int i, j, k;
  int cn = c->cn;
  int rn = c->rn;
  double *e = c->e;
  for (k = 0; k < rn - 1; k++)
  {              /* eliminujemy (zerujemy) kolumnę nr k */
    int piv = k; /* wybór elementu dominującego - maks. z k-tej kol., poniżej diag */
    for (i = k + 1; i < rn; i++)
      if (fabs(*(e + i * cn + k)) > fabs(*(e + piv * cn + k)))
        piv = i;
    if (piv != k)
    { /* jeśli diag. nie jest pivtem - wymień wiersze */
      xchg_rows(c, piv, k);
    }
    for (i = k + 1; i < rn; i++)
    { /* pętla po kolejnych
                                           wierszach poniżej diagonalii k,k */
      double d = *(e + i * cn + k) / *(e + k * cn + k);
      for (j = k; j < cn; j++)
        *(e + i * cn + j) -= d * *(e + k * cn + j);
    }
  }
}
void gradient(matrix_t *eqs){
    fprintf(stderr, "TU");
    double **A = malloc( eqs->rn *sizeof *A );
    double *b = malloc (eqs->rn *sizeof *b);

    if(A==NULL || b==NULL)
      printf("PAMIEC\n");
    for(int i=0; i<eqs->cn; i++){
      if((A[i]= malloc(eqs->cn *sizeof *A[i]))==NULL)
        printf("PAMIEC\n");
    }

    int row, col;

    for(row = 0; row < eqs->rn; row++){
      for(col = 0; col < eqs->cn;col++ ){
          A[row][col] = eqs->e[row * eqs->cn + col];
          printf("%lf ", A[row][col]);
      }
      printf("\n");
    }
}

matrix_t *
symm_pivot_ge_matrix(matrix_t *a, int *row_per)
{
  matrix_t *c = copy_matrix(a);
  if (c != NULL)
  {
    int i, j, k;
    int cn = c->cn;
    int rn = c->rn;
    double *e = c->e;
    for (i = 0; i < rn; i++)
      row_per[i] = i;
    for (k = 0; k < rn - 1; k++)
    {              /* eliminujemy (zerujemy) kolumnę nr k */
      int piv = k; /* wybór eleemntu dominującego - maks. z k-tej kol., poniżej diag */
      for (i = k + 1; i < rn; i++)
        if (fabs(*(e + i * cn + k)) > fabs(*(e + piv * cn + k)))
          piv = i;
      if (piv != k)
      { /* jeśli diag. nie jest pivtem - wymień wiersze */
        int tmp;
        xchg_rows(c, piv, k);
        xchg_cols(c, piv, k);
        tmp = row_per[k];
        row_per[k] = row_per[piv];
        row_per[piv] = tmp;
      }
      for (i = k + 1; i < rn; i++)
      { /* pętla po kolejnych
                                           wierszach poniżej diagonalii k,k */
        double d = *(e + i * cn + k) / *(e + k * cn + k);
        for (j = k; j < cn; j++)
          *(e + i * cn + j) -= d * *(e + k * cn + j);
      }
    }
  }
  return c;
}

int *pivot_get_inv_per(matrix_t *m, int *row_per)
{
  /* odtwarzamy oryginalną kolejność niewiadowmych */
  int *iper = malloc(m->rn * sizeof *iper);
  int i;

  for (i = 0; i < m->rn; i++)
    iper[row_per[i]] = i;

  return iper;
}
