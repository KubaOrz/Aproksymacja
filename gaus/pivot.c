#include "matrix.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

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
    macierz A, x ,r, p ,b;
    A.rn = eqs -> rn;
    A.cn = eqs -> cn - 1;
    A.e = malloc( eqs->rn *sizeof *A.e );
    b.rn = eqs -> rn;
    b.cn = 1;
    b.e = malloc(eqs->rn * sizeof *b.e);
    //double *b = malloc (eqs->rn *sizeof *b);
    x.rn = eqs->rn;
    x.cn = 1;
    x.e = malloc(x.rn * sizeof *x.e);
    //double **x = malloc(eqs -> rn * sizeof *x);
    r.rn = eqs -> rn;
    r.cn = 1;
    r.e = malloc(eqs->rn * sizeof *r.e);
    p.rn = eqs -> rn;
    p.cn = 1;
    p.e = malloc(eqs -> rn * sizeof *p.e);
    //double *wiersz = malloc((eqs->cn-1) *sizeof *wiersz);
    int k =0;
    double alfa;

    if(A.e==NULL || x.e==NULL || r.e==NULL || p.e ==NULL)
      printf("PAMIEC\n");

    //memset(x.e, 0, eqs -> rn);

    for(int i=0; i<eqs->rn; i++){
      if((A.e[i]= malloc(A.cn *sizeof *A.e[i]))==NULL)
        printf("PAMIEC\n");
      if((x.e[i]= malloc(x.cn *sizeof *x.e[i]))==NULL)
        printf("PAMIEC\n");
      if((r.e[i]= malloc(r.cn *sizeof *r.e[i]))==NULL)
        printf("PAMIEC\n");
      if((p.e[i]= malloc(p.cn *sizeof *p.e[i]))==NULL)
        printf("PAMIEC\n");
      if((b.e[i]= malloc(b.cn *sizeof *p.e[i]))==NULL)
        printf("PAMIEC\n");
    }
    

    for (int i = 0; i < x.rn; i++)
      x.e[i][0] = 0;

    int row, col;

    for(row = 0; row < eqs->rn; row++){
      for(col = 0; col < eqs->cn - 1;col++ ){
          A.e[row][col] = eqs->e[row * eqs->cn + col];
          printf("%lf ", A.e[row][col]);
      }
      b.e[row][0] = eqs -> e[row * eqs -> cn + eqs -> cn - 1];
      printf("\n");
    }

    while(1){ /* kolejne iteracje */
      //zapisanie residuum
      /*for (int i = 0; i < eqs -> rn; i++) {
          r[i] = A[i][eqs->cn] - pomnoz(A, x, i, eqs->cn);
          p[i] = r[i];
      }*/
      r = *odejmij(b, *pomnoz(A, x));
      p = r;
      for(int i =0 ; i<r.rn; i++)
        fprintf(stderr, "%lf ",r.e[i][0]);

    //alfa liczymy
    for(int i =0 ; i<r.rn; i++)
        fprintf(stderr, "%lf ",r.e[i][0]);
    alfa = (*pomnoz(*transponuj(r), r)).e[0][0] / (*pomnoz(*transponuj(p), *pomnoz(A,p))).e[0][0];
    printf("ALFA%f\n",alfa);
    break;
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
