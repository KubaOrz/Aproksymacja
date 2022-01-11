#include "piv_ge_solver.h"
#include <stdlib.h>

int
piv_ge_solver (matrix_t * eqs)
{
  if (eqs != NULL) {
    gradient(eqs);
    //pivot_ge_in_situ_matrix (eqs); // zamiana na nasza 
    if (bs_matrix (eqs) == 0) {
      return 0;
    }
    else {
      return 1;
    }
  }
  else
    return 1;
}
