#include "two-matrix-ifc.h"
#include "pinocchio.h"

void outsource(struct bank_Input *input, struct bank_Output *output)
{
  int i, j, k;
  int t;
  for (i=0; i<SIZE; i+=1)
    {
      for (j=0; j<SIZE; j+=1)
        {
          t =0;
          for (k=0; k<SIZE; k+=1)
            {
              t += MAT(input->a, i, k) * MAT(input->b, k, j);
            }
		  MAT(output->r, i, j) =
#if TRUNCATE==0
			  t;
#else	
			  bound(t, 0, 1000); // safely ensuring that we do not overflow (proving fails if we do)
#endif
        }
    }
}

