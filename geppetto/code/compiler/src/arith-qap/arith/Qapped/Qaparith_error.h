#ifndef __ARITH_ERROR_H
#define __ARITH_ERROR_H

#define ERR_SUCCESS                 1  /* No error */

#define ERR_OUT_OF_MEMORY          -1 /* Out of memory, *alloc failed */
#define ERR_INVALID_LIMB_SIZE      -2 /* Limb size for the modulus is <= 0 */
#define ERR_PARAMETERS_NOT_SET     -3 /* Initial parameters for the modulus have not been set. */
#define ERR_EVEN_MODULUS           -4 /* The modulus is even. */

#endif /* __ARITH_ERROR_H */
