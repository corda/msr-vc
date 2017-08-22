//#include "Types.h"
#ifndef NDEBUG
  #define UNDEFINENDEBUG
  #define NDEBUG
#endif
#include <msbignum.h>
#ifdef UNDEFINENDEBUG
  #undef NDEBUG
  #undef UNDEFINENDEBUG
#endif


/* From ext4bignum.c in the QTD: Enigma/samples/bignum/expon/ */
#ifndef NT_SUCCESS
#define NT_SUCCESS(Status) ((NTSTATUS)(Status) >= 0)
#endif

#ifdef __cplusplus
extern "C" {
#endif

BOOL 
WINAPI 
random_bytes(
    BYTE* buffer, 
    const size_t count, 
    PVOID pContext)
{
//    bigctx_t                *pBignumCtxt = (bigctx_t*)pContext;
    BOOL                    fRet = FALSE;

    unsigned int            i = 0;

    NTSTATUS    status;

    //if (NULL == pBignumCtxt)
    //{
    //    fRet = FALSE;
    //    goto cleanup;
    //}

    status = BCryptGenRandom(   NULL,
                                buffer,
                                (ULONG)count,
                                BCRYPT_USE_SYSTEM_PREFERRED_RNG);
    if(!NT_SUCCESS(status))
    {
        fRet = FALSE;
        goto cleanup;
    }

    fRet = TRUE;

cleanup:
    return fRet;
}

// Allocation functions for bignum
void* WINAPI mp_alloc_temp (
  __in      DWORDREGC cb,
  __in_opt  LPCSTR    pszSource_info,
  __in      bigctx_t  *pCtx)
{
    void    *pVoid;

    UNREFERENCED_PARAMETER(pszSource_info);

    pVoid = malloc((DWORD)cb);

    if (NULL == pVoid)
    {
        SetMpErrno_clue(MP_ERRNO_NO_MEMORY, "mp_alloc_temp", pCtx);
    }

    memset(pVoid, 0, cb);

    return pVoid;
}


void WINAPI mp_free_temp(
  __in                  void     *pVoid,
  __in_opt              LPCSTR    pszSource_info,
  __in                  bigctx_t  *pCtx)
{
    UNREFERENCED_PARAMETER(pszSource_info);
    UNREFERENCED_PARAMETER(pCtx);

    free(pVoid);
}


void WINAPI SetMpErrno(mp_errno_tc code, bigctx_t *pBignumCtxt)
{
    if (NULL != pBignumCtxt)
    {
        pBignumCtxt->latest_errno = code;
    }
}

void WINAPI SetMpErrno_clue1(__in mp_errno_tc code, __in_opt const char *hint, PBIGCTX_ARG)
{
    UNREFERENCED_PARAMETER(hint);

    SetMpErrno(code, PBIGCTX_PASS);
}

#ifdef __cplusplus
}
#endif