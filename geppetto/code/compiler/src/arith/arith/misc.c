/* misc.c
 * by Joppe W. Bos and Michael Naehrig (mnaehrig@microsoft.com), 
 * Cryptography Research Group, MSR Redmond, 2014
 * 
 * This file is part of the ARITH library version 0.01
 *
 * DISCLAIMER:
 * This is RESEARCH code.
 * Parts of the code might be incomplete.
 * Please use at your own risk.
 * DO NOT USE IN PRODUCTION CODE!!!
 * This code is still under active development.
 */

#include "misc.h"

void show_cpu_info () {
  int CPUInfo[4] = {-1};
  char CPUBrandString[0x40];
  unsigned int nExIds, i;

  /* Calling __cpuid with 0x80000000 as the InfoType argument
   * gets the number of valid extended IDs. */
  __cpuid(CPUInfo, 0x80000000);
  nExIds = CPUInfo[0];
  memset(CPUBrandString, 0, sizeof(CPUBrandString));

  /* Get the information associated with each extended ID. */
  for (i=0x80000000; i<=nExIds; ++i) {
    __cpuid(CPUInfo, i);
        
    /* Interpret CPU brand string and cache information. */
    if  (i == 0x80000002)
      memcpy(CPUBrandString, CPUInfo, sizeof(CPUInfo));
    else if  (i == 0x80000003)
      memcpy(CPUBrandString + 16, CPUInfo, sizeof(CPUInfo));
    else if  (i == 0x80000004)
      memcpy(CPUBrandString + 32, CPUInfo, sizeof(CPUInfo));
  }

  if  (nExIds >= 0x80000004) 
    printf ("CPU information:\n\t%s\n", CPUBrandString);
}

void check_heap(char *file, int line) {
    static char *lastOkFile = "here";
    static int lastOkLine = 0;
    static int heapOK = 1;

    if (!heapOK) return;

    if (_heapchk() == _HEAPOK)
    {
    printf ("All ok\n");
        lastOkFile = file;
        lastOkLine = line;
       return;
    }

    heapOK = 0;
    printf("Heap corruption detected at %s (%d)\n", file, line);
    printf("Last OK at %s (%d)\n", lastOkFile, lastOkLine);
}
