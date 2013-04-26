// -*- c++ -*-
//
// Copyright (c) 2002-2003 Indiana University.  All rights reserved.
// Copyright (c) 1996, 1997, 1998, 2000 University of Notre Dame.
//                         All rights reserved.
//
// This file is part of the OOMPI software package.  For license
// information, see the LICENSE file in the top level directory of the
// OOMPI source distribution.
//
// $Id$
//
// OOMPI Class library
// Utility functions
//

#ifndef _OOMPI_UTIL_H_
#define _OOMPI_UTIL_H_

#include "oompi-config.h"


typedef enum
{
  SEM_INIT,
  SEM_FINAL,
  SEM_BUFFER,
  SEM_ENVIRONMENT,
  SEM_ERR_HANDLER,
  SEM_DATATYPE,
  SEM_MAX
} OOMPI_Semaphore;


class OOMPI_Util
{
public:

  // JMS Extremely weak semaphore

  inline OOMPI_Util(void)
  {
    if (util_sem == 1)
    {
      util_sem = 0;
      for (int i = 0; i < SEM_MAX; i++)
        sem_array[i] = 1;
    }
  }

  inline void Get_sem(OOMPI_Semaphore semtype);
  inline void Release_sem(OOMPI_Semaphore semtype);

protected:
  static int util_sem;
  static int sem_array[];

private:
};


//
// Inline definitions
//

void OOMPI_Util::Get_sem(OOMPI_Semaphore semtype)
{
#if OOMPI_HAVE_PTHREADS
  // put pthread stuff here....
#else
  while (sem_array[semtype] != 1)
    continue;
  sem_array[semtype] = 0;
#endif
}


void
OOMPI_Util::Release_sem(OOMPI_Semaphore semtype)
{
  if (sem_array[semtype] == 0)
    sem_array[semtype] = 1;
}

#endif
