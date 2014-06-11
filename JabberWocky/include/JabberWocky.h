
#ifndef JABBER_WOCKY_H
#define JABBER_WOCKY_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <mpi.h>
#include <assert.h>
#include <string.h>
#include <netinet/in.h>
#include <dirent.h>
#include <stdarg.h>

enum ADIOS_DATATYPES {adios_unknown = -1             /* (size) */

                     ,adios_byte = 0                 /* (1) */
                     ,adios_short = 1                /* (2) */
                     ,adios_integer = 2              /* (4) */
                     ,adios_long = 4                 /* (8) */

                     ,adios_unsigned_byte = 50       /* (1) */
                     ,adios_unsigned_short = 51      /* (2) */
                     ,adios_unsigned_integer = 52    /* (4) */
                     ,adios_unsigned_long = 54       /* (8) */

                     ,adios_real = 5                 /* (4) */
                     ,adios_double = 6               /* (8) */
                     ,adios_long_double = 7          /* (16) */

                     ,adios_string = 9               /* (?) */
                     ,adios_complex = 10             /* (8) */
                     ,adios_double_complex = 11      /* (16) */
                     };

enum ADIOS_FLAG {adios_flag_unknown = 0
                ,adios_flag_yes = 1
                ,adios_flag_no = 2
                };

enum ADIOS_BUFFER_ALLOC_WHEN {ADIOS_BUFFER_ALLOC_UNKNOWN
                             ,ADIOS_BUFFER_ALLOC_NOW
                             ,ADIOS_BUFFER_ALLOC_LATER
                             };


#endif
