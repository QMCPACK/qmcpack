
#include "include/utils.h"
#include <stdarg.h>

void adios_error (enum ADIOS_ERRCODES errcode, char *fmt, ...)
{
    va_list ap;
    adios_errno = (int)errcode;
    va_start(ap, fmt);
    (void) vsnprintf(aerr, ERRMSG_MAXLEN, fmt, ap);
    va_end(ap);
    log_error("%s", aerr);
}


