/* This routine returns the elapsed cpu time in seconds.  The standard unix command "clock" wraps 
   after 36 minutes, but this routine doesn't (stolen from Keith Refson).

   input: none

   output: elapsed cpu time in seconds (return value) */

#include "auxiliary.h"

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>

double cpu_time(void) {
    struct rusage ru;
    //int getrusage();

    (void) getrusage(RUSAGE_SELF, &ru);

    return (ru.ru_utime.tv_sec + ru.ru_stime.tv_sec
            + 1.0e-6 * (ru.ru_utime.tv_usec + ru.ru_stime.tv_usec));
}
