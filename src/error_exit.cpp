/* Error handling routine. Used exactly like printf. */

#include <iostream>
#include <stddef.h>
#include <stdlib.h>
#include <stdarg.h>
#include <signal.h>

void error_exit(const char *error_msg, ...) {
  va_list args;

  va_start (args, error_msg);
  /* Print error message to standard output and exit. */
  vfprintf (stderr, error_msg, args);
  va_end (args);

  exit(1);
}

void warning(const char *warning_msg) {
  std::cout << "  WARNING: " << warning_msg << std::endl;
}


