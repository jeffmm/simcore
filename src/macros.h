#ifndef _MACROS_H
#define _MACROS_H

#define SQR(x)             ((x) * (x))
#define CUBE(x)            ((x) * (x) * (x))
#define NINT(x)            ((x) < 0.0 ? (int) ((x) - 0.5) : (int) ((x) + 0.5))
#define ABS(x)             ((x) < 0 ? -(x) : (x))
#define MAX(x,y)           ((x) > (y) ? (x) : (y))
#define MIN(x,y)           ((x) < (y) ? (x) : (y))
#define SIGN(a,b)          ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SIGNOF(x)          ((x) >= 0.0 ? 1 : -1)

// Debugging print macros!
#ifdef NDEBUG
#define DPRINTF(M, ...)
#else
#define DPRINTF(M, ...) do { fprintf(stdout, "" M "", ##__VA_ARGS__); } while(0)
#endif

#endif //_MACROS_H
