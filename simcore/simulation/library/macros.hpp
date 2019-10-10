#ifndef _MACROS_H_
#define _MACROS_H_

#define SQR(x) ((x) * (x))
#define CUBE(x) ((x) * (x) * (x))
#define NINT(x) ((x) < 0.0 ? (int)((x)-0.5) : (int)((x) + 0.5))
#define ABS(x) ((x) < 0 ? -(x) : (x))
#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SIGNOF(x) ((x) >= 0.0 ? 1 : -1)

#endif  //_MACROS_H_
