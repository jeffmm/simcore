#ifndef _SIMCORE_WRITE_BMP_H_
#define _SIMCORE_WRITE_BMP_H_

#include <stdint.h>

#define ERROR 1
#define SUCCESS 0

/* WARNING : THESE DECLARATIONS ARE ONLY FOR 64bit PROCESSORS */
struct template_bmpheader {
    uint32_t bfsize;
    uint16_t reserved1;
    uint16_t reserved2;
    uint32_t headersize;
};

struct template_bmpfile {
    uint32_t bmpsize;
    int32_t bmpwidth;
    int32_t bmpheight;
    uint16_t biplanes;
    uint16_t bitcount;
    uint32_t bicompression;
    uint32_t bisizeimage;
    int32_t bixpelspermeter;
    int32_t biypelspermeter;
    uint32_t biclrused;
    uint32_t biclrimportant;
};

#endif
