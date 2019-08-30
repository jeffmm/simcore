#ifndef NOGRAPH

#include "writebmp.hpp"
#include "graphics.hpp"

/* Local Functions */
void write_bmp24(FILE *fp, int width, int height, GLvoid *pixels);
void write_bmpheader(FILE *fp, int width, int height);
void write_bmpimage(FILE *fp, int width, int height, GLvoid *pixels);
void checkvalues(int width, int height, GLvoid *pixels);

/* Global Variables */
int pad, rasterwidth;
struct template_bmpheader bmpheader;
struct template_bmpfile bmpfile;

void write_bmp24(FILE *fp, int width, int height, GLvoid *pixels) {
  /* Write header */
  write_bmpheader(fp, width, height);

  /* Write image data */
  write_bmpimage(fp, width, height, pixels);
}

/* Support Functions ------------------------------------------------------ */

void write_bmpheader(FILE *fp, int width, int height) {
  int bmpnumbytes;

  /* Initialize header structure */
  bmpheader.reserved1 = 0;
  bmpheader.reserved2 = 0;

  /* BMP header defined to be 54 bytes by spec */
  bmpheader.headersize = 54;

  /* Determine image width if it is not a multiple of 4. This is, in
     24-bit bitmap files each raster line should be a multiple of 4 bytes */
  if ((width % 4) == 0)
    pad = 0;
  else
    pad = 4 - (width % 4);

  /* Padded new width */
  rasterwidth = width + pad;

  /* Number of bytes in image */
  bmpnumbytes = 3 * rasterwidth * height;

  /* File size in bytes */
  bmpheader.bfsize = 54 + bmpnumbytes;

  /* Start writing the file. First 2 bytes always the chars B and M for a bmp
   * file */
  fputc('B', fp);
  fputc('M', fp);

  /* Write bitmap header info. First section is 12 bytes */
  fwrite(&bmpheader, 12, 1, fp);

  /* Bitmap Info. sizeof file header is always 40 bytes */
  bmpfile.bmpsize = 40;
  bmpfile.bmpwidth = rasterwidth;
  bmpfile.bmpheight = height;
  bmpfile.biplanes = 1;
  bmpfile.bitcount = 24;
  bmpfile.bicompression = 0;
  bmpfile.bisizeimage = 0;
  bmpfile.bixpelspermeter = 0;
  bmpfile.biypelspermeter = 0;
  bmpfile.biclrused = 0;
  bmpfile.biclrimportant = 0;

  /* Write second header. File pointer should now be at pixel array */
  fwrite(&bmpfile, 40, 1, fp);
}

void write_bmpimage(FILE *fp, int width, int height, GLvoid *pixels) {
  int i, m, n;
  unsigned char zero = 0;
  unsigned char *curpix;

  /* Assign pointer to memory image */
  curpix = (unsigned char *)pixels;

  /* Save image */
  if (pad == 0) {
    fwrite(curpix, 1, 3 * width * height, fp);
  } else {
    for (i = 0; i < height; i++) {
      /* for (j = 0; j < width; j++) */
      /*     /\* Position pointer to the end of next pixels *\/ */
      /*     for (k = 0; k < 3; k++) */
      /*             fwrite(curpix++, 1, 3, fp); */

      /* Write horizontal line to output file */
      fwrite(curpix, 1, 3 * width, fp);
      curpix += 3 * width;

      /* Fill out remaining zero */
      if (pad > 0)
        for (m = 0; m < pad; m++)
          for (n = 0; n < 3; n++) fwrite(&zero, 1, 1, fp);
    }
  }
}

/* Miscellaneous Functions ------------------------------------------------ */

void checkvalues(int width, int height, GLvoid *pixels) {
  int i, j, numpixels;
  unsigned char *curpix;

  curpix = (unsigned char *)pixels;

  numpixels = width * height;

  for (i = 0; i < numpixels; i++) {
    printf("pixel num : %d \n", i);

    for (j = 0; j < 3; j++) printf("%d \n", *curpix++);

    printf("\n");
  }
}

#endif
