#include "simcore/auxiliary.hpp"
#include "simcore/graphics.hpp"

#ifndef NOGRAPH

/* Functions in writegif.c */
/* extern void writegif(FILE *fp, XImage *image, int width, int height, int
 * numcol, unsigned char rgb[][256]); */

/* Functions in writebmp.c */
extern void write_bmp24(FILE *fp, int width, int height, GLvoid *pixels);

/* Local Functions */
int getelem2d(int i, int j, int k, int width, int height, int depth,
              GLubyte *bmpmem);
void putelem2d(int fwdvalue, int i, int j, int k, int width, int height,
               int depth, GLubyte *array3d);

void grabber(int width, int height, std::string fname, int framenum) {
  int pixelsqty;
  int i, j, k;
  int fwdvalue;
  char bmpfname[2048];
  FILE *fp_bmp;
  GLubyte *pixels;
  GLubyte *revpixels;
  int length, zeronum;
  char sframenum[6], szero[6];

  /* Determine pixels quantity */
  pixelsqty = 3 * width * height;

  /* Set up Alignment format for bitmaps */
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);

  /* Define pixel array for a no extra byte. Mode 1. */
  pixels = (GLubyte *)malloc(pixelsqty);
  if (pixels == NULL) {
    Logger::Error("Not enough memory in grabber");
  }

  /* Define revpixel array for a no extra byte. Mode 1. */
  revpixels = (GLubyte *)malloc(pixelsqty);
  if (revpixels == NULL) {
    Logger::Error("Not enough memory in grabber");
  }

  /* Store Full Image from the back buffer */
  glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);

  /* Store image in memory */
  glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels);

  /* Save pixels' RGB values inverse to an array */
  for (i = 0; i < width; i++)
    for (j = 0; j < height; j++)
      for (k = 0; k < 3; k++) {
        fwdvalue = getelem2d(i, j, k, width, height, 3, pixels);
        putelem2d(fwdvalue, i, j, (2 - k), width, height, 3, revpixels);
      }

  /* Create 24-Bit bitmap filename. */
  sprintf(sframenum, "%d", framenum);
  length = strlen(sframenum);
  zeronum = 5 - length;
  strcpy(szero, "");
  for (i = 0; i < zeronum; ++i) strcat(szero, "0");
  sprintf(bmpfname, "%s_%s%s.bmp", fname.c_str(), szero, sframenum);
  // printf("%s\n", bmpfname);

  /* Open file */
  fp_bmp = fopen(bmpfname, "wb");
  if (fp_bmp == NULL) {
    Logger::Error("Unable to open output file in grabber");
  }

  /* Write 24-bit BMP file */
  write_bmp24(fp_bmp, width, height, revpixels);

  /* Close file and free memory */
  fclose(fp_bmp);
  free(pixels);
  free(revpixels);
}

int getelem2d(int i, int j, int k, int width, int height, int depth,
              GLubyte *array3d) {
  int value;

  value = array3d[k + depth * (j + height * i)];

  return value;
}

void putelem2d(int fwdvalue, int i, int j, int k, int width, int height,
               int depth, GLubyte *array3d) {
  array3d[k + depth * (j + height * i)] = fwdvalue;
}

#endif
