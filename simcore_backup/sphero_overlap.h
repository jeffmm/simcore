#ifndef _SIMCORE_SPHERO_OVERLAP_H_
#define _SIMCORE_SPHERO_OVERLAP_H_

int sphero_overlap(int n_dim, int n_periodic, double **h, double skin, double diameter_12,
                   double *r_1, double *s_1, double *u_1, double length_1,
                   double *r_2, double *s_2, double *u_2, double length_2);

#endif
