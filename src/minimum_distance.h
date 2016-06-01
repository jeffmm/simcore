#ifndef _MINIMUM_DISTANCE_H
#define _MINIMUM_DISTANCE_H

void min_distance_point_point(int n_dim, int n_periodic, double **unit_cell, 
                              double const * const r1, double const * const s1, 
                              double const * const r2, double const * const s2, 
                              double *dr, double *dr_mag2);

void min_distance_point_carrier_line(int n_dim, int n_periodic, double **h,
                                     double *r_point, double *s_point,
                                     double *r_line, double *s_line, double *u_line,
                                     double length, double *dr, double *r_contact);

void min_distance_sphero(int n_dim, int n_periodic, double **h,
                         double const * const r_1, double const * const s_1,
                         double const * const u_1, double const length_1,
                         double const * const r_2, double const * const s_2, 
                         double const * const u_2, double const length_2,
                         double *r_min, double *r_min_mag2, 
                         double *contact1, double *contact2);

void min_distance_sphero_dr(int n_dim, int n_periodic, double **h,
                            double *r_1, double *s_1, double *u_1, double length_1,
                            double *r_2, double *s_2, double *u_2, double length_2,
                            double *dr, double *r_min, double *r_min_mag2, 
                            double *lambda, double *mu);

void min_distance_sphere_sphero(int n_dim, int n_periodic, double **h,
                                double const * const r_1, double const * const s_1,
                                double const * const r_2, double const * const s_2, 
                                double const * const u_2, double const length_2,
                                double *r_min, double *r_min_mag2,
                                double *contact2);

void min_distance_sphero_plane(double *r_bond, double *u_bond, double length,
                               double *r_plane, double *n_plane, double *lambda,
                               double *r_min_mag2, double *r_min);


void min_distance_carrier_lines(int n_dim, int n_periodic, double **h,
                                double *r_1, double *s_1, double *u_1, 
                                double *r_2, double *s_2, double *u_2, 
                                double *r_min, double *r_min_mag2, 
                                double *lambda, double *mu);

#endif
