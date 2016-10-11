/* This routine checks for overlap between two spherocylinders, which could be nonadditive
   (a different overlap criterion can be applied to every distinct pair of types).
   If skin > 0.0, "near" overlaps are also detected (this feature is used
   in constructing the neighbor lists).

input: number of spatial dimensions (n_dim)
       number of periodic dimensions (n_periodic)
       unit cell matrix (h)
       neighbor list skin, if applicable (skin)
       overlap distance (diameter_12)
       real position of first spherocylinder (r_1)
       scaled position of first spherocylinder (s_1)
       director of first spherocylinder (u_1)
       length of first spherocylinder (length_1)
       real position of second spherocylinder (r_2)
       scaled position of second spherocylinder (s_2)
       director of second spherocylinder (u_2)
       length of second spherocylinder (length_2)

output: flag indicating whether spherocylinders are overlapping;
        1 if so, 0 if not (return value). */

#define SMALL 1.0e-12

#include "auxiliary.h"

int sphero_overlap(int n_dim, int n_periodic, double **h, double skin, double diameter_12,
                   double *r_1, double *s_1, double *u_1, double length_1,
                   double *r_2, double *s_2, double *u_2, double length_2)
{
    int i, j;
    double half_length_1, half_length_2,
        dr2, dr_dot_u_1, dr_dot_u_2, u_1_dot_u_2, denom,
        lambda, mu, lambda_a, lambda_b, mu_a, mu_b, lambda_mag, 
        mu_mag, r_min2, r_min2_a, r_min2_b;
    double ds[3], dr[3], r_min[3], r_min_a[3], r_min_b[3];
    
    /* Compute various constants. */
    half_length_1 = 0.5 * length_1;
    half_length_2 = 0.5 * length_2;

    /* Compute pair separation vector. */
    for (i = 0; i < n_periodic; ++i) {  /* First handle periodic subspace. */
        ds[i] = s_2[i] - s_1[i];
        ds[i] -= NINT(ds[i]);
    }
    for (i = 0; i < n_periodic; ++i) {
        dr[i] = 0.0;
        for (j = 0; j < n_periodic; ++j)
            dr[i] += h[i][j] * ds[j];
    }
    for (i = n_periodic; i < n_dim; ++i)        /* Then handle free subspace. */
        dr[i] = r_2[i] - r_1[i];

    /* Compute squared pair separation. */
    dr2 = 0.0;
    for (i = 0; i < n_dim; ++i)
        dr2 += SQR(dr[i]);

    /* If the pair separation exceeds the spherocylinder length plus 
       diameter plus the skin thickness, the spherocylinders don't overlap. */
    if (dr2 > SQR(half_length_1 + half_length_2 + diameter_12 + skin))
        return 0;

    /* Compute minimum distance (see Allen et al., Adv. Chem. Phys. 86, 1 (1993)).
       First consider two infinitely long lines. */
    dr_dot_u_1 = dr_dot_u_2 = u_1_dot_u_2 = 0.0;
    for (i = 0; i < n_dim; ++i) {
        dr_dot_u_1 += dr[i] * u_1[i];
        dr_dot_u_2 += dr[i] * u_2[i];
        u_1_dot_u_2 += u_1[i] * u_2[i];
    }
    denom = 1.0 - SQR(u_1_dot_u_2);
    if (denom < SMALL) {
        lambda = dr_dot_u_1 / 2.0;
        mu = -dr_dot_u_2 / 2.0;
    } else {
        lambda = (dr_dot_u_1 - u_1_dot_u_2 * dr_dot_u_2) / denom;
        mu = (-dr_dot_u_2 + u_1_dot_u_2 * dr_dot_u_1) / denom;
    }
    lambda_mag = ABS(lambda);
    mu_mag = ABS(mu);

    /* Now take into account the fact that the two line segments are of finite length. */
    if (lambda_mag > half_length_1 && mu_mag > half_length_2) {

        /* Calculate first possible case. */
        lambda_a = SIGN(half_length_1, lambda);
        mu_a = -dr_dot_u_2 + lambda_a * u_1_dot_u_2;
        mu_mag = ABS(mu_a);
        if (mu_mag > half_length_2)
            mu_a = SIGN(half_length_2, mu_a);

        /* Calculate minimum distance between two spherocylinders. */
        r_min2_a = 0.0;
        for (i = 0; i < n_dim; ++i) {
            r_min_a[i] = dr[i] - lambda_a * u_1[i] + mu_a * u_2[i];
            r_min2_a += SQR(r_min_a[i]);
        }

        /* Calculate second possible case. */
        mu_b = SIGN(half_length_2, mu);
        lambda_b = dr_dot_u_1 + mu_b * u_1_dot_u_2;
        lambda_mag = ABS(lambda_b);
        if (lambda_mag > half_length_1)
            lambda_b = SIGN(half_length_1, lambda_b);

        /* Calculate minimum distance between two spherocylinders. */
        r_min2_b = 0.0;
        for (i = 0; i < n_dim; ++i) {
            r_min_b[i] = dr[i] - lambda_b * u_1[i] + mu_b * u_2[i];
            r_min2_b += SQR(r_min_b[i]);
        }

        /* Choose the minimum minimum distance. */
        r_min2 = MIN(r_min2_a, r_min2_b);
    } else if (lambda_mag > half_length_1) {

        /* Adjust lambda and mu. */
        lambda = SIGN(half_length_1, lambda);
        mu = -dr_dot_u_2 + lambda * u_1_dot_u_2;
        mu_mag = ABS(mu);
        if (mu_mag > half_length_2)
            mu = SIGN(half_length_2, mu);

        /* Calculate minimum distance between two spherocylinders. */
        r_min2 = 0.0;
        for (i = 0; i < n_dim; ++i) {
            r_min[i] = dr[i] - lambda * u_1[i] + mu * u_2[i];
            r_min2 += SQR(r_min[i]);
        }
    } else if (mu_mag > half_length_2) {

        /* Adjust lambda and mu. */
        mu = SIGN(half_length_2, mu);
        lambda = dr_dot_u_1 + mu * u_1_dot_u_2;
        lambda_mag = ABS(lambda);
        if (lambda_mag > half_length_1)
            lambda = SIGN(half_length_1, lambda);

        /* Calculate minimum distance between two spherocylinders. */
        r_min2 = 0.0;
        for (i = 0; i < n_dim; ++i) {
            r_min[i] = dr[i] - lambda * u_1[i] + mu * u_2[i];
            r_min2 += SQR(r_min[i]);
        }
    } else {

        /* Calculate minimum distance between two spherocylinders. */
        r_min2 = 0.0;
        for (i = 0; i < n_dim; ++i) {
            r_min[i] = dr[i] - lambda * u_1[i] + mu * u_2[i];
            r_min2 += SQR(r_min[i]);
        }
    }

    /* Return 1 if spherocylinders overlap, 0 if not. */
    return (r_min2 <= SQR(diameter_12 + skin));
}

#undef SMALL
