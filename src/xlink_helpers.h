#ifndef _SIMCORE_XLINK_HELPERS_H_
#define _SIMCORE_XLINK_HELPERS_H_

#include <gsl/gsl_integration.h>

// Helper routines for xlinks
namespace xlh {

  struct _xlink_params {
    double alpha;
    double r0;
    double y0;
  };
  typedef struct _xlink_params xlink_params;

  double integrand_1_2(double x, void *params) {
    xlink_params *myparams = (xlink_params*)params;
    double alpha = myparams->alpha;
    double r0 = myparams->r0;
    double y0 = myparams->y0;

    double exponent = sqrt(x*x + y0*y0) - r0;
    exponent *= -alpha * exponent;
    return exp(exponent);
  }

  double prob_1_2(std::vector<double> &x, void *params) {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

    xlink_params *inparams = (xlink_params*)params;
    xlink_params gslparams;
    gslparams.alpha = inparams->alpha;
    gslparams.r0 = inparams->r0;
    gslparams.y0 = x[1];

    gsl_function F;
    F.function = &integrand_1_2;
    F.params = &gslparams;

    double result, error;
    if (x[0] > 0)
      gsl_integration_qags(&F, 0, x[0], 0, 1e-7, 1000,
          w, &result, &error);
    else
      result = 0;

    gsl_integration_workspace_free(w);

    return result;
  }
}

#endif
