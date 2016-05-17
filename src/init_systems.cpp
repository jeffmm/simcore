#ifdef _CYTOSCORE_BROWNIAN_DIMER_H_
for (int i=0; i<params_.n_dimer; ++i) {
  BrownianDimer bdmr(params_.n_dim, params_.delta, gsl_rng_get(rng.r));
  bdmr.Init(&params_);
  dimers_.push_back(bdmr);
}
#endif
