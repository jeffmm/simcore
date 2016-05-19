#ifdef _CYTOSCORE_BROWNIAN_DIMER_H_
for (int i=0; i<params_.n_dimer; ++i) {
  SpeciesBase * bdmr = new BrownianDimer(params_.n_dim, params_.delta, gsl_rng_get(rng.r));
  bdmr->Init(&params_, space.GetStruct());
  species_.push_back(bdmr);
}
#endif
