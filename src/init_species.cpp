#ifdef _CYTOSCORE_BROWNIAN_DIMER_H_
SpeciesBase * spcs = new Species<BrownianDimer>(params_.n_dimer, &params_, space_.GetStruct(), gsl_rng_get(rng_.r));
spcs->Init();
species_.push_back(spcs);
#endif
