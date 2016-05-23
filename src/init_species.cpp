// TODO: There are certainly much better ways of doing this. Object container and object params structures most likely, for more clear inputs.
#ifdef _CYTOSCORE_BROWNIAN_DIMER_H_
if (params_.n_dimer > 0) {
  SpeciesBase * spcs = new BrownianDimerSpecies(params_.n_dimer, &params_, space_.GetStruct(), gsl_rng_get(rng_.r));
  spcs->Init();
  species_.push_back(spcs);
}
#endif
#ifdef _CYTOSCORE_BROWNIAN_BEAD_H_
if (params_.n_br_bead > 0) {
  SpeciesBase * spcs = new BrownianBeadSpecies(params_.n_br_bead, &params_, space_.GetStruct(), gsl_rng_get(rng_.r));
  spcs->Init();
  species_.push_back(spcs);
}
#endif
#ifdef _CYTOSCORE_MD_BEAD_H_
if (params_.n_md_bead > 0) {
  SpeciesBase * spcs = new MDBeadSpecies(params_.n_md_bead, &params_, space_.GetStruct(), gsl_rng_get(rng_.r));
  spcs->Init();
  species_.push_back(spcs);
}
#endif
