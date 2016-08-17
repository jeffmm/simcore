// TODO: There are certainly much better ways of doing this. Object container and object params structures most likely, for more clear inputs.
#ifdef _SIMCORE_BR_DIMER_H_
if (params_.n_dimer > 0) {
  SpeciesBase * spcs = new BrDimerSpecies(params_.n_dimer, &params_, space_.GetStruct(), gsl_rng_get(rng_.r));
  spcs->Init();
  species_.push_back(spcs);
}
#endif
#ifdef _SIMCORE_BR_BEAD_H_
if (params_.n_br_bead > 0) {
  SpeciesBase * spcs = new BrBeadSpecies(params_.n_br_bead, &params_, space_.GetStruct(), gsl_rng_get(rng_.r));
  spcs->Init();
  species_.push_back(spcs);
}
#endif
#ifdef _SIMCORE_NEON_H_
if (params_.n_neon > 0) {
  SpeciesBase * spcs = new NeonSpecies(params_.n_neon, &params_, space_.GetStruct(), gsl_rng_get(rng_.r));
  spcs->Init();
  species_.push_back(spcs);
}
#endif
#ifdef _SIMCORE_ARGON_H_
if (params_.n_argon > 0) {
  SpeciesBase * spcs = new ArgonSpecies(params_.n_argon, &params_, space_.GetStruct(), gsl_rng_get(rng_.r));
  spcs->Init();
  species_.push_back(spcs);
}
#endif
#ifdef _SIMCORE_MD_BEAD_H_
if (params_.n_md_bead > 0) {
  SpeciesBase * spcs = new MDBeadSpecies(params_.n_md_bead, &params_, space_.GetStruct(), gsl_rng_get(rng_.r));
  spcs->Init();
  species_.push_back(spcs);
}
#endif
#ifdef _SIMCORE_BR_ROD_H_
if (params_.n_rod > 0) {
  SpeciesBase * spcs = new BrRodSpecies(params_.n_rod, &params_, space_.GetStruct(), gsl_rng_get(rng_.r));
  //spcs->Init();
  // Special stuff for a rod
  BrRodSpecies* prspec = dynamic_cast<BrRodSpecies*>(spcs);
  if (strcmp(params_.config_file, "random") == 0) {
    prspec->Init();
  } else {
    prspec->ConfiguratorRod();
  }
  species_.push_back(spcs);
}
#endif
#ifdef _SIMCORE_MD_KMC_BEAD_H_
if (params_.n_md_kmc_bead > 0) {
  SpeciesBase * spcs = new MDKMCBeadSpecies(params_.n_md_kmc_bead, &params_, space_.GetStruct(), gsl_rng_get(rng_.r));
  spcs->Init();
  species_.push_back(spcs);
}
#endif
#ifdef _SIMCORE_BR_IMPLICIT_H_
if (params_.n_br_implicit > 0) {
  SpeciesBase * spcs = new BrImplicitBeadSpecies(params_.n_br_implicit, &params_, space_.GetStruct(), gsl_rng_get(rng_.r));
  spcs->Init();
  species_.push_back(spcs);
}
#endif
#ifdef _SIMCORE_BR_WALKER_H_
if (params_.n_br_walker > 0) {
  SpeciesBase * spcs = new BrWalkerSpecies(params_.n_br_walker, &params_, space_.GetStruct(), gsl_rng_get(rng_.r));
  spcs->Init();
  species_.push_back(spcs);
}
#endif
#ifdef _SIMCORE_XLINK_H_
if (params_.n_xlink > 0) {
  SpeciesBase * spcs = new XlinkSpecies(params_.n_xlink, &params_, space_.GetStruct(), gsl_rng_get(rng_.r));
  XlinkSpecies *pxspec = dynamic_cast<XlinkSpecies*>(spcs);
  if (strcmp(params_.config_file, "random") == 0) {
    pxspec->Init();
  } else {
    pxspec->ConfiguratorXlink();
  }
  species_.push_back(spcs);
}
#endif
#ifdef _SIMCORE_FILAMENT_H_
if (params_.n_filament > 0) {
  SpeciesBase * spcs = new FilamentSpecies(params_.n_filament, &params_, space_.GetStruct(), gsl_rng_get(rng_.r));
  spcs->Init();
  species_.push_back(spcs);
}
#endif
