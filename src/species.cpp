#include "species.h"

void SpeciesBase::Init(system_parameters *params, space_struct *space, long seed) {
  params_ = params;
  sparams_ = &(params->species);
  space_ = space;
  rng_.Init(seed);
}

void SpeciesBase::InitPositFile(std::string run_name) {
  std::string sid_str = sid_._to_string();
  std::string posit_file_name = run_name + "_" + sid_str + ".posit";
  oposit_file_.open(posit_file_name, std::ios::out | std::ios::binary ); 
  if (!oposit_file_.is_open()) {
    std::cout<<"ERROR: Output file "<< posit_file_name <<" did not open\n";
    exit(1);
  }
  oposit_file_.write(reinterpret_cast<char*> (&params_->n_steps), sizeof(int));
  oposit_file_.write(reinterpret_cast<char*> (&sparams_->n_posit), sizeof(int));
  oposit_file_.write(reinterpret_cast<char*> (&params_->delta), sizeof(double));
}

void SpeciesBase::InitPositFileInput(std::string run_name) {
  std::string sid_str = sid_._to_string();
  std::string posit_file_name = run_name + "_" + sid_str + ".posit";
  iposit_file_.open(posit_file_name, std::ios::in | std::ios::binary ); 
  if (!iposit_file_.is_open()) {
    std::cout<<"ERROR: Input file "<< posit_file_name <<" did not open\n";
    exit(1);
  }
  int n_steps, n_posit;
  double delta;
  iposit_file_.read(reinterpret_cast<char*> (&n_steps), sizeof(int));
  iposit_file_.read(reinterpret_cast<char*> (&n_posit), sizeof(int));
  iposit_file_.read(reinterpret_cast<char*> (&delta), sizeof(double));
  if (n_steps != params_->n_steps || n_posit != sparams_->n_posit || delta != params_->delta) {
    std::cout << "ERROR: Input file " << posit_file_name << " does not match parameter file\n";
    printf("%d %d, %d %d, %2.5f %2.5f\n",n_steps, params_->n_steps, n_posit, sparams_->n_posit, delta, params_->delta);
    exit(1);
  }
  ReadPosits();
}


void SpeciesBase::InitSpecFile(std::string run_name) {
  std::string sid_str = sid_._to_string();
  std::string spec_file_name = run_name + "_" + sid_str + ".spec";
  ospec_file_.open(spec_file_name, std::ios::out | std::ios::binary ); 
  if (!ospec_file_.is_open()) {
    std::cout<<"ERROR: Output file "<< spec_file_name <<" did not open\n";
    exit(1);
  }
  ospec_file_.write(reinterpret_cast<char*> (&params_->n_steps), sizeof(int));
  ospec_file_.write(reinterpret_cast<char*> (&sparams_->n_spec), sizeof(int));
  ospec_file_.write(reinterpret_cast<char*> (&params_->delta), sizeof(double));
}

void SpeciesBase::InitSpecFileInput(std::string run_name) {
  std::string sid_str = sid_._to_string();
  std::string spec_file_name = run_name + "_" + sid_str + ".spec";
  ispec_file_.open(spec_file_name, std::ios::in | std::ios::binary ); 
  if (!ispec_file_.is_open()) {
    std::cout<<"ERROR: Input file "<< spec_file_name <<" did not open\n";
    exit(1);
  }
  int n_steps, n_spec;
  double delta;
  ispec_file_.read(reinterpret_cast<char*> (&n_steps), sizeof(int));
  ispec_file_.read(reinterpret_cast<char*> (&n_spec), sizeof(int));
  ispec_file_.read(reinterpret_cast<char*> (&delta), sizeof(double));
  if (n_steps != params_->n_steps || n_spec != sparams_->n_spec || delta != params_->delta) {
    std::cout<< "ERROR: Input file " << spec_file_name << " does not match parameter file\n";
    exit(1);
  }
  ReadSpecs();
}

void SpeciesBase::InitOutputFiles(std::string run_name) {
  if (sparams_->posit_flag) 
    InitPositFile(run_name);
  if (sparams_->spec_flag) 
    InitSpecFile(run_name);
  if (sparams_->checkpoint_flag) {
    std::string sid_str = sid_._to_string();
    checkpoint_file_ = run_name + "_" + sid_str + ".checkpoint";
  }
}

void SpeciesBase::InitCheckpoints(std::string run_name) {
  std::string sid_str = sid_._to_string();
  checkpoint_file_ = run_name + "_" + sid_str + ".checkpoint";
  if (!sparams_->checkpoint_flag) {
    std::cout << "ERROR: Checkpoint file " << checkpoint_file_ << " not available for parameter file!\n";
    exit(1);
  }
  ReadCheckpoints();
  std::cout << "WARNING: Loading checkpoint. SimCORE will overwrite any present output files from previous simulation \"" << run_name << "\". Proceed? (y/N)\n";
  std::string response;
  std::cin >> response;
  if (response.compare("y") != 0 && response.compare("Y") != 0)
    exit(0);
  // XXX Need to fix appending to output file in the correct position
  if (sparams_->posit_flag) {
    std::string posit_file_name = run_name + "_" + sid_str + ".posit";
    oposit_file_.open(posit_file_name, std::ios::out | std::ios::binary ); 
    if (!oposit_file_.is_open()) {
      std::cout<<"ERROR: Output file "<< posit_file_name <<" did not open for appending\n";
      exit(1);
    }
    oposit_file_.write(reinterpret_cast<char*> (&params_->n_steps), sizeof(int));
    oposit_file_.write(reinterpret_cast<char*> (&sparams_->n_posit), sizeof(int));
    oposit_file_.write(reinterpret_cast<char*> (&params_->delta), sizeof(double));
  }
  if (sparams_->spec_flag) {
    std::string spec_file_name = run_name + "_" + sid_str + ".spec";
    ospec_file_.open(spec_file_name, std::ios::out | std::ios::binary ); 
    if (!ospec_file_.is_open()) {
      std::cout<<"ERROR: Output file "<< spec_file_name <<" did not open for appending\n";
      exit(1);
    }
    ospec_file_.write(reinterpret_cast<char*> (&params_->n_steps), sizeof(int));
    ospec_file_.write(reinterpret_cast<char*> (&sparams_->n_spec), sizeof(int));
    ospec_file_.write(reinterpret_cast<char*> (&params_->delta), sizeof(double));
  }
}

void SpeciesBase::InitInputFiles(std::string run_name, bool posits_only) {
  if (sparams_->posit_flag) 
    InitPositFileInput(run_name);
  if (!posits_only && sparams_->spec_flag) 
    InitSpecFileInput(run_name);
}

void SpeciesBase::CloseFiles() { 
  if (oposit_file_.is_open())
    oposit_file_.close(); 
  if (iposit_file_.is_open())
    iposit_file_.close(); 
  if (ospec_file_.is_open())
    ospec_file_.close(); 
  if (ispec_file_.is_open())
    ispec_file_.close(); 
  FinalizeAnalysis();
}

std::vector<Object*> SpeciesBase::GetInteractors() {
  std::vector<Object*> ix;
  return ix;
}
