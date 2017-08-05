#include "species.h"

void Species::Init(system_parameters *params, space_struct *space, long seed) {
  spec_name_ = "";
  params_ = params;
  sparams_ = &(params->species);
  space_ = space;
  rng_.Init(seed);
}

void Species::InitPositFile(std::string run_name) {
  std::string posit_file_name = run_name + "_" + spec_name_ + ".posit";
  oposit_file_.open(posit_file_name, std::ios::out | std::ios::binary ); 
  if (!oposit_file_.is_open()) {
    std::cout<<"ERROR: Output file "<< posit_file_name <<" did not open\n";
    exit(1);
  }
  oposit_file_.write(reinterpret_cast<char*> (&params_->n_steps), sizeof(int));
  oposit_file_.write(reinterpret_cast<char*> (&sparams_->n_posit), sizeof(int));
  oposit_file_.write(reinterpret_cast<char*> (&params_->delta), sizeof(double));
}

void Species::InitPositFileInput(std::string run_name) {
  std::string posit_file_name = run_name + "_" + spec_name_ + ".posit";
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


void Species::InitSpecFile(std::string run_name) {
  std::string spec_file_name = run_name + "_" + spec_name_ + ".spec";
  ospec_file_.open(spec_file_name, std::ios::out | std::ios::binary ); 
  if (!ospec_file_.is_open()) {
    std::cout<<"ERROR: Output file "<< spec_file_name <<" did not open\n";
    exit(1);
  }
  ospec_file_.write(reinterpret_cast<char*> (&params_->n_steps), sizeof(int));
  ospec_file_.write(reinterpret_cast<char*> (&sparams_->n_spec), sizeof(int));
  ospec_file_.write(reinterpret_cast<char*> (&params_->delta), sizeof(double));
}

void Species::InitSpecFileInput(std::string run_name) {
  std::string spec_file_name = run_name + "_" + spec_name_ + ".spec";
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

void Species::InitOutputFiles(std::string run_name) {
  if (sparams_->posit_flag) 
    InitPositFile(run_name);
  if (sparams_->spec_flag) 
    InitSpecFile(run_name);
  if (sparams_->checkpoint_flag) {
    checkpoint_file_ = run_name + "_" + spec_name_ + ".checkpoint";
  }
}

void Species::InitCheckpoints(std::string run_name) {
  checkpoint_file_ = run_name + "_" + spec_name_ + ".checkpoint";
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
    std::string posit_file_name = run_name + "_" + spec_name_ + ".posit";
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
    std::string spec_file_name = run_name + "_" + spec_name_ + ".spec";
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

void Species::InitInputFiles(std::string run_name, bool posits_only) {
  if (sparams_->posit_flag) 
    InitPositFileInput(run_name);
  if (!posits_only && sparams_->spec_flag) 
    InitSpecFileInput(run_name);
}

void Species::CloseFiles() { 
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

void Species::AddMember() {
  Object* newmember = new Object;
  newmember->Init();
  members_.push_back(newmember);
  n_members_++;
}
 
void Species::AddMember(Object* newmem) {
  members_.push_back(newmem);
  n_members_++;
}
 
void Species::PopMember() {
  delete members_[n_members_-1];
  members_.pop_back();
  n_members_--;
}

void Species::PopAll() {
  while (n_members_ > 0) {
    PopMember();
  }
}
 
void Species::Draw(std::vector<graph_struct*> * graph_array) {
  for (auto it=members_.begin(); it!=members_.end(); ++it)
    (*it)->Draw(graph_array);
}

void Species::UpdatePositions() {
  for (auto it=members_.begin(); it!=members_.end(); ++it)
    (*it)->UpdatePosition();
}

std::vector<Object*> Species::GetInteractors() {
  std::vector<Object*> interactors;
  for (auto it=members_.begin(); it!=members_.end(); ++it) {
    std::vector<Object*> ix_vec = (*it)->GetInteractors();
    interactors.insert(interactors.end(), ix_vec.begin(), ix_vec.end());
  }
  return interactors;
}

double Species::GetKineticEnergy() {
  double ke=0;
  for (auto it=members_.begin(); it!=members_.end(); ++it)
    ke+=(*it)->GetKineticEnergy();
  return ke;
}

double Species::GetPotentialEnergy() {
  double pe=0;
  for (auto it=members_.begin(); it!=members_.end(); ++it)
    pe+=(*it)->GetPotentialEnergy();
    /* The total potential energy is going to be half of the
     potential energy felt by each particle. Potential energy 
     is shared, so I need to avoid double counting. */
  return 0.5*pe;
}

void Species::ZeroForces() {
  for (auto it=members_.begin(); it!=members_.end(); ++it) 
    (*it)->ZeroForce();
}

void Species::Report() {
  for (auto it=members_.begin(); it!=members_.end(); ++it) 
    (*it)->Report();
}

int Species::GetCount() {
  int count = 0;
  for (auto it=members_.begin(); it!=members_.end(); ++it) 
    count += (*it)->GetCount();
  return count;
}
 
void Species::WritePosits() {
  int size = members_.size();
  oposit_file_.write(reinterpret_cast<char*>(&size), sizeof(size));
  for (auto it=members_.begin(); it!=members_.end(); ++it) 
    (*it)->WritePosit(oposit_file_);
}
 
void Species::WriteSpecs() {
  int size = members_.size();
  ospec_file_.write(reinterpret_cast<char*>(&size), sizeof(size));
  for (auto it=members_.begin(); it!=members_.end(); ++it) 
    (*it)->WriteSpec(ospec_file_);
}

void Species::WriteCheckpoints() {
  int size = members_.size();
  std::fstream ocheck_file(checkpoint_file_,std::ios::out | std::ios::binary);
  if (!ocheck_file.is_open()) {
    std::cout<<"ERROR: Output "<< checkpoint_file_ <<" file did not open\n";
    exit(1);
  }
  ocheck_file.write(reinterpret_cast<char*>(&size), sizeof(size));
  for (auto it=members_.begin(); it!=members_.end(); ++it) 
    (*it)->WriteCheckpoint(ocheck_file);
  ocheck_file.close();
}

void Species::ReadPosits() {
  if (iposit_file_.eof()) {
    printf("  EOF reached\n");
    early_exit = true;
    return;
  }
  if (! iposit_file_.is_open()) {
    printf(" ERROR: Posit file unexpectedly not open! Exiting.\n");
    early_exit=true;
    return;
  }
  int size = -1;
  Object *member;
  iposit_file_.read(reinterpret_cast<char*>(&size), sizeof(size));
  // Hacky workaround FIXME
  // This prevents strange errors that occasionally crop up when reading inputs
  // from early exits
  if (size < 0) {
    early_exit = true;
    return;
  }
  if (size != n_members_) {
    member = new Object;
    members_.resize(size, member);
    delete member;
    n_members_ = size;
  }
  for (auto it=members_.begin(); it!=members_.end(); ++it) 
    (*it)->ReadPosit(iposit_file_);
}
 
void Species::ReadCheckpoints() {
  std::fstream icheck_file(checkpoint_file_, std::ios::in | std::ios::binary);
  if (!icheck_file.is_open()) {
    std::cout<<"  ERROR: Output "<< checkpoint_file_ <<" file did not open\n";
    exit(1);
  }
  int size = 0;
  Object *member;
  icheck_file.read(reinterpret_cast<char*>(&size), sizeof(size));
  member = new Object;
  members_.resize(size, member);
  for (auto it=members_.begin(); it!=members_.end(); ++it) 
    (*it)->ReadCheckpoint(icheck_file);
  icheck_file.close();
}
 
void Species::ReadSpecs() {
  if (ispec_file_.eof()) {
    printf("  EOF reached\n");
    early_exit = true;
    return;
  }
  if (! ispec_file_.is_open()) {
    printf(" ERROR: Spec file unexpectedly not open! Exiting.\n");
    early_exit = true;
    return;
  }
  int size = -1;
  Object *member;
  ispec_file_.read(reinterpret_cast<char*>(&size), sizeof(size));
  // Hacky workaround FIXME
  if (size < 0) {
    early_exit = true;
    return;
  }
  if (size != n_members_) {
    member = new Object;
    members_.resize(size, member);
    delete member;
    n_members_ = size;
  }
  for (auto it=members_.begin(); it!=members_.end(); ++it) 
    (*it)->ReadSpec(ispec_file_);
}

void Species::ScalePositions() {
  for (auto it=members_.begin(); it!=members_.end(); ++it) 
    (*it)->ScalePosition();
}

void Species::CleanUp() {
  for (auto it=members_.begin(); it!=members_.end(); ++it)
    delete (*it);
}

void Species::ArrangeMembers() {
  if (GetInsertionType().compare("interactor_crystal") == 0)
    CrystalArrangement();
  else if (GetInsertionType().compare("centered_oriented") == 0 )
    CenteredOrientedArrangement();
  else
    warning("Arrangement not recognized and ArrangeMembers not overwritten by species!\n");
}

void Species::CenteredOrientedArrangement() {
  // This is redundant for filaments, since they have already inserted themselves properly.
  double pos[3] = {0,0,0};
  double u[3] = {0,0,0};
  u[params_->n_dim-1] = 1.0;
  for (auto it=members_.begin(); it!=members_.end(); ++it) {
    (*it)->SetPosition(pos);
    (*it)->SetOrientation(u);
  }
}

void Species::CrystalArrangement() {
  double d = members_[0]->GetDiameter();
  double l = members_[0]->GetLength();
  int n_dim = params_->n_dim;
  double R = space_->radius;
  double pos[3] = {0,0,0};
  double u[3] = {0,0,0};
  u[n_dim-1] = 1;
  double nX_max = floor((2.0*R/d)-1);
  double nY_max = floor((2.0*R-d)/(l+d));
  double nZ_max = (n_dim == 3 ? nX_max : 1);
  double N = nX_max*nY_max*nZ_max;
  if (n_members_ > N) {
    error_exit("Number of members in species exceeds maximum possible for crystal in system radius! Max possible: %d", (int) N);
  }
  double fraction = N/n_members_;
  if (fraction < 1) fraction = 1;
  if (l == 0) {
    if (n_dim == 2)
      nY_max = floor(pow(n_members_, 1.0/n_dim));
    else if (n_dim == 3)
      nY_max = ceil(pow(n_members_, 1.0/n_dim));
  }
  if (n_dim == 2) {
    nX_max = ceil(n_members_/nY_max);
  }
  else if (n_dim == 3) {
    nX_max = floor(sqrt(n_members_/nY_max));
    if (nX_max == 0) nX_max = 1;
    nZ_max = ceil(n_members_/(nX_max*nY_max));
  }
  for (int i=0; i<n_dim; ++i) {
    pos[i] = -R+0.5*d;
  }
  pos[n_dim-1] += 0.5*l;
  int inserted = 0;
  int insert_x = 0;
  int insert_y = 0;
  int insert_z = 0;
  bool shift = false;
  double diff_y = (2*R - nY_max*(l+d))/nY_max;
  double diff_x = (2*R - nX_max*d)/nX_max;
  double diff_z = (2*R - nZ_max*d)/nZ_max;
  int u0 = 1;
  for (auto it=members_.begin(); it!=members_.end(); ++it) {
    // Check crystal orientation type
    if (params_->uniform_crystal == 0) {
      // Random orientations
      u[n_dim-1] = (gsl_rng_uniform_int(rng_.r,2) == 0 ? 1 : -1);
    }
    else if (params_->uniform_crystal == 2) {
      // Stagger orientations
      u[n_dim-1] = -u[n_dim-1];
    }
    // Insert object
    (*it)->InsertAt(pos,u);
    inserted++;
    // Update next position
    pos[n_dim-1] += l+d+diff_y;
    if (++insert_y == nY_max) {
      if (params_->uniform_crystal == 2 && shift) {
        // Stagger orientations row-wise
        u[n_dim-1] = u0;
        u0 = -u0;
      }
      // Stagger positions column-wise
      shift = !shift;
      insert_y = 0;
      pos[n_dim-1] = -R+0.5*(l+d) + (shift ? 0.5*(l+d+diff_y) : 0);
      pos[0] += d + diff_x;
      if (++insert_x == nX_max) {
        if (inserted < n_members_ && n_dim == 2) {
          error_exit("Ran out of room while arranging crystal in 2D! Arranged %d/%d",inserted,n_members_);
        }
        else if (n_dim == 2) continue;
        insert_x = 0;
        pos[0] = -R+0.5*d;
        pos[1] += d + diff_z;
        if (++insert_z == nZ_max && inserted < n_members_) {
          error_exit("Ran out of room while arranging crystal in 3D! Arranged %d/%d",inserted,n_members_);
        }
      }
    }
  }
}

