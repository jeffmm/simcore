#include "object.h"

Object::Object() {
  // Initialize object ID
  oid_ = ++next_oid_;
  // Initialize RNG and set new object seed
  rng_.Init(_seed_);
  _seed_ = 1000*gsl_rng_get(rng_.r)-gsl_rng_get(rng_.r);
  // Set some defaults
  std::fill(position_,position_+3,0.0);
  std::fill(prev_position_,prev_position_+3,0.0);
  std::fill(prev_orientation_,prev_orientation_+3,0.0);
  std::fill(scaled_position_,scaled_position_+3,0.0);
  std::fill(orientation_,orientation_+3,0.0);
  std::fill(force_,force_+3,0.0);
  std::fill(torque_, torque_+3, 0.0);
  std::fill(dr_zero_, dr_zero_+3, 0.0);
  draw_ = draw_type::orientation;
  color_ = 0;
  diameter_ = 1;
  length_ = 0;
  p_energy_ = 0;
  interacting_ = false;
  is_mesh_ = false;
  dr_tot_ = 0;
  mesh_id_ = 0;
  polar_order_ = 0;
  contact_number_ = 0;
  n_contact_ = 0;
  has_overlap_ = false;
  in_flock_ = 0;
  flock_change_state_ = 0;
}

int Object::next_oid_ = 0;
long Object::_seed_ = 7777777;
system_parameters * Object::params_ = nullptr;
space_struct * Object::space_ = nullptr;
int Object::n_dim_ = 0;
double Object::delta_ = 0;

void Object::SetSeed(long seed) {
  _seed_ = seed;
}
void Object::SetParams(system_parameters * params) {
  params_ = params;
}
void Object::SetSpace(space_struct * space) {
  space_ = space;
}
void Object::SetNDim(int n_dim) {
  n_dim_ = n_dim;
}
void Object::SetDelta(double delta) {
  delta_ = delta;
}
// Trivial Get/Set functions
int const Object::GetOID() const {
  return oid_;
}
void Object::SetPosition(double const *const pos) {
  std::copy(pos, pos+n_dim_, position_);
}
void Object::SetScaledPosition(double const *const spos) {
  std::copy(spos, spos+n_dim_, scaled_position_);
}
void Object::SetOrientation(double const * const u) {
  std::copy(u, u+n_dim_, orientation_);
}
void Object::SetPrevPosition(double const * const ppos) {
  std::copy(ppos, ppos+n_dim_, prev_position_);
}
void Object::SetPrevOrientation(double const * const pu) {
  std::copy(pu, pu+n_dim_, prev_orientation_);
}
void Object::SetDiameter(double new_diameter) {
  diameter_ = new_diameter;
}
void Object::SetLength(double new_length) {
  length_ = new_length;
}
void Object::AddForce(double const * const f) {
  for (int i=0; i<n_dim_; ++i)
    force_[i]+=f[i];
}
void Object::SubForce(double const * const f) {
  for (int i=0; i<n_dim_; ++i)
    force_[i]-=f[i];
}
void Object::SetForce(double const * const f) {
  std::copy(f, f+n_dim_, force_);
}
void Object::AddTorque(double const * const t) {
  for (int i=0; i<3; ++i)
    torque_[i]+=t[i];
}
void Object::SubTorque(double const * const t) {
  for (int i=0; i<3; ++i)
    torque_[i]-=t[i];
}
void Object::SetTorque(double const * const t) {
  std::copy(t, t+3, torque_);
}
void Object::AddPotential(double const p) {
  p_energy_ += p;
}
void Object::AddPolarOrder(double const po) {
  //std::lock_guard<std::mutex> lk(mtx_);
  polar_order_ += po;
}
void Object::AddContactNumber(double const cn) {
  //std::lock_guard<std::mutex> lk(mtx_);
  contact_number_ += cn;
}
void Object::CalcPolarOrder() {
  if (contact_number_ > 0) {
    polar_order_ /= contact_number_;
  }
  else {
    polar_order_ = 0;
  }
  //if (polar_order_ < -1 || polar_order_ > 1) {
    //std::cout << "error 2: " << polar_order_ << " " << contact_number_ << "\n";
  //}
}
void Object::ZeroPolarOrder() {
  contact_number_ = 0;
  polar_order_ = 0;
}
void Object::SetInteractor(bool ix) {
  interacting_ = ix;
}
double const * const Object::GetPosition() {
  return position_;
}
double const * const Object::GetPrevPosition() {
  return prev_position_;
}
double const * const Object::GetScaledPosition() {
  return scaled_position_;
}
double const * const Object::GetOrientation() {
  return orientation_;
}
double const * const Object::GetForce() {
  return force_;
}
double const * const Object::GetTorque() {
  return torque_;
}
double const Object::GetDiameter() {
  return diameter_;
}
double const Object::GetLength() {
  return length_;
}
double const Object::GetPotentialEnergy() {
  return p_energy_;
}
double const Object::GetPolarOrder() {
  return polar_order_;
}
double const Object::GetContactNumber() {
  return contact_number_;
}
bool const Object::IsInteractor() {
  return interacting_;
}
bool const Object::IsMesh() {
  return is_mesh_;
}
void Object::HasOverlap(bool overlap) {
  has_overlap_=overlap;
}
void Object::SetFlockType(int in_flock) {
  in_flock_ = in_flock;
}
int Object::GetFlockType() {
  return in_flock_;
}
void Object::SetFlockChangeState(int fcs) {
  flock_change_state_ = fcs;
}
int Object::GetFlockChangeState() {
  return flock_change_state_;
}
int const Object::GetMeshID() const {
  return mesh_id_;
}
void Object::SetMeshID(int mid) {
  mesh_id_ = mid;
}
void Object::ToggleIsMesh() {
  is_mesh_ = !is_mesh_;
}
species_id const Object::GetSID() {
  return sid_;
}
void Object::SetSID(species_id sid) {
  sid_ = sid;
}

// Virtual functions
// Returns true if object is out of bounds
bool Object::CheckBounds(double buffer) {
  if (space_->type == +boundary_type::none) return false;
  double const * const r = GetInteractorPosition();
  double const * const u = GetInteractorOrientation();
  double const l = GetInteractorLength();
  double const d = GetInteractorDiameter();
  double r_mag = 0.0;
  double z0 = 0.0;
  double r_boundary = space_->radius;
  buffer += 0.5*diameter_;
  if (buffer > r_boundary) {
    error_exit("Cannot insert object into system!\n");
  }
  int sign = ( l>0 ? SIGNOF(dot_product(n_dim_,r,u)) : 0 );
  if (space_->type == +boundary_type::box)  {
    if (space_->n_periodic == n_dim_) return false;
    for (int j=space_->n_periodic; j<n_dim_; ++j) {
      double r_far = r[j] + sign*0.5*l*u[j];
      if (ABS(r_far) > (r_boundary - buffer)) return true;
    }
    return false;
  }
  if (space_->type == +boundary_type::budding && 
      r[n_dim_-1] > space_->bud_neck_height) {
    z0 = space_->bud_height;
    r_boundary = space_->bud_radius;
  }
  for (int i=0;i<n_dim_-1;++i) {
    r_mag += SQR(r[i] + sign*0.5*l*u[i]);
  }
  r_mag += SQR(r[n_dim_-1] + sign*0.5*l*u[n_dim_-1] - z0);
  return (r_mag > SQR(r_boundary - buffer));
}

//bool Object::CheckBounds() {
  //double const * const pos = GetPosition();
  //double buffer = diameter_;
  //if (space_->n_periodic == n_dim_)
    //buffer = 0;
  //double R = space_->radius;
  //bool out_of_bounds = false;
  //switch (space_->type._to_integral()) {
    //case 0: // none
      //break;
    //case 2: // sphere
      //{
        //double rmag = 0;
        //for (int i=0; i<n_dim_; ++i) {
          //rmag += pos[i]*pos[i];
        //}
        //// Check that r^2 <= R^2
        //rmag=sqrt(rmag);
        //if (rmag+buffer>R) {
          //out_of_bounds = true;
        //}
        //break;
      //}
    //case 1: // box
      //// Make sure each dimension of position, x, y etc is within the box radius
      //for (int i=0; i<n_dim_; ++i) {
        //if (ABS(pos[i]) + buffer > R) {
          //out_of_bounds = true;
        //}
      //}
      //break;
    //case 3: // budding
      //// FIXME Add checks to make sure object is placed within budding boundary...
      //// right now, just trust user knows what they are doing.
      //break;
    //default:
      //error_exit("ERROR: Boundary type not recognized.\n");
  //}
  //return out_of_bounds;
//}

void Object::InsertRandom() {
  double mag;
  double buffer = diameter_;
  if (space_->n_periodic == n_dim_)
    buffer = 0;
  double R = space_->radius;
  // XXX Check to make sure object fits inside unit cell if non-periodic
  if (R - buffer < 0) {
    error_exit("Object #%d is too large to place in system.\n system radius: "
        "%2.2f, buffer: %2.2f", GetOID(), R, buffer);
  }
  switch (space_->type) {
    // If no boundary, insert wherever
    case +boundary_type::none: // none
      for (int i=0; i<n_dim_; ++i) {
        position_[i] = (2.0*gsl_rng_uniform_pos(rng_.r)-1.0) * (R - buffer);
      }
      break;
    // box type boundary
    case +boundary_type::box: // box
      for (int i=0; i<n_dim_; ++i) {
        position_[i] = (2.0*gsl_rng_uniform_pos(rng_.r)-1.0) * (R - buffer);
      }
      break;
    // spherical boundary
    case +boundary_type::sphere: // sphere
      generate_random_unit_vector(n_dim_, position_, rng_.r);
      mag = gsl_rng_uniform_pos(rng_.r) * (R - buffer);
      for (int i=0; i<n_dim_; ++i) {
        position_[i] *= mag;
      }
      break;
    // budding yeast boundary type
    case +boundary_type::budding: // budding
      {
        double r = space_->bud_radius;
        double roll = gsl_rng_uniform_pos(rng_.r);
        double v_ratio = 0;
        if (n_dim_ == 2) {
          v_ratio = SQR(r) / (SQR(r)+SQR(R));
        }
        else {
          v_ratio = CUBE(r) / (CUBE(r)+CUBE(R));
        }
        mag = gsl_rng_uniform_pos(rng_.r);
        generate_random_unit_vector(n_dim_, position_, rng_.r);
        if (roll < v_ratio) {
          // Place coordinate in daughter cell
          mag *= (r - buffer);
          for (int i=0; i<n_dim_; ++i) {
            position_[i] *= mag;
          }
          position_[n_dim_-1] += space_->bud_height;
        }
        else {
          mag *= (R - buffer);
          for (int i=0; i<n_dim_; ++i) {
            position_[i] *= mag;
          }
        }
        break;
      }
    default:
      error_exit("ERROR: Boundary type unrecognized!\n");
  }
  generate_random_unit_vector(n_dim_, orientation_, rng_.r);
  UpdatePeriodic();
}

void Object::InsertRandomOriented(double *u){
  InsertRandom();
  normalize_vector(u,n_dim_);
  SetOrientation(u);
}

void Object::InsertAt(double *pos, double *u) {
  // Check to make sure position is inside unit cell if non-periodic
  double R = space_->radius;
  bool out_of_bounds = false;
  switch (space_->type._to_integral()) {
    case 0: // none
      break;
    case 2: // sphere
      {
        double rsq = 0;
        for (int i=0; i<n_dim_; ++i) {
          rsq += pos[i]*pos[i];
        }
        // Check that r^2 <= R^2
        if (rsq>R*R) {
          out_of_bounds = true;
        }
        break;
      }
    case 1: // box
      // Make sure each dimension of position, x, y etc is within the box radius
      for (int i=0; i<n_dim_; ++i) {
        if (ABS(pos[i]) > R) {
          out_of_bounds = true;
        }
      }
      break;
    case 3: // budding
      // FIXME Add checks to make sure object is placed within budding boundary...
      // right now, just trust user knows what they are doing.
      break;
    default:
      error_exit("ERROR: Boundary type not recognized.\n");
  }
  if (out_of_bounds) {
    error_exit("Object %d placed outside of system unit cell! System radius: " 
        "%2.2f, pos: {%2.2f %2.2f %2.2f}\n", GetOID(), R, pos[0], pos[1],
        pos[2]);
  }
  SetPosition(pos);
  normalize_vector(u,n_dim_);
  SetOrientation(u);
  UpdatePeriodic();
}

void Object::ZeroForce() {
  std::fill(force_,force_+3,0.0);
  std::fill(torque_,torque_+3,0.0);
  p_energy_ = 0.0;
}

void Object::Draw(std::vector<graph_struct*> * graph_array) {
  std::copy(scaled_position_,scaled_position_+3, g_.r);
  for (int i=space_->n_periodic; i<n_dim_; ++i) {
    g_.r[i] = position_[i];
  }
  std::copy(orientation_, orientation_+3, g_.u);
  g_.color = color_;
  g_.diameter = diameter_;
  g_.length = length_;
  g_.draw = draw_;
  graph_array->push_back(&g_);
}

// Updates scaled position, leaving position fixed
void Object::UpdatePeriodic() {
  double s[3], r[3];
  std::copy(scaled_position_, scaled_position_+3, s);
  periodic_boundary_conditions(space_->n_dim, space_->n_periodic, 
      space_->unit_cell, space_->unit_cell_inv, position_, s);
  SetScaledPosition(s);
}
void Object::UpdatePositionMP() {
  error_exit("UpdatePositionMP() needs to be overwritten. Exiting!");
}
void Object::SetColor(double const c, draw_type dtype) {
  color_ = c;
  draw_ = dtype;
}
void Object::ScalePosition() {
  for (int i=0; i<n_dim_; ++i) {
    position_[i] = 0;
    for (int j=0;j<n_dim_; ++j) {
      position_[i]+=space_->unit_cell[n_dim_*i+j]*scaled_position_[j];
    }
  }
}
int Object::GetCount() {
  return 1;
}
std::vector<Object*> Object::GetInteractors() {
  return interactors_;
}
double const * const Object::GetInteractorPosition() {
  return GetPosition();
}
double const * const Object::GetInteractorPrevPosition() {
  return GetPrevPosition();
}
double const * const Object::GetInteractorScaledPosition() {
  return GetScaledPosition();
}
double const * const Object::GetInteractorOrientation() {
  return GetOrientation();
}
double const Object::GetInteractorDiameter() {
  return GetDiameter();
}
double const Object::GetInteractorLength() {
  return GetLength();
}
double const Object::GetVolume() {
  if (n_dim_ == 2) {
    if (length_ == 0) {
      return 0.25*M_PI*diameter_*diameter_;
    }
    else {
      return length_*diameter_+0.25*M_PI*diameter_*diameter_;
    }
  }
  if (n_dim_ == 3) {
    if (length_ == 0) {
      return 1.0/6.0 * M_PI * diameter_*diameter_*diameter_;
    }
    else {
      return 1.0/6.0 * M_PI * diameter_*diameter_*diameter_ + 0.25*M_PI*length_*diameter_*diameter_;
    }
    return -1;
  }
}
double const Object::GetDrTot() {
  return dr_tot_;
}
void Object::ZeroDrTot() {
  std::copy(position_,position_+3,dr_zero_);
  dr_tot_ = 0;
}
void Object::UpdateDrTot() {
  for (int i=0;i<n_dim_;++i) {
    double dr = position_[i] - dr_zero_[i];
    dr_tot_ += dr*dr;
  }
}
bool Object::HasNeighbor(int other_oid) {
  // Generic objects are not assumed to have neighbors
  return false;
}
void Object::GiveInteraction(Interaction *ix) {
  ixs_.push_back(ix);
}
std::vector<Interaction*> * Object::GetInteractions() {
  return &ixs_;
}
void Object::ClearInteractions() {
  ixs_.clear();
}
void Object::Cleanup() {}

// Object I/O functions
void Object::WriteCheckpoint(std::fstream &ocheck) {
  void * rng_state = gsl_rng_state(rng_.r);
  size_t rng_size = gsl_rng_size(rng_.r);
  ocheck.write(reinterpret_cast<char*>(&rng_size), sizeof(size_t));
  ocheck.write(reinterpret_cast<char*>(rng_state), rng_size);
  WriteSpec(ocheck);
}
void Object::ReadCheckpoint(std::fstream &icheck) {
  if (icheck.eof()) return;
  void * rng_state = gsl_rng_state(rng_.r);
  size_t rng_size;
  icheck.read(reinterpret_cast<char*>(&rng_size), sizeof(size_t));
  icheck.read(reinterpret_cast<char*>(rng_state), rng_size);
  ReadSpec(icheck);
}
void Object::WritePosit(std::fstream &oposit){
  for(auto& pos : position_)
    oposit.write(reinterpret_cast<char*>(&pos), sizeof(pos));
  for(auto& spos : scaled_position_)
    oposit.write(reinterpret_cast<char*>(&spos), sizeof(spos));
  for(auto& u : orientation_)
    oposit.write(reinterpret_cast<char*>(&u), sizeof(u));
  oposit.write(reinterpret_cast<char*>(&diameter_), sizeof(diameter_));
  oposit.write(reinterpret_cast<char*>(&length_), sizeof(length_));
}
void Object::ReadPosit(std::fstream &iposit){
  if (iposit.eof()) return;
  for(auto& pos : position_)
    iposit.read(reinterpret_cast<char*>(&pos), sizeof(pos));
  for(auto& spos : scaled_position_)
    iposit.read(reinterpret_cast<char*>(&spos), sizeof(spos));
  for(auto& u : orientation_)
    iposit.read(reinterpret_cast<char*>(&u), sizeof(u));
  iposit.read(reinterpret_cast<char*>(&diameter_), sizeof(diameter_));
  iposit.read(reinterpret_cast<char*>(&length_), sizeof(length_));
}
void Object::WriteSpec(std::fstream &ospec) {
  WritePosit(ospec);
}
void Object::ReadSpec(std::fstream &ispec) {
  ReadPosit(ispec);
}
void Object::ReadPositFromSpec(std::fstream &ispec) {
  ReadPosit(ispec);
}
void Object::GetAvgPosition(double * ap) {
  std::copy(position_,position_+3,ap);
}
void Object::GetAvgOrientation(double * au) {
  std::copy(orientation_,orientation_+3, au);
}
void Object::SetAvgPosition() {
  // Nothing to be done
}
void Object::SetRNGState(const std::string& filename) {
  // Load the rng state from binary file
  FILE* pfile = fopen(filename.c_str(), "r");
  auto retval = gsl_rng_fread(pfile, rng_.r);
  if (retval != 0) {
    std::cout << "Reading rng state failed " << retval << std::endl;
  }
}

void Object::Report() {
  fprintf(stderr, "    OID: %d\n",GetOID());
  fprintf(stderr, "      r: {%2.2f %2.2f %2.2f}\n", position_[0], 
      position_[1], position_[2]);
  fprintf(stderr, "      u: {%2.2f %2.2f %2.2f}\n", orientation_[0], 
      orientation_[1], orientation_[2]);
  fprintf(stderr, "      l: %2.2f\n",length_);
  fprintf(stderr, "      d: %2.2f\n",diameter_);
}

