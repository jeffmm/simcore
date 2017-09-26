#include "object.h"
#include "minimum_distance.h"

Object::Object() {
  // Initialize object ID
  oid_ = ++next_oid_;
  // Initialize RNG and set new object seed
  rng_.Init(_seed_);
  _seed_ = gsl_rng_get(rng_.r);
  // Set some defaults
  std::fill(position_,position_+3,0.0);
  std::fill(prev_position_,prev_position_+3,0.0);
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
void Object::SetDiameter(double new_diameter) {
  diameter_ = new_diameter;
}
void Object::SetLength(double new_length) {
  length_ = new_length;
}
void Object::AddForce(double const * const f) {
  for (int i=0; i<3; ++i)
    force_[i]+=f[i];
}
void Object::SubForce(double const * const f) {
  for (int i=0; i<3; ++i)
    force_[i]-=f[i];
}
void Object::SetForce(double const * const f) {
  for (int i=0; i<3; ++i)
    force_[i]=f[i];
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
  for (int i=0; i<3; ++i)
    torque_[i]=t[i];
}
void Object::AddPotential(double const p) {
  p_energy_ += p;
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
bool const Object::IsInteractor() {
  return interacting_;
}
bool const Object::IsMesh() {
  return is_mesh_;
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

// Virtual functions
bool Object::CheckBounds(double const * pos, double buffer) {
  if (space_->n_periodic == n_dim_)
    buffer = 0;
  double R = space_->radius;
  bool out_of_bounds = false;
  switch (space_->type._to_integral()) {
    case 0: // none
      break;
    case 2: // sphere
      {
        double rmag = 0;
        for (int i=0; i<n_dim_; ++i) {
          rmag += pos[i]*pos[i];
        }
        // Check that r^2 <= R^2
        rmag=sqrt(rmag);
        if (rmag+buffer>R) {
          out_of_bounds = true;
        }
        break;
      }
    case 1: // box
      // Make sure each dimension of position, x, y etc is within the box radius
      for (int i=0; i<n_dim_; ++i) {
        if (ABS(pos[i]) + buffer > R) {
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
  return out_of_bounds;
}

void Object::InsertRandom() {
  double mag;
  double buffer = diameter_;
  if (space_->n_periodic == n_dim_)
    buffer = 0;
  double R = space_->radius;
  // XXX Check to make sure object fits inside unit cell if non-periodic
  if (R - buffer < 0)
    error_exit("Object #%d is too large to place in system.\n system radius: %2.2f, buffer: %2.2f",GetOID(), R, buffer);
  switch (space_->type._to_integral()) {
    // If no boundary, insert wherever
    case 0: // none
      for (int i=0; i<n_dim_; ++i) {
        position_[i] = (2.0*gsl_rng_uniform_pos(rng_.r)-1.0) * (R - buffer);
      }
      break;
    // spherical boundary
    case 2: // sphere
      generate_random_unit_vector(n_dim_, position_, rng_.r);
      mag = gsl_rng_uniform_pos(rng_.r) * (R - buffer);
      for (int i=0; i<n_dim_; ++i) {
        position_[i] *= mag;
      }
      break;
    // box type boundary
    case 1: // box
      for (int i=0; i<n_dim_; ++i) {
        position_[i] = (2.0*gsl_rng_uniform_pos(rng_.r)-1.0) * (R - buffer);
      }
      break;
    // budding yeast boundary type
    case 3: // budding
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
    error_exit("Object %d placed outside of system unit cell! System radius: %2.2f, pos: {%2.2f %2.2f %2.2f}\n",GetOID(),R,pos[0],pos[1],pos[2]);
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
  periodic_boundary_conditions(space_->n_dim, space_->n_periodic, space_->unit_cell, space_->unit_cell_inv, position_, s);
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
    return 0.25*M_PI*diameter_*diameter_;
  }
  if (n_dim_ == 3) {
    return 1.0/6.0 * M_PI * diameter_*diameter_*diameter_;
  }
  return -1;
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
  fprintf(stderr, "      r: {%2.2f %2.2f %2.2f}\n",position_[0],position_[1],position_[2]);
  fprintf(stderr, "      u: {%2.2f %2.2f %2.2f}\n",orientation_[0],orientation_[1],orientation_[2]);
  fprintf(stderr, "      l: %2.2f\n",length_);
  fprintf(stderr, "      d: %2.2f\n",diameter_);
}

// Find the minimum distance beween two particles
void MinimumDistance(Object* o1, Object* o2, Interaction *ix, space_struct *space) {
  double const * const r1 = o1->GetInteractorPosition();
  double const * const s1 = o1->GetInteractorScaledPosition();
  double const * const u1 = o1->GetInteractorOrientation();
  double const l1 = o1->GetInteractorLength();
  double const d1 = o1->GetInteractorDiameter();
  double const * const r2 = o2->GetInteractorPosition();
  double const * const s2 = o2->GetInteractorScaledPosition();
  double const * const u2 = o2->GetInteractorOrientation();
  double const l2 = o2->GetInteractorLength();
  double const d2 = o2->GetInteractorDiameter();
  /* TODO: Think about how best to do this for general shapes, like 2d
     polygons that can represent the local surface of more complex 3d
     shapes. Perhaps assume all local surface to be triangular polygons.*/
  ix->dr_mag2 = 0;
  std::fill(ix->dr, ix->dr+3, 0.0);
  std::fill(ix->contact1, ix->contact1+3, 0.0);
  std::fill(ix->contact2, ix->contact2+3, 0.0);
  ix->buffer_mag = 0.5*(d1+d2);
  ix->buffer_mag2 = ix->buffer_mag*ix->buffer_mag;
  if (l1 == 0 && l2 == 0) 
    min_distance_point_point(space->n_dim, space->n_periodic, space->unit_cell,
                 r1, s1, r2, s2, ix->dr, &ix->dr_mag2);
  else if (l1 == 0 && l2 > 0)
    min_distance_sphere_sphero(space->n_dim, space->n_periodic, space->unit_cell,
                   r1, s1, r2, s2, u2, l2,
                   ix->dr, &ix->dr_mag2, ix->contact2);
  else if (l1 > 0 && l2 == 0) {
    min_distance_sphere_sphero(space->n_dim, space->n_periodic, space->unit_cell,
                   r2, s2, r1, s1, u1, l1,
                   ix->dr, &ix->dr_mag2, ix->contact1);
    for (int i=0;i<3;++i) {
      ix->dr[i] = -ix->dr[i];
    }
  }
  else if (l1 > 0 && l2 > 0) 
    min_distance_sphero(space->n_dim, space->n_periodic, space->unit_cell,
              r1, s1, u1, l1, r2, s2, u2, l2,
              ix->dr, &ix->dr_mag2, ix->contact1, ix->contact2);
}

