#include "simcore/struct_analysis.hpp"

// Assumes 2d for now. In the future, could do 2d projections of 3d data
void StructAnalysis::Init(system_parameters *params) {
  count_ = 0;
  n_objs_ = 0;
  n_overlaps_ = 0;
  n_crossings_init_ = 0;
  n_crossings_complete_ = 0;
  params_ = params;
  n_dim_ = params_->n_dim;
  local_order_analysis_ = params_->local_order_analysis;
  polar_order_analysis_ = params_->polar_order_analysis;
  overlap_analysis_ = params_->overlap_analysis;
  density_analysis_ = params_->density_analysis;
  if (n_dim_ != 2) {
    Logger::Error("3D structure analysis is not yet implemented");
  }
  local_bin_width_ = params_->local_order_bin_width;
  density_bin_width_ = params_->density_bin_width;
  if (local_order_analysis_) {
    /* Array width is the width of the pdf (etc) histograms in units of sigma
       and the bin width is also in units of sigma, the smaller, the higher the
       total resolution */
    local_bins_1d_ = (int)params_->local_order_width / local_bin_width_;
    n_bins_1d_ = local_bins_1d_;
    // Guarantee n_bins_1d_ is even, for sanity sake
    if (n_bins_1d_ % 2 != 0) {
      n_bins_1d_++;
    }
    n_bins_ = n_bins_1d_ * n_bins_1d_;
    /* Issue a RAM warning if the arrays are very large */
    if (6.0 * n_bins_ * 4.0 / 1e9 > 1) {
      Logger::Warning("Local structure content to exceed %2.2f GB of RAM!",
                      6.0 * n_bins_ * 4.0 / 1e9);
    }
    /* initialize arrays */
    pdf_array_ = new float[n_bins_];
    nematic_array_ = new float[n_bins_];
    polar_array_ = new float[n_bins_];
    pdf_array_temp_ = new float[n_bins_];
    nematic_array_temp_ = new float[n_bins_];
    polar_array_temp_ = new float[n_bins_];
    std::fill(pdf_array_, pdf_array_ + n_bins_, 0.0);
    std::fill(nematic_array_, nematic_array_ + n_bins_, 0.0);
    std::fill(polar_array_, polar_array_ + n_bins_, 0.0);
    std::fill(pdf_array_temp_, pdf_array_temp_ + n_bins_, 0.0);
    std::fill(nematic_array_temp_, nematic_array_temp_ + n_bins_, 0.0);
    std::fill(polar_array_temp_, polar_array_temp_ + n_bins_, 0.0);
  }
  if (density_analysis_) {
    /* Array width is the width of the pdf (etc) histograms in units of sigma
       and the bin width is also in units of sigma, the smaller, the higher the
       total resolution */
    density_bins_1d_ = (int)2.0 * params_->system_radius / density_bin_width_;
    n_bins_1d_ = density_bins_1d_;
    // Guarantee n_bins_1d_ is even, for sanity sake
    // if (n_bins_1d_%2 != 0) {
    // n_bins_1d_++;
    //}
    n_bins_ = n_bins_1d_ * n_bins_1d_;
    /* Issue a RAM warning if the arrays are very large */
    if (6.0 * n_bins_ * 4.0 / 1e9 > 1) {
      Logger::Warning("Local structure content to exceed %2.2f GB of RAM!",
                      6.0 * n_bins_ * 4.0 / 1e9);
    }
    /* initialize arrays */
    density_array_ =
        (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n_bins_);
    fft_array_ = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n_bins_);
    ZeroDensityArray();
    fft_plan =
        fftw_plan_dft_2d(density_bins_1d_, density_bins_1d_, density_array_,
                         fft_array_, FFTW_FORWARD, FFTW_MEASURE);
    structure_array_ = new float[n_bins_];
    std::fill(structure_array_, structure_array_ + n_bins_, 0.0);
  }
  if (overlap_analysis_) {
    std::string overlap_file_name = params_->run_name + ".overlaps";
    overlap_file_.open(overlap_file_name, std::ios::out);
    overlap_file_ << "time n_instant_bond_overlaps n_total_crossings_init "
                     "n_total_crossings_complete\n";
  }
}

void StructAnalysis::ZeroDensityArray() {
  for (int i = 0; i < density_bins_1d_; ++i) {
    for (int j = 0; j < density_bins_1d_; ++j) {
      density_array_[i * density_bins_1d_ + j][0] = 0; // real value
      density_array_[i * density_bins_1d_ + j][1] = 0; // complex value
    }
  }
}

void StructAnalysis::SetNumObjs(int nobj) {
  // if (n_objs_ > 0) {
  // Logger::Error("Object number in struct analysis reset after initialization.
  // This is unexpected.");
  //}
  n_objs_ = nobj;
  // polar_order_ = new float[n_objs_];
  // contact_number_ = new float[n_objs_];
  // std::fill(polar_order_,polar_order_+n_objs_,0.0);
  // std::fill(contact_number_,contact_number_+n_objs_,0.0);
}

void StructAnalysis::Clear() {
  if (local_order_analysis_) {
    WriteStructData();
    delete[] pdf_array_;
    delete[] nematic_array_;
    delete[] polar_array_;
    delete[] pdf_array_temp_;
    delete[] nematic_array_temp_;
    delete[] polar_array_temp_;
  }
  if (overlap_file_.is_open()) {
    overlap_file_.close();
  }
  if (density_analysis_) {
    WriteDensityData();
    delete[] structure_array_;
    fftw_destroy_plan(fft_plan);
    fftw_free(density_array_);
    fftw_free(fft_array_);
  }

  // if (n_objs_ > 0) {
  // delete[] polar_order_;
  // delete[] contact_number_;
  //}
}

void StructAnalysis::WriteStructData() {
  if (count_ < 1 || n_objs_ < 1) {
    Logger::Warning(
        "Time average count and/or object number not initialized in "
        "StructAnalysis. Skipping local order output.");
    return;
  }
  std::string pdf_file_name = params_->run_name + ".local_pdf";
  std::string nematic_file_name = params_->run_name + ".local_nematic";
  std::string polar_file_name = params_->run_name + ".local_polar";
  pdf_file_.open(pdf_file_name, std::ios::out);
  if (!pdf_file_.is_open()) {
    Logger::Error("Output file %s did not open", pdf_file_name.c_str());
  }
  float vol = 4 * params_->system_radius * params_->system_radius;
  float rhoNinv = vol / (n_objs_ * n_objs_);
  n_bins_1d_ = local_bins_1d_;
  for (int i = 0; i < n_bins_1d_; ++i) {
    for (int j = 0; j < n_bins_1d_; ++j) {
      pdf_file_ << " " << rhoNinv * pdf_array_[i * n_bins_1d_ + j] / count_;
    }
    pdf_file_ << "\n";
  }
  pdf_file_.close();

  nematic_file_.open(nematic_file_name, std::ios::out);
  if (!nematic_file_.is_open()) {
    Logger::Error("Output file %s did not open", nematic_file_name.c_str());
  }
  for (int i = 0; i < n_bins_1d_; ++i) {
    for (int j = 0; j < n_bins_1d_; ++j) {
      // float weight = pdf_array_[i*n_bins_1d_+j];
      // weight = (weight > 0 ? weight : 1);
      // weight = pdf_array_[i*n_bins_1d_+j]/count_;
      nematic_file_ << " " << nematic_array_[i * n_bins_1d_ + j] / count_;
    }
    nematic_file_ << "\n";
  }
  nematic_file_.close();

  polar_file_.open(polar_file_name, std::ios::out);
  if (!polar_file_.is_open()) {
    Logger::Error("Output file %s did not open", polar_file_name.c_str());
  }
  for (int i = 0; i < n_bins_1d_; ++i) {
    for (int j = 0; j < n_bins_1d_; ++j) {
      // float weight = pdf_array_[i*n_bins_1d_+j];
      // weight = (weight > 0 ? weight : 1);
      polar_file_ << " " << polar_array_[i * n_bins_1d_ + j] / count_;
    }
    polar_file_ << "\n";
  }
  polar_file_.close();
}

void StructAnalysis::WriteDensityData() {
  if (count_ < 1 || n_objs_ < 1) {
    Logger::Warning("Time average count and/or object number not initialized\
        in StructAnalysis. Skipping density output.");
    return;
  }
  std::string structure_file_name = params_->run_name + ".structure";
  structure_file_.open(structure_file_name, std::ios::out);
  n_bins_1d_ = density_bins_1d_;
  for (int i = 0; i < n_bins_1d_; ++i) {
    for (int j = 0; j < n_bins_1d_; ++j) {
      structure_file_ << " "
                      << structure_array_[i * n_bins_1d_ + j] / SQR(n_objs_) /
                             count_;
    }
    structure_file_ << "\n";
  }
  structure_file_.close();
}

void StructAnalysis::CalculateStructurePair(
    std::vector<Interaction>::iterator pix) {
  if (local_order_analysis_) {
    CalculateLocalOrderPair(pix);
  }
  if (polar_order_analysis_) {
    CalculatePolarOrderPair(pix);
  }
  if (overlap_analysis_) {
    CountOverlap(pix);
  }
}

void StructAnalysis::CalculatePolarOrderPair(
    std::vector<Interaction>::iterator pix) {
  Object *obj1 = pix->obj1;
  Object *obj2 = pix->obj2;
  /* For now, ignore intra-filament correlations */
  if (obj1->GetMeshID() == obj2->GetMeshID()) {
    return;
  }
  double dr2 = pix->dr_mag2;
  double expdr2 = exp(-dr2);
  double const *const u1 = obj1->GetInteractorOrientation();
  double const *const u2 = obj2->GetInteractorOrientation();
  double u1_dot_u2 = dot_product(n_dim_, u1, u2);
  // if (u1_dot_u2 > 1 || u1_dot_u2 < -1) {
  // std::cout << "error 1: " << u1_dot_u2 << "\n";
  //}
  pix->polar_order = u1_dot_u2 * expdr2;
  pix->contact_number = expdr2;
  if (pix->polar_order != pix->polar_order ||
      pix->contact_number != pix->contact_number) {
    fprintf(stderr, "%2.2f %2.2f %2.2f %2.2f %2.2f %2.2f %2.2f\n",
        u1[0], u1[1], u2[0], u2[1], dr2, u1_dot_u2, expdr2);
  }
  // obj1->AddPolarOrder(u1_dot_u2*expdr2);
  // obj2->AddPolarOrder(u1_dot_u2*expdr2);
  // obj1->AddContactNumber(expdr2);
  // obj2->AddContactNumber(expdr2);
}

void StructAnalysis::BinDensity(Object *obj) {
  double const *const r1 = obj->GetInteractorScaledPosition();
  double const *const u1 = obj->GetInteractorOrientation();
  double const l1 = obj->GetInteractorLength();

  /* Check to see if we are only interested in COMs of objects */
  if (params_->density_com_only) {
    int index[2];
    for (int i = 0; i < 2; ++i) {
      index[i] =
          (int)(floor(r1[i] / density_bin_width_) + 0.5 * density_bins_1d_);
    }
    BinArray(index[0], index[1], 0, false);
  }

  /* Otherwise, we want to bin the entire length of the object, assuming linear
     shape */
  else {
    double site1[2] = {0, 0};
    double site2[2] = {0, 0};

    for (int i = 0; i < n_dim_; ++i) {
      site1[i] = 2.0 * params_->system_radius * r1[i] - 0.5 * l1 * u1[i];
      site2[i] = 2.0 * params_->system_radius * r1[i] + 0.5 * l1 * u1[i];
    }
    PopulateLineDensity(site1, site2);
  }
}

void StructAnalysis::CalculateLocalOrderPair(
    std::vector<Interaction>::iterator pix) {
  Object *obj1 = pix->obj1;
  Object *obj2 = pix->obj2;

  /* For now, ignore intra-filament correlations */
  if (obj1->GetMeshID() == obj2->GetMeshID()) {
    return;
  }

  double const *const r1 = obj1->GetInteractorPosition();
  double const *const r2 = obj2->GetInteractorPosition();
  double const *const u1 = obj1->GetInteractorOrientation();
  double const *const u2 = obj2->GetInteractorOrientation();
  double const l1 = obj1->GetInteractorLength();
  double const l2 = obj2->GetInteractorLength();

  /* Rotation angles that are used to rotate into the
     reference frame of the corresponding objects when
     aligned along the y-axis */
  double rotate_angle_1 = SIGNOF(u1[0]) * acos(u1[1]);
  double rotate_angle_2 = SIGNOF(u2[0]) * acos(u2[1]);

  /* Calculate dot product for nematic, polar order */

  double u1_dot_u2 = dot_product(n_dim_, u1, u2);

  // The vector pointing from object 1 center to object 2 center
  double r12[2] = {0, 0};
  // The end sites of object 1 relative to object 2
  double site11[2] = {0, 0};
  double site12[2] = {0, 0};
  // The end sites of object 2 relative to object 1
  double site21[2] = {0, 0};
  double site22[2] = {0, 0};

  for (int i = 0; i < n_dim_; ++i) {
    r12[i] = r2[i] - r1[i];
    site21[i] = r12[i] - 0.5 * l2 * u2[i];
    site22[i] = r12[i] + 0.5 * l2 * u2[i];
    site11[i] = -r12[i] - 0.5 * l1 * u1[i];
    site12[i] = -r12[i] + 0.5 * l1 * u1[i];
  }

  /* Need to rotate object 2 sites about origin by
     rotation angle 1 to put them into the reference
     frame of object 1 when aligned along the y axis */

  double temp[2];
  double cos1 = cos(rotate_angle_1);
  double sin1 = sin(rotate_angle_1);
  temp[0] = cos1 * site21[0] - sin1 * site21[1];
  temp[1] = sin1 * site21[0] + cos1 * site21[1];
  for (int i = 0; i < n_dim_; ++i) {
    site21[i] = temp[i];
  }
  temp[0] = cos1 * site22[0] - sin1 * site22[1];
  temp[1] = sin1 * site22[0] + cos1 * site22[1];
  for (int i = 0; i < n_dim_; ++i) {
    site22[i] = temp[i];
  }

  /* Now, I need to populate the bins of pdf etc histograms
     that contain the line joining site21 and site22 in the
     reference frame of object 1 when aligned along the
     y-axis */
  PopulateLineLocal(site21, site22, u1_dot_u2);

  /* Need to rotate object 1 sites about origin by
     rotation angle 2 */
  double cos2 = cos(rotate_angle_2);
  double sin2 = sin(rotate_angle_2);
  temp[0] = cos2 * site11[0] - sin2 * site11[1];
  temp[1] = sin2 * site11[0] + cos2 * site11[1];
  for (int i = 0; i < n_dim_; ++i) {
    site11[i] = temp[i];
  }
  temp[0] = cos2 * site12[0] - sin2 * site12[1];
  temp[1] = sin2 * site12[0] + cos2 * site12[1];
  for (int i = 0; i < n_dim_; ++i) {
    site12[i] = temp[i];
  }

  /* Populate the bins that contain line joining site 11
     and site 12 in reference frame of object 2 when
     aligned along y-axis */
  PopulateLineLocal(site11, site12, u1_dot_u2);
}

/* PopulateLine uses Bresenham line drawing algorithm to
   populate the bins that define a line joining site1 and
   site2 in a discretized (pixelized) space */
void StructAnalysis::PopulateLineLocal(double *site1, double *site2,
                                       double dotprod) {
  // Convert positions into bin numbers
  int index1[2], index2[2];
  for (int i = 0; i < 2; ++i) {
    index1[i] =
        (int)(floor(site1[i] / local_bin_width_) + 0.5 * local_bins_1d_);
    index2[i] =
        (int)(floor(site2[i] / local_bin_width_) + 0.5 * local_bins_1d_);
  }
  BinLine(index1[0], index1[1], index2[0], index2[1], dotprod, true);
}

/* PopulateLine uses Bresenham line drawing algorithm to
   populate the bins that define a line joining site1 and
   site2 in a discretized (pixelized) space */
void StructAnalysis::PopulateLineDensity(double *site1, double *site2) {
  // Convert positions into bin numbers
  int index1[2], index2[2];
  for (int i = 0; i < 2; ++i) {
    index1[i] =
        (int)(floor(site1[i] / density_bin_width_) + 0.5 * density_bins_1d_);
    index2[i] =
        (int)(floor(site2[i] / density_bin_width_) + 0.5 * density_bins_1d_);
  }
  BinLine(index1[0], index1[1], index2[0], index2[1], 0.0, false);
}

/* Draw line using Bresenham line drawing algorithm
   using integer algebra. Only draws line upward and
   to the right, so we need to reorient the line
   accordingly */
void StructAnalysis::BinLine(int x0, int y0, int x1, int y1, double dotprod,
                             bool is_local) {
  /* Checks if line is less than 45 degrees */
  if (ABS(y1 - y0) < ABS(x1 - x0)) {
    /* Ensures that slope is only increasing */
    if (x0 > x1) {
      BinLineLow(x1, y1, x0, y0, dotprod, is_local);
    } else {
      BinLineLow(x0, y0, x1, y1, dotprod, is_local);
    }
  } else {
    /* If line is more than 45 degrees, switch x and y axis */
    if (y0 > y1) {
      /* Ensures that slope is only increasing */
      BinLineHigh(x1, y1, x0, y0, dotprod, is_local);
    } else {
      BinLineHigh(x0, y0, x1, y1, dotprod, is_local);
    }
  }
}

/* The simple Bresenham line algorithm using integer algebra.
   Assumes that the line is less than 45 degrees and the site
   positions are integers */
void StructAnalysis::BinLineLow(int x0, int y0, int x1, int y1, double dotprod,
                                bool is_local) {
  int dx = x1 - x0;
  int dy = y1 - y0;
  int yi = 1;
  if (dy < 0) {
    yi = -1;
    dy = -dy;
  }
  int D = 2 * dy - dx;
  int y = y0;
  for (int x = x0; x < x1; ++x) {
    BinArray(x, y, dotprod, is_local);
    if (D > 0) {
      y = y + yi;
      D = D - 2 * dx;
    }
    D = D + 2 * dy;
  }
}

/* The Bresenham line drawing algorithm with inverted x and y
   axes, for the case that the line slope is greater than 45
   degrees */
void StructAnalysis::BinLineHigh(int x0, int y0, int x1, int y1, double dotprod,
                                 bool is_local) {
  int dx = x1 - x0;
  int dy = y1 - y0;
  int xi = 1;
  if (dx < 0) {
    xi = -1;
    dx = -dx;
  }
  int D = 2 * dx - dy;
  int x = x0;
  for (int y = y0; y < y1; ++y) {
    BinArray(x, y, dotprod, is_local);
    if (D > 0) {
      x = x + xi;
      D = D - 2 * dy;
    }
    D = D + 2 * dx;
  }
}

/* Increments x,y bin of the pdf histograms. Converts x,y to a
   1d coordinate and incorporates the dotproduct alignment of the
   two objects for nematic and polar order */
void StructAnalysis::BinArray(int x, int y, double dotprod, bool is_local) {
  /* If the coordinate is out of range of our region of interest
     for whatever reason, ignore the coordinate */
  n_bins_1d_ = (is_local ? local_bins_1d_ : density_bins_1d_);
  if (x < 0)
    x += n_bins_1d_;
  if (y < 0)
    y += n_bins_1d_;
  if (x > n_bins_1d_ - 1)
    x -= n_bins_1d_;
  if (y > n_bins_1d_ - 1)
    y -= n_bins_1d_;
  if (x < 0 || y < 0 || x > n_bins_1d_ - 1 || y > n_bins_1d_ - 1) {
    Logger::Warning("Invalid value found in BinArray in struct_analysis.cpp");
    return;
  }
  int index = n_bins_1d_ * y + x;
  std::lock_guard<std::mutex> lk(mtx_);
  if (is_local) {
    pdf_array_temp_[index] += 1.0;
    nematic_array_temp_[index] += (2 * dotprod * dotprod - 1);
    polar_array_temp_[index] += dotprod;
  } else {
    density_array_[index][0] += 1.0;
  }
}

void StructAnalysis::AverageStructure() {
  if (local_order_analysis_) {
    for (int i = 0; i < n_bins_; ++i) {
      float weight = pdf_array_temp_[i];
      weight = (weight > 0 ? weight : 1);
      pdf_array_[i] += pdf_array_temp_[i];
      nematic_array_[i] += nematic_array_temp_[i] / weight;
      polar_array_[i] += polar_array_temp_[i] / weight;
      pdf_array_temp_[i] = 0;
      nematic_array_temp_[i] = 0;
      polar_array_temp_[i] = 0;
    }
  }
  if (overlap_analysis_) {
    crossing_list_.CompareLists(&n_crossings_init_, &n_crossings_complete_);
    overlap_file_ << ((params_->i_step) * params_->delta) << " " << n_overlaps_
                  << " " << 0.5 * n_crossings_init_ << " "
                  << 0.5 * n_crossings_complete_ << "\n";
    n_overlaps_ = 0;
  }
  if (density_analysis_) {
    // Turn the histogram into a density by averaging by the number of objects
    // in the system
    for (int i = 0; i < density_bins_1d_; ++i) {
      for (int j = 0; j < density_bins_1d_; ++j) {
        density_array_[i * density_bins_1d_ + j][0] /= n_objs_;
      }
    }
    // Get the fft of the instantaneous density
    fftw_execute(fft_plan);
    for (int i = 0; i < density_bins_1d_; ++i) {
      for (int j = 0; j < density_bins_1d_; ++j) {
        structure_array_[i * density_bins_1d_ + j] +=
            (SQR(fft_array_[i * density_bins_1d_ + j][0]) +
             SQR(fft_array_[i * density_bins_1d_ + j][1]));
      }
    }
    ZeroDensityArray();
  }
}

/* In order to count the number of overlaps, I examine the vector pointing
 * between the heads of the two objects and then the vector pointing between
 * the tails of the two objects. If these objects point in opposite (relative)
 * directions (e.g. if their dot product is negative) then the filaments must
 * be crossing each other */
void StructAnalysis::CountOverlap(std::vector<Interaction>::iterator pix) {
  Object *obj1 = pix->obj1;
  Object *obj2 = pix->obj2;
  if (obj1->GetMeshID() == obj2->GetMeshID()) {
    return;
  }
  if (ABS(pix->dr_mag2) < 1e-8) {
    obj1->HasOverlap(true);
    obj2->HasOverlap(true);
    AddOverlap();
    // CountOverlapEvents(obj1->GetMeshID(), obj2->GetMeshID());
    crossing_list_.AddNeighbors(obj1->GetMeshID(), obj2->GetMeshID());
  }
}

// void StructAnalysis::CountOverlapEvents(int mid1, int mid2, bool
// is_overlapping) {
//// If pair was previously overlapping
// if (crossing_list_.AreNeighbors(mid1, mid2)) {
//// Still overlapping. No new information
// if (is_overlapping) {
// return;
//}
//// No longer overlapping. Crossing event has ended.
// else {
// crossing_list_.RemoveNeighbors(mid1,mid2);
//}
//}
//// Pair not previously overlapping
// else {
//// Fresh crossing event. Count crossing event and track pair
// if (is_overlapping) {
// crossing_list_.AddNeighbors(mid1,mid2);
//}
//// Still not overlapping.
// else {
// return;
//}
//}
//}

void StructAnalysis::AddOverlap() {
  std::lock_guard<std::mutex> lk(mtx_);
  n_overlaps_++;
}
void StructAnalysis::AddCrossingComplete() {
  // std::lock_guard<std::mutex> lk(mtx_);
  n_crossings_complete_++;
}
void StructAnalysis::AddCrossingInit() {
  // std::lock_guard<std::mutex> lk(mtx_);
  n_crossings_init_++;
}
