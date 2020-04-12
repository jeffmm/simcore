#include "simcore/filament_curvature_cluster_analysis.hpp"

int Cluster::n_dim_ = -1;
double *Cluster::curvature_radii_ = nullptr;
double **Cluster::curvature_centers_ = nullptr;
double **Cluster::cc_vectors_ = nullptr;
int *Cluster::handedness_ = nullptr;
double Cluster::color_ = 0;
const double Cluster::GetSqrDistanceFilament(int i) {
  return GetSqrDistance(curvature_centers_[i]);
}

const double Cluster::GetSqrDistance(double *pos) {
  double dr[3] = {0};
  double mid[3] = {0};
  double dr_mag2;
  MinimumDistance mindist;
  mindist.PointPoint(position_, position_, pos, pos, dr, &dr_mag2, mid);
  return dr_mag2;
}
void Cluster::CalcPosition() {

  std::copy(position_, position_ + 3, prev_position_);
  std::fill(position_, position_ + 3, 0);
  int size = GetSize();
  if (size == 0) {
    return;
  }
  auto first = members_.begin();
  double *pos = curvature_centers_[*first];
  for (int i = 0; i < n_dim_; ++i) {
    position_[i] = pos[i];
  }
  double tot[3] = {0};
  avg_radius_ = 0;
  avg_handedness_ = 0;
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    double rad = curvature_radii_[*it];
    int hand = handedness_[*it];
    avg_radius_ += rad;
    avg_handedness_ += hand;
    if (rad > radius_) {
      radius_ = curvature_radii_[*it];
    }
    if (it == members_.begin())
      continue;
    pos = curvature_centers_[*it];
    double ds[3] = {0};
    for (int i = 0; i < n_dim_; ++i) {
      ds[i] = position_[i] - pos[i];
      ds[i] -= NINT(ds[i]);
      tot[i] += ds[i];
    }
  }
  for (int i = 0; i < n_dim_; ++i) {
    position_[i] -= tot[i] / size;
    position_[i] -= NINT(position_[i]);
  }
  avg_radius_ /= size;
  avg_handedness_ /= size;
}

void Cluster::Merge(Cluster &other) {
  for (auto it = other.members_.begin(); it != other.members_.end(); ++it) {
    AddMember(*it);
  }
  other.members_.clear();
  CalcPosition();
}

// Check if particle is close enough to be considered still within the cluster
const bool Cluster::CheckInCluster(int i) {
  double dr[3] = {0};
  double mid[3] = {0};
  double dr_mag2;
  MinimumDistance mindist;
  double *pos = curvature_centers_[i];
  mindist.PointPoint(position_, position_, pos, pos, dr, &dr_mag2, mid);
  /* First check that the vector to the curvature center points the right way,
     and if not, enforce that dr must be less than the avg_radius */
  if (dot_product(n_dim_, cc_vectors_[i], dr) > 0 &&
      dr_mag2 > SQR(avg_radius_)) {
    return false;
  }
  // Otherwise use the maximum radius
  if (dr_mag2 > SQR(radius_)) {
    return false;
  }
  return true;
}

void Cluster::Draw(std::vector<graph_struct *> &graph_array) {
  std::copy(position_, position_ + 3, g_.r);
  g_.color = color_ + label_ * 0.1 * M_PI * M_PI;
  g_.diameter = 2;
  g_.length = 0;
  g_.draw = draw_type::fixed;
  graph_array.push_back(&g_);
}

std::pair<double, double> Cluster::GetMeanSqrDistance() {
  int size = GetSize();
  if (size == 0) {
    return std::make_pair(-1, -1);
  }
  double mean_sqr_dist = 0;
  double mean_sqr_dist_sqr = 0;
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    const double dist_sqr = GetSqrDistanceFilament((*it));
    mean_sqr_dist += dist_sqr;
    mean_sqr_dist_sqr += dist_sqr * dist_sqr;
  }
  mean_sqr_dist_sqr /= size;
  mean_sqr_dist /= size;
  return std::make_pair(mean_sqr_dist,
                        sqrt((mean_sqr_dist_sqr - SQR(mean_sqr_dist)) / size));
}

void CurvatureClusterAnalysis::InitOutput() {
  SetAnalysisName("curve_cluster");
  Analysis::InitOutput();
}
void CurvatureClusterAnalysis::InitAnalysis() {
  n_dim_ = space_->n_dim;
  Cluster::SetNDim(n_dim_);
  Cluster::SetColor(sparams_->color);
  curvature_centers_ = new double *[n_members_];
  cc_vectors_ = new double *[n_members_];
  curvature_radii_ = new double[n_members_];
  handedness_ = new int[n_members_];
  std::fill(curvature_radii_, curvature_radii_ + 3, 0.0);
  std::fill(handedness_, handedness_ + 3, 0.0);
  for (int i = 0; i < n_members_; ++i) {
    curvature_centers_[i] = new double[3];
    cc_vectors_[i] = new double[3];
    std::fill(curvature_centers_[i], curvature_centers_[i] + 3, 0.0);
    std::fill(cc_vectors_[i], cc_vectors_[i] + 3, 0.0);
  }
  Cluster::InitDataStructures(curvature_centers_, cc_vectors_, curvature_radii_,
                              handedness_);
  /* Requires MinimumDistance to be initialized */
  cluster_by_handedness_ = sparams_->cluster_by_handedness;
  RequireInteractionAnalysis();
  output_
      << "time cluster_label n_filaments pos_x pos_y pos_z avg_radius"
         " max_radius handedness mean_sqr_distance mean_sqr_distance_stderr\n";
}
void CurvatureClusterAnalysis::RunAnalysis() {
  if (params_->i_step < sparams_->n_equil) {
    return;
  }
  double dr[3] = {0};
  double mid[3] = {0};
  double dr_mag2 = 0;
  double avg_pos[3] = {0};
  MinimumDistance mindist;
  for (int i = 0; i < n_members_; ++i) {
    curvature_radii_[i] =
        (*members_)[i].GetCenterOfCurvature(curvature_centers_[i]);
    (*members_)[i].ClearCluster();
    handedness_[i] = (*members_)[i].GetHandedness();
    (*members_)[i].GetAvgScaledPosition(avg_pos);
    mindist.PointPoint(avg_pos, avg_pos, curvature_centers_[i],
                       curvature_centers_[i], cc_vectors_[i], &dr_mag2, mid);
  }
  for (int i = 0; i < n_members_ - 1; ++i) {
    double *pos_i = curvature_centers_[i];
    int hand = handedness_[i];
    for (int j = i + 1; j < n_members_; ++j) {
      // Only consider clustering filaments with the same handedness for now
      if (cluster_by_handedness_ && handedness_[j] != hand)
        continue;
      /* Check if their centers of curvature are within the largest radii of
         curvature */
      double *pos_j = curvature_centers_[j];
      double rad = MAX(curvature_radii_[i], curvature_radii_[j]);
      mindist.PointPoint(pos_i, pos_i, pos_j, pos_j, dr, &dr_mag2, mid);
      if (dr_mag2 < SQR(rad)) {
        /* If so, add them to the same cluster */
        ClusterFilaments(i, j);
      }
    }
  }
  CalculateClusterPositions();
  CheckNoCluster();
  DeleteEmptyClusters();
  GetClusterOutputs();
  output_ << std::flush;
}

void CurvatureClusterAnalysis::EndAnalysis() {
  for (int i = 0; i < n_members_; ++i) {
    delete[] curvature_centers_[i];
  }
  delete[] curvature_centers_;
  delete[] curvature_radii_;
}

void CurvatureClusterAnalysis::CreateNewCluster(int i, int j) {
  int label = cluster_label_++;
  clusters_[label] = Cluster();
  clusters_[label].SetLabel(label);
  clusters_[label].AddMember(i);
  clusters_[label].AddMember(j);
  (*members_)[i].AssignCluster(label);
  (*members_)[j].AssignCluster(label);
  clusters_[label].CalcPosition();
  if (debug_) {
    output_ << time_ << " Debug: Filaments " << i << " and " << j
            << " formed new cluster " << label << "\n";
  }
}

void CurvatureClusterAnalysis::ClusterFilaments(int i, int j) {
  const int prev_ci = (*members_)[i].GetPrevCluster();
  const int prev_cj = (*members_)[j].GetPrevCluster();
  const int ci = (*members_)[i].GetCluster();
  const int cj = (*members_)[j].GetCluster();

  if (prev_ci == 0 && prev_cj == 0 && ci == 0 && cj == 0) {
    /* If they're not clustered already, create a new cluster centered at the
       average of their centers of curvature */
    CreateNewCluster(i, j);
  } else if (prev_ci == 0 && ci == 0) {
    /* Particle i was not already clustered, so add to particle j's cluster */
    int cluster_j = (cj > 0 ? cj : prev_cj);
    clusters_[cluster_j].AddMember(i);
    (*members_)[i].AssignCluster(cluster_j);
    if (debug_) {
      output_ << time_ << " Debug: Filament " << i << " joined cluster "
              << cluster_j << "\n";
    }
  } else if (prev_cj == 0 && cj == 0) {
    /* Particle j was not already clustered, so add to particle i's cluster */
    int cluster_i = (ci > 0 ? ci : prev_ci);
    clusters_[cluster_i].AddMember(j);
    (*members_)[j].AssignCluster(cluster_i);
    if (debug_) {
      output_ << time_ << " Debug: Filament " << j << " joined cluster "
              << cluster_i << "\n";
    }
  } else {
    /* Both are/were in clusters. */
    int cluster_i = (ci > 0 ? ci : prev_ci);
    int cluster_j = (cj > 0 ? cj : prev_cj);
    /* Check if they're already in the same cluster */
    if (cluster_i == cluster_j) {
      // Check that we aren't just mindlessly reassigning them
      if (cluster_i == prev_ci && cluster_j == prev_cj) {
        double dr_i = clusters_[cluster_i].GetSqrDistanceFilament(i);
        double dr_j = clusters_[cluster_j].GetSqrDistanceFilament(j);
        double cluster_radius = clusters_[cluster_i].GetRadius();
        if (dr_i > SQR(cluster_radius) && dr_j > SQR(cluster_radius)) {
          clusters_[cluster_i].RemoveMember(i);
          clusters_[cluster_j].RemoveMember(j);
          if (debug_) {
            output_ << time_ << " Debug: Filaments " << i << " and " << j
                    << " left cluster " << cluster_i
                    << " to form a new cluster\n";
          }
          CreateNewCluster(i, j);
          return;
        }
        if (dr_i < SQR(cluster_radius)) {
          (*members_)[i].AssignCluster(cluster_i);
        }
        if (dr_j < SQR(cluster_radius)) {
          (*members_)[j].AssignCluster(cluster_j);
        }
      }
      return;
    }
    /* Otherwise, check the distance between each of their centers of
       curvature and the centers of their respective clusters and switch only
       if the other cluster is closer */
    const double drij = clusters_[cluster_i].GetSqrDistanceFilament(j);
    const double drjj = clusters_[cluster_j].GetSqrDistanceFilament(j);
    if (drij < drjj) {
      /* Filament j switches to same cluster as filament i */
      clusters_[cluster_i].AddMember(j);
      clusters_[cluster_j].RemoveMember(j);
      (*members_)[j].AssignCluster(cluster_i);
      if (debug_) {
        output_ << time_ << " Debug: Filament " << j << " left cluster "
                << cluster_j << " to join cluster " << cluster_i << "\n";
      }
    } else {
      (*members_)[j].AssignCluster(cluster_j);
    }
    double drji = clusters_[cluster_j].GetSqrDistanceFilament(i);
    double drii = clusters_[cluster_i].GetSqrDistanceFilament(i);
    if (drji < drii) {
      /* Filament i switches to same cluster as filament j */
      clusters_[cluster_j].AddMember(i);
      clusters_[cluster_i].RemoveMember(i);
      (*members_)[i].AssignCluster(cluster_j);
      if (debug_) {
        output_ << time_ << " Debug: Filament " << i << " left cluster "
                << cluster_i << " to join cluster " << cluster_j << "\n";
      }
    } else {
      (*members_)[i].AssignCluster(cluster_i);
    }
  }
}

void CurvatureClusterAnalysis::DeleteEmptyClusters() {
  for (auto it = clusters_.begin(); it != clusters_.end();) {
    if (it->second.GetSize() == 1) {
      if (debug_) {
        output_ << time_ << " Debug: Deleting size-one cluster "
                << it->second.GetLabel() << " with lone filament  "
                << it->second.GetLastMember() << "\n";
      }
      (*members_)[it->second.GetLastMember()].AssignCluster(0);
      it = clusters_.erase(it);
    } else if (it->second.GetSize() == 0) {
      if (debug_) {
        output_ << time_ << " Debug: Deleting empty cluster "
                << it->second.GetLabel() << "\n";
      }
      it = clusters_.erase(it);
    } else {
      ++it;
    }
  }
}
void CurvatureClusterAnalysis::CheckNoCluster() {
  /* Find any filaments that are no longer in a cluster and remove them from
     any previous clusters */
  for (int i = 0; i < n_members_; ++i) {
    if ((*members_)[i].GetCluster() == 0) {
      const int prev_cluster = (*members_)[i].GetPrevCluster();
      if (prev_cluster > 0) {
        if (clusters_[prev_cluster].CheckInCluster(i)) {
          (*members_)[i].AssignCluster(prev_cluster);
        } else {
          if (debug_) {
            output_ << time_ << " Debug: Filament " << i << " left cluster "
                    << prev_cluster << "\n";
          }
          clusters_[prev_cluster].RemoveMember(i);
        }
      }
    } else {
      const int cluster = (*members_)[i].GetCluster();
      if (!clusters_[cluster].CheckInCluster(i)) {
        if (debug_) {
          output_ << time_ << " Debug: Filament " << i
                  << " got too far away from cluster " << cluster << "\n";
        }
        (*members_)[i].AssignCluster(0);
        clusters_[cluster].RemoveMember(i);
      }
    }
  }
}

void CurvatureClusterAnalysis::CalculateClusterPositions() {
  for (auto it = clusters_.begin(); it != clusters_.end(); ++it) {
    it->second.CalcPosition();
  }
  CheckClusterMerge();
}
void CurvatureClusterAnalysis::GetClusterOutputs() {
  for (auto it = clusters_.begin(); it != clusters_.end(); ++it) {
    if (it->second.GetLifetime() > sparams_->cluster_lifetime_min) {
      const double *const pos = it->second.GetPosition();
      std::pair<double, double> msdist = it->second.GetMeanSqrDistance();
      output_ << time_ << " " << it->second.GetLabel() << " "
              << it->second.GetSize() << " " << pos[0] << " " << pos[1] << " "
              << pos[2] << " " << it->second.GetAvgRadius() << " "
              << it->second.GetRadius() << " " << it->second.GetHandedness()
              << " " << msdist.first << " " << msdist.second << "\n";
    }
  }
}

void CurvatureClusterAnalysis::CheckClusterMerge() {
  MinimumDistance mindist;
  double dr[3] = {0};
  double mid[3] = {0};
  std::unordered_map<int, Cluster>::iterator it;
  std::unordered_map<int, Cluster>::iterator jt;
  std::unordered_map<int, Cluster>::iterator endm1 = clusters_.end();
  std::advance(endm1, -1);
  for (it = clusters_.begin(); it != endm1; ++it) {
    const double *const pos_i = it->second.GetPosition();
    double rad_i = it->second.GetAvgRadius();
    rad_i = SQR(rad_i);
    const int hand = it->second.GetHandedness();
    for (jt = it, std::advance(jt, 1); jt != clusters_.end(); ++jt) {
      if (cluster_by_handedness_ && jt->second.GetHandedness() != hand) {
        continue;
      }
      double dr_mag2 = 0;
      const double *const pos_j = jt->second.GetPosition();
      mindist.PointPoint(pos_i, pos_i, pos_j, pos_j, dr, &dr_mag2, mid);
      double rad_j = jt->second.GetAvgRadius();
      if (dr_mag2 < rad_i || dr_mag2 < SQR(rad_j)) {
        if (it->second.GetSize() > jt->second.GetSize()) {
          if (debug_) {
            output_ << time_ << " Debug: Merging cluster "
                    << jt->second.GetLabel() << " into cluster "
                    << it->second.GetLabel() << "\n";
          }
          const std::set<int> cluster_members = jt->second.GetMembers();
          const int label = it->second.GetLabel();
          for (auto kt = cluster_members.begin(); kt != cluster_members.end();
               ++kt) {
            (*members_)[*kt].AssignCluster(label);
          }
          it->second.Merge(jt->second);
        } else {
          if (debug_) {
            output_ << time_ << " Debug: Merging cluster "
                    << it->second.GetLabel() << " into cluster "
                    << jt->second.GetLabel() << "\n";
          }
          const std::set<int> cluster_members = it->second.GetMembers();
          const int label = jt->second.GetLabel();
          for (auto kt = cluster_members.begin(); kt != cluster_members.end();
               ++kt) {
            (*members_)[*kt].AssignCluster(label);
          }
          jt->second.Merge(it->second);
        }
      }
    }
  }
}

void CurvatureClusterAnalysis::Draw(std::vector<graph_struct *> &graph_array) {
  for (auto it = clusters_.begin(); it != clusters_.end(); ++it) {
    it->second.Draw(graph_array);
  }
}
