#ifndef _SIMCORE_FILAMENT_CURVATURE_CLUSTER_ANALYSIS_H_
#define _SIMCORE_FILAMENT_CURVATURE_CLUSTER_ANALYSIS_H_

#include "analysis.hpp"
#include "cell_list.hpp"
#include "filament.hpp"
#include "minimum_distance.hpp"
#include <unordered_map>

class Cluster {
private:
  static int n_dim_;
  static double **curvature_centers_;
  static double **cc_vectors_;
  static double *curvature_radii_;
  static int *handedness_;
  static double color_;

  graph_struct g_;
  std::set<int> members_;
  double position_[3] = {0};
  double prev_position_[3] = {0};
  int label_ = -1;
  int lifetime_ = 0;
  double radius_ = 0;
  double avg_radius_ = 0;
  int avg_handedness_ = 0;

public:
  Cluster() = default;
  static void SetNDim(int ndim) { n_dim_ = ndim; }
  static void SetColor(double c) { color_ = c; }
  static void InitDataStructures(double **c_centers, double **cc_vec,
                                 double *c_radii, int *hand) {
    curvature_centers_ = c_centers;
    cc_vectors_ = cc_vec;
    curvature_radii_ = c_radii;
    handedness_ = hand;
  }
  void AddMember(int i) {
    members_.insert(i);
    if (curvature_radii_[i] > radius_) {
      radius_ = curvature_radii_[i];
    }
  }
  void RemoveMember(int i) { members_.erase(members_.find(i)); }
  const int GetSize() { return members_.size(); }
  void SetLabel(int label) { label_ = label; }
  const int GetLabel() { return label_; }
  void CalcPosition();
  const double GetSqrDistanceFilament(int i);
  const double GetSqrDistance(double *pos);
  const int GetLastMember() { return *(members_.begin()); }
  const int GetLifetime() { return ++lifetime_; }
  const double *const GetPosition() { return position_; }
  const double GetRadius() { return radius_; }
  const double GetAvgRadius() { return avg_radius_; }
  const int GetHandedness() { return avg_handedness_; }
  void Merge(Cluster &other);
  const bool CheckInCluster(int i);
  void Draw(std::vector<graph_struct *> &graph_array);
  std::pair<double, double> GetMeanSqrDistance();
  const std::set<int> &GetMembers() const { return members_; }
};

class CurvatureClusterAnalysis
    : public Analysis<Filament, species_id::filament> {

protected:
  bool debug_ = false; // For verbose debugging of outputs
  std::unordered_map<int, Cluster> clusters_;
  bool cluster_by_handedness_;
  double **curvature_centers_;
  double *curvature_radii_;
  double **cc_vectors_;
  int *handedness_;
  int cluster_label_ = 1;
  int n_dim_ = -1;

  void InitOutput();
  void InitAnalysis();
  void RunAnalysis();
  void EndAnalysis();
  void ClusterFilaments(int i, int j);
  void CreateNewCluster(int i, int j);
  void DeleteEmptyClusters();
  void CheckNoCluster();
  void CalculateClusterPositions();
  void CheckClusterMerge();
  void GetClusterOutputs();

public:
  void Draw(std::vector<graph_struct *> &graph_array);
};

#endif // _SIMCORE_FILAMENT_CURVATURE_CLUSTER_ANALYSIS_H_
