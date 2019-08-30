#include "rigid_mesh.hpp"

void RigidMesh::SetInteractorLength(double const len) { ix_length_ = len; }
void RigidMesh::SetInteractorDiameter(double const d) { ix_diameter_ = d; }
void RigidMesh::SetInteractorPosition(double const* const pos) {
  std::copy(pos, pos + 3, ix_position_);
}
void RigidMesh::SetInteractorPrevPosition(double const* const ppos) {
  std::copy(ppos, ppos + 3, ix_prev_position_);
}
void RigidMesh::SetInteractorScaledPosition(double const* const spos) {
  std::copy(spos, spos + 3, ix_scaled_position_);
}
void RigidMesh::SetInteractorOrientation(double const* const u) {
  std::copy(u, u + 3, ix_orientation_);
}
double const RigidMesh::GetInteractorLength() { return ix_length_; }
double const RigidMesh::GetInteractorDiameter() { return ix_diameter_; }
double const* const RigidMesh::GetInteractorPosition() { return ix_position_; }
double const* const RigidMesh::GetInteractorPrevPosition() {
  return ix_prev_position_;
}
double const* const RigidMesh::GetInteractorScaledPosition() {
  return ix_scaled_position_;
}
double const* const RigidMesh::GetInteractorOrientation() {
  return ix_orientation_;
}
void RigidMesh::Draw(std::vector<graph_struct*>* graph_array) {
  Mesh::Draw(graph_array);
}
