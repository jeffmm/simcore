#include "bond.h"
void Bond::Init() {}

void Bond::Draw(std::vector<graph_struct*> * graph_array) {
  for (int i=0; i<space_->n_periodic; ++i)
    g_.r[i] = scaled_position_[i];
  for (int i=space_->n_periodic; i<n_dim_; ++i)
    g_.r[i] = position_[i];
  std::copy(orientation_, orientation_+3, g_.u);
  g_.color = color_;
  //g_.length = (length_-diameter_ > 0 ? length_-diameter_ : 0);
  g_.length = length_;
  if (graph_diameter_ == 0) 
    g_.diameter = diameter_;
  else 
    g_.diameter = graph_diameter_;
  g_.draw_type = draw_type_;
  graph_array->push_back(&g_);
}
