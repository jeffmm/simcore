#include "bond.h"
void Bond::Init() {}

void Bond::Draw(std::vector<graph_struct*> * graph_array) {
  std::copy(position_, position_+3, g_.r);
  std::copy(orientation_, orientation_+3, g_.u);
  std::copy(color_, color_+4, g_.color);
  //g_.length = (length_-diameter_ > 0 ? length_-diameter_ : 0);
  g_.length = length_;
  if (graph_diameter_ == 0) 
    g_.diameter = diameter_;
  else 
    g_.diameter = graph_diameter_;
  g_.draw_type = draw_type_;
  graph_array->push_back(&g_);
}
