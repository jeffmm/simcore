#include "spindle_pole_body.h"

// SPB specifics
void SpindlePoleBody::InitConfigurator(const double r,
                                       const double theta,
                                       const double phi,
                                       const double diameter,
                                       const double attach_diameter) {
  diameter_ = diameter;
  attach_diameter_ = attach_diameter;
  std::cout << "Inserting SPB at\n" << std::setprecision(16)
    << "\tr: " << r << ", "
    << "theta: " << theta << ", "
    << "phi: " << phi << "\n";
  double r_anchor[3] = {0.0, 0.0, 0.0};
  r_anchor[0] = r * sin(theta) * cos(phi);
  r_anchor[1] = r * sin(theta) * sin(phi);
  r_anchor[2] = r * cos(theta);
  SetPosition(r_anchor);
  SetPrevPosition(r_anchor);

  UpdatePeriodic();
}



// Species specifics
void SpindlePoleBodySpecies::Configurator() {
  char *filename = params_->config_file;
  std::cout << "SpindlePoleBody species\n";

  YAML::Node node = YAML::LoadFile(filename);

  std::cout << " Generic Properties:\n";
  std::string insertion_type;
  insertion_type = node["spb"]["properties"]["insertion_type"].as<std::string>();
  std::cout << "   insertion type: " << insertion_type << std::endl;
  bool can_overlap = node["spb"]["properties"]["overlap"].as<bool>();
  std::cout << "   can overlap:    " << (can_overlap ? "true" : "false") << std::endl;

  // Coloring
  double color[4] = {1.0, 0.0, 0.0, 1.0};
  int draw_type = 0; // default to orientation
  if (node["spb"]["properties"]["color"]) {
    for (int i = 0; i < 4; ++i) {
      color[i] = node["spb"]["properties"]["color"][i].as<double>();
    }
    std::cout << "   color: [" << color[0] << ", " << color[1] << ", " << color[2] << ", "
      << color[3] << "]\n";
  }
  if (node["spb"]["properties"]["draw_type"]) {
    std::string draw_type_s = node["spb"]["properties"]["draw_type"].as<std::string>();
    std::cout << "   draw_type: " << draw_type_s << std::endl;
    if (draw_type_s.compare("flat") == 0) {
      draw_type = 0;
    }
  }

  if (insertion_type.compare("rtp") == 0) {
    if (!can_overlap) {
      std::cout << "Warning, location insertion overrides overlap\n";
      can_overlap = true;
    }
    int nspbs = (int)node["spb"]["spb"].size();
    std::cout << "   nspbs: " << nspbs << std::endl;
    for (int ispb = 0; ispb < nspbs; ++ispb) {
      double diameter         = node["spb"]["spb"][ispb]["diameter"].as<double>();
      double attach_diameter  = node["spb"]["spb"][ispb]["attach_diameter"].as<double>();
      double spb_theta        = node["spb"]["spb"][ispb]["theta"].as<double>();
      double spb_phi          = node["spb"]["spb"][ispb]["phi"].as<double>();
      SpindlePoleBody *member = new SpindlePoleBody(params_, space_, gsl_rng_get(rng_.r), GetSID());
      member->InitConfigurator(0.5 * space_->unit_cell[0][0], spb_theta, spb_phi, diameter, attach_diameter);
      member->SetColor(color, draw_type);
      member->Dump();
      members_.push_back(member);
    }
  } else {
    std::cout << "Insertion type " << insertion_type << " not implemented yet\n";
    exit(1);
  }
}
