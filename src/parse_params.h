#ifndef _SIMCORE_PARSE_PARAMS_H_
#define _SIMCORE_PARSE_PARAMS_H_

#include <iostream>
#include <string>
#include "yaml-cpp/yaml.h"
#include "parameters.h"

void parse_params(YAML::Node node, system_parameters *params) {

  std::string param_name, param_value, struct_name;
  for (YAML::const_iterator it=node.begin(); it!= node.end(); ++it) {
    if (it->second.IsMap()) {
      struct_name = it->first.as<std::string>();
      if (false) {}
      else if (struct_name.compare("species") == 0) {
        for (YAML::const_iterator jt=it->second.begin(); jt!= it->second.end(); ++jt) {
          param_name = jt->first.as<std::string>();
          if (false) {}
          else if (param_name.compare("num")==0) {
            params->species.num = jt->second.as<int>();
          }
          else if (param_name.compare("insertion_type")==0) {
            params->species.insertion_type = jt->second.as<std::string>();
          }
          else if (param_name.compare("overlap")==0) {
            params->species.overlap = jt->second.as<int>();
          }
          else if (param_name.compare("draw_type")==0) {
            params->species.draw_type = jt->second.as<std::string>();
          }
          else if (param_name.compare("color")==0) {
            params->species.color = jt->second.as<double>();
          }
          else if (param_name.compare("posit_flag")==0) {
            params->species.posit_flag = jt->second.as<int>();
          }
          else if (param_name.compare("spec_flag")==0) {
            params->species.spec_flag = jt->second.as<int>();
          }
          else if (param_name.compare("checkpoint_flag")==0) {
            params->species.checkpoint_flag = jt->second.as<int>();
          }
          else if (param_name.compare("n_posit")==0) {
            params->species.n_posit = jt->second.as<int>();
          }
          else if (param_name.compare("n_spec")==0) {
            params->species.n_spec = jt->second.as<int>();
          }
          else if (param_name.compare("n_checkpoint")==0) {
            params->species.n_checkpoint = jt->second.as<int>();
          }
          else {
            std::cout << "  WARNING! Unrecognized " << struct_name <<" parameter: '" << param_name << "'\n";
          }
        }
      }
      else if (struct_name.compare("filament") == 0) {
        for (YAML::const_iterator jt=it->second.begin(); jt!= it->second.end(); ++jt) {
          param_name = jt->first.as<std::string>();
          if (false) {}
          else if (param_name.compare("diameter")==0) {
            params->filament.diameter = jt->second.as<double>();
          }
          else if (param_name.compare("length")==0) {
            params->filament.length = jt->second.as<double>();
          }
          else if (param_name.compare("persistence_length")==0) {
            params->filament.persistence_length = jt->second.as<double>();
          }
          else if (param_name.compare("max_length")==0) {
            params->filament.max_length = jt->second.as<double>();
          }
          else if (param_name.compare("min_length")==0) {
            params->filament.min_length = jt->second.as<double>();
          }
          else if (param_name.compare("max_child_length")==0) {
            params->filament.max_child_length = jt->second.as<double>();
          }
          else if (param_name.compare("spiral_flag")==0) {
            params->filament.spiral_flag = jt->second.as<int>();
          }
          else if (param_name.compare("driving_factor")==0) {
            params->filament.driving_factor = jt->second.as<double>();
          }
          else if (param_name.compare("friction_ratio")==0) {
            params->filament.friction_ratio = jt->second.as<double>();
          }
          else if (param_name.compare("dynamic_instability_flag")==0) {
            params->filament.dynamic_instability_flag = jt->second.as<int>();
          }
          else if (param_name.compare("force_induced_catastrophe_flag")==0) {
            params->filament.force_induced_catastrophe_flag = jt->second.as<int>();
          }
          else if (param_name.compare("f_shrink_to_grow")==0) {
            params->filament.f_shrink_to_grow = jt->second.as<double>();
          }
          else if (param_name.compare("f_shrink_to_pause")==0) {
            params->filament.f_shrink_to_pause = jt->second.as<double>();
          }
          else if (param_name.compare("f_pause_to_grow")==0) {
            params->filament.f_pause_to_grow = jt->second.as<double>();
          }
          else if (param_name.compare("f_pause_to_shrink")==0) {
            params->filament.f_pause_to_shrink = jt->second.as<double>();
          }
          else if (param_name.compare("f_grow_to_pause")==0) {
            params->filament.f_grow_to_pause = jt->second.as<double>();
          }
          else if (param_name.compare("f_grow_to_shrink")==0) {
            params->filament.f_grow_to_shrink = jt->second.as<double>();
          }
          else if (param_name.compare("metric_forces")==0) {
            params->filament.metric_forces = jt->second.as<int>();
          }
          else if (param_name.compare("v_poly")==0) {
            params->filament.v_poly = jt->second.as<double>();
          }
          else if (param_name.compare("v_depoly")==0) {
            params->filament.v_depoly = jt->second.as<double>();
          }
          else if (param_name.compare("theta_analysis")==0) {
            params->filament.theta_analysis = jt->second.as<int>();
          }
          else if (param_name.compare("num")==0) {
            params->filament.num = jt->second.as<int>();
          }
          else if (param_name.compare("insertion_type")==0) {
            params->filament.insertion_type = jt->second.as<std::string>();
          }
          else if (param_name.compare("overlap")==0) {
            params->filament.overlap = jt->second.as<int>();
          }
          else if (param_name.compare("draw_type")==0) {
            params->filament.draw_type = jt->second.as<std::string>();
          }
          else if (param_name.compare("color")==0) {
            params->filament.color = jt->second.as<double>();
          }
          else if (param_name.compare("posit_flag")==0) {
            params->filament.posit_flag = jt->second.as<int>();
          }
          else if (param_name.compare("spec_flag")==0) {
            params->filament.spec_flag = jt->second.as<int>();
          }
          else if (param_name.compare("checkpoint_flag")==0) {
            params->filament.checkpoint_flag = jt->second.as<int>();
          }
          else if (param_name.compare("n_posit")==0) {
            params->filament.n_posit = jt->second.as<int>();
          }
          else if (param_name.compare("n_spec")==0) {
            params->filament.n_spec = jt->second.as<int>();
          }
          else if (param_name.compare("n_checkpoint")==0) {
            params->filament.n_checkpoint = jt->second.as<int>();
          }
          else {
            std::cout << "  WARNING! Unrecognized " << struct_name <<" parameter: '" << param_name << "'\n";
          }
        }
      }
      else if (struct_name.compare("hard_rod") == 0) {
        for (YAML::const_iterator jt=it->second.begin(); jt!= it->second.end(); ++jt) {
          param_name = jt->first.as<std::string>();
          if (false) {}
          else if (param_name.compare("diameter")==0) {
            params->hard_rod.diameter = jt->second.as<double>();
          }
          else if (param_name.compare("length")==0) {
            params->hard_rod.length = jt->second.as<double>();
          }
          else if (param_name.compare("min_length")==0) {
            params->hard_rod.min_length = jt->second.as<double>();
          }
          else if (param_name.compare("max_length")==0) {
            params->hard_rod.max_length = jt->second.as<double>();
          }
          else if (param_name.compare("max_child_length")==0) {
            params->hard_rod.max_child_length = jt->second.as<double>();
          }
          else if (param_name.compare("driving_factor")==0) {
            params->hard_rod.driving_factor = jt->second.as<double>();
          }
          else if (param_name.compare("num")==0) {
            params->hard_rod.num = jt->second.as<int>();
          }
          else if (param_name.compare("insertion_type")==0) {
            params->hard_rod.insertion_type = jt->second.as<std::string>();
          }
          else if (param_name.compare("overlap")==0) {
            params->hard_rod.overlap = jt->second.as<int>();
          }
          else if (param_name.compare("draw_type")==0) {
            params->hard_rod.draw_type = jt->second.as<std::string>();
          }
          else if (param_name.compare("color")==0) {
            params->hard_rod.color = jt->second.as<double>();
          }
          else if (param_name.compare("posit_flag")==0) {
            params->hard_rod.posit_flag = jt->second.as<int>();
          }
          else if (param_name.compare("spec_flag")==0) {
            params->hard_rod.spec_flag = jt->second.as<int>();
          }
          else if (param_name.compare("checkpoint_flag")==0) {
            params->hard_rod.checkpoint_flag = jt->second.as<int>();
          }
          else if (param_name.compare("n_posit")==0) {
            params->hard_rod.n_posit = jt->second.as<int>();
          }
          else if (param_name.compare("n_spec")==0) {
            params->hard_rod.n_spec = jt->second.as<int>();
          }
          else if (param_name.compare("n_checkpoint")==0) {
            params->hard_rod.n_checkpoint = jt->second.as<int>();
          }
          else {
            std::cout << "  WARNING! Unrecognized " << struct_name <<" parameter: '" << param_name << "'\n";
          }
        }
      }
      else if (struct_name.compare("br_bead") == 0) {
        for (YAML::const_iterator jt=it->second.begin(); jt!= it->second.end(); ++jt) {
          param_name = jt->first.as<std::string>();
          if (false) {}
          else if (param_name.compare("diameter")==0) {
            params->br_bead.diameter = jt->second.as<double>();
          }
          else if (param_name.compare("driving_factor")==0) {
            params->br_bead.driving_factor = jt->second.as<double>();
          }
          else if (param_name.compare("packing_fraction")==0) {
            params->br_bead.packing_fraction = jt->second.as<double>();
          }
          else if (param_name.compare("num")==0) {
            params->br_bead.num = jt->second.as<int>();
          }
          else if (param_name.compare("insertion_type")==0) {
            params->br_bead.insertion_type = jt->second.as<std::string>();
          }
          else if (param_name.compare("overlap")==0) {
            params->br_bead.overlap = jt->second.as<int>();
          }
          else if (param_name.compare("draw_type")==0) {
            params->br_bead.draw_type = jt->second.as<std::string>();
          }
          else if (param_name.compare("color")==0) {
            params->br_bead.color = jt->second.as<double>();
          }
          else if (param_name.compare("posit_flag")==0) {
            params->br_bead.posit_flag = jt->second.as<int>();
          }
          else if (param_name.compare("spec_flag")==0) {
            params->br_bead.spec_flag = jt->second.as<int>();
          }
          else if (param_name.compare("checkpoint_flag")==0) {
            params->br_bead.checkpoint_flag = jt->second.as<int>();
          }
          else if (param_name.compare("n_posit")==0) {
            params->br_bead.n_posit = jt->second.as<int>();
          }
          else if (param_name.compare("n_spec")==0) {
            params->br_bead.n_spec = jt->second.as<int>();
          }
          else if (param_name.compare("n_checkpoint")==0) {
            params->br_bead.n_checkpoint = jt->second.as<int>();
          }
          else {
            std::cout << "  WARNING! Unrecognized " << struct_name <<" parameter: '" << param_name << "'\n";
          }
        }
      }
      else if (struct_name.compare("md_bead") == 0) {
        for (YAML::const_iterator jt=it->second.begin(); jt!= it->second.end(); ++jt) {
          param_name = jt->first.as<std::string>();
          if (false) {}
          else if (param_name.compare("diameter")==0) {
            params->md_bead.diameter = jt->second.as<double>();
          }
          else if (param_name.compare("mass")==0) {
            params->md_bead.mass = jt->second.as<double>();
          }
          else if (param_name.compare("driving_factor")==0) {
            params->md_bead.driving_factor = jt->second.as<double>();
          }
          else if (param_name.compare("num")==0) {
            params->md_bead.num = jt->second.as<int>();
          }
          else if (param_name.compare("insertion_type")==0) {
            params->md_bead.insertion_type = jt->second.as<std::string>();
          }
          else if (param_name.compare("overlap")==0) {
            params->md_bead.overlap = jt->second.as<int>();
          }
          else if (param_name.compare("draw_type")==0) {
            params->md_bead.draw_type = jt->second.as<std::string>();
          }
          else if (param_name.compare("color")==0) {
            params->md_bead.color = jt->second.as<double>();
          }
          else if (param_name.compare("posit_flag")==0) {
            params->md_bead.posit_flag = jt->second.as<int>();
          }
          else if (param_name.compare("spec_flag")==0) {
            params->md_bead.spec_flag = jt->second.as<int>();
          }
          else if (param_name.compare("checkpoint_flag")==0) {
            params->md_bead.checkpoint_flag = jt->second.as<int>();
          }
          else if (param_name.compare("n_posit")==0) {
            params->md_bead.n_posit = jt->second.as<int>();
          }
          else if (param_name.compare("n_spec")==0) {
            params->md_bead.n_spec = jt->second.as<int>();
          }
          else if (param_name.compare("n_checkpoint")==0) {
            params->md_bead.n_checkpoint = jt->second.as<int>();
          }
          else {
            std::cout << "  WARNING! Unrecognized " << struct_name <<" parameter: '" << param_name << "'\n";
          }
        }
      }
      else {
        std::cout << "  WARNING! Unrecognized struct parameter '" << struct_name << "'\n";
      }
    }
    else {
      param_name = it->first.as<std::string>();
      if (false) {}
      else if (param_name.compare("seed")==0) {
        params->seed = it->second.as<long>();
      }
      else if (param_name.compare("n_runs")==0) {
        params->n_runs = it->second.as<int>();
      }
      else if (param_name.compare("n_random")==0) {
        params->n_random = it->second.as<int>();
      }
      else if (param_name.compare("run_name")==0) {
        params->run_name = it->second.as<std::string>();
      }
      else if (param_name.compare("n_dim")==0) {
        params->n_dim = it->second.as<int>();
      }
      else if (param_name.compare("n_periodic")==0) {
        params->n_periodic = it->second.as<int>();
      }
      else if (param_name.compare("boundary_type")==0) {
        params->boundary_type = it->second.as<int>();
      }
      else if (param_name.compare("system_radius")==0) {
        params->system_radius = it->second.as<double>();
      }
      else if (param_name.compare("n_steps")==0) {
        params->n_steps = it->second.as<int>();
      }
      else if (param_name.compare("delta")==0) {
        params->delta = it->second.as<double>();
      }
      else if (param_name.compare("cell_length")==0) {
        params->cell_length = it->second.as<double>();
      }
      else if (param_name.compare("n_update_cells")==0) {
        params->n_update_cells = it->second.as<int>();
      }
      else if (param_name.compare("graph_flag")==0) {
        params->graph_flag = it->second.as<int>();
      }
      else if (param_name.compare("n_graph")==0) {
        params->n_graph = it->second.as<int>();
      }
      else if (param_name.compare("graph_diameter")==0) {
        params->graph_diameter = it->second.as<double>();
      }
      else if (param_name.compare("graph_background")==0) {
        params->graph_background = it->second.as<int>();
      }
      else if (param_name.compare("load_checkpoint")==0) {
        params->load_checkpoint = it->second.as<int>();
      }
      else if (param_name.compare("insertion_type")==0) {
        params->insertion_type = it->second.as<std::string>();
      }
      else if (param_name.compare("target_pressure")==0) {
        params->target_pressure = it->second.as<double>();
      }
      else if (param_name.compare("movie_flag")==0) {
        params->movie_flag = it->second.as<int>();
      }
      else if (param_name.compare("movie_directory")==0) {
        params->movie_directory = it->second.as<std::string>();
      }
      else if (param_name.compare("time_flag")==0) {
        params->time_flag = it->second.as<int>();
      }
      else if (param_name.compare("bud_height")==0) {
        params->bud_height = it->second.as<double>();
      }
      else if (param_name.compare("bud_radius")==0) {
        params->bud_radius = it->second.as<double>();
      }
      else if (param_name.compare("lj_epsilon")==0) {
        params->lj_epsilon = it->second.as<double>();
      }
      else if (param_name.compare("wca_eps")==0) {
        params->wca_eps = it->second.as<double>();
      }
      else if (param_name.compare("wca_sig")==0) {
        params->wca_sig = it->second.as<double>();
      }
      else if (param_name.compare("ss_a")==0) {
        params->ss_a = it->second.as<double>();
      }
      else if (param_name.compare("ss_rs")==0) {
        params->ss_rs = it->second.as<double>();
      }
      else if (param_name.compare("ss_eps")==0) {
        params->ss_eps = it->second.as<double>();
      }
      else if (param_name.compare("f_cutoff")==0) {
        params->f_cutoff = it->second.as<double>();
      }
      else if (param_name.compare("max_overlap")==0) {
        params->max_overlap = it->second.as<int>();
      }
      else if (param_name.compare("constant_volume")==0) {
        params->constant_volume = it->second.as<int>();
      }
      else if (param_name.compare("constant_pressure")==0) {
        params->constant_pressure = it->second.as<int>();
      }
      else if (param_name.compare("pressure_time")==0) {
        params->pressure_time = it->second.as<int>();
      }
      else if (param_name.compare("target_radius")==0) {
        params->target_radius = it->second.as<double>();
      }
      else if (param_name.compare("compressibility")==0) {
        params->compressibility = it->second.as<double>();
      }
      else if (param_name.compare("stoch_flag")==0) {
        params->stoch_flag = it->second.as<int>();
      }
      else if (param_name.compare("thermo_flag")==0) {
        params->thermo_flag = it->second.as<int>();
      }
      else if (param_name.compare("n_thermo")==0) {
        params->n_thermo = it->second.as<int>();
      }
      else if (param_name.compare("interaction_flag")==0) {
        params->interaction_flag = it->second.as<int>();
      }
      else {
        std::cout << "  WARNING! Unrecognized parameter '" <<  param_name << "'\n";
      }
    }
  }
}
#endif // _SIMCORE_PARSE_PARAMS_H_