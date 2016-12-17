#include "output_manager.h"

OutputManager::OutputManager(){}

void OutputManager::Init(system_parameters *params, 
    std::vector<graph_struct*> *graph_array, int *i_step, std::string run_name) {
  run_name_ = run_name;
  params_ = params;
  i_step_ = i_step;
  graph_array_ = graph_array;
  node_ = YAML::LoadFile(params_->datafile);
  Clear();
}

void OutputManager::WriteOutputs(){
  CalcOutputs();
  if ( (*i_step_ % params_->n_posit == 0) && params_->posit_flag) 
    WriteSpeciesPosits();

  WriteTotalThermo();
  WriteSpeciesThermo();
}

void OutputManager::WriteSpeciesPosits(){
  for (auto spec_it : species_)
    spec_it.second->WritePosits();
}

void OutputManager::AddSpecie(SpeciesBase *spec){
  SID sid = spec->GetSID();
  species_[sid] = spec;
  if (params_->posit_flag == 1){
    spec->InitOutputFile(run_name_);
  }
}

void OutputManager::AddSpecies(std::vector<SpeciesBase*> *species){
  for (auto spec_it : (*species)) 
    AddSpecie(spec_it);

  if (IsMovie())
    InitPositInput();
}

void OutputManager::Close() {
  for (auto spec : species_){
    spec.second->Close();
  }
  if (thermo_file_.is_open()) thermo_file_.close();
}

void OutputManager::Clear() {
  for (auto spec : species_)
    spec.second->ClearThermo();

  //Clear virial tensor
  memset(tot_virial_, 0, sizeof(tot_virial_));
  std::fill(tot_direct_, tot_direct_+3, 0);
  std::fill(tot_pol_direct_, tot_pol_direct_+3, 0);
  tot_energy_ =0;
}

//TODO Put in safe guard to make sure species that are not in configure file are not run
//TODO Have better initializers for each of this
void OutputManager::InitPositInput() {
  //New species map so only those with posit files will remain
  std::map<SID, SpeciesBase*> specs;

  for (auto pos_it : posit_files_){
    int nchar, n_posit, n_steps;
    std::fstream ip;

    //Open binary file of all the positions
    ip.open(pos_it, std::ios::binary | std::ios::in );
    if (!ip.is_open()){
      std::cout<<"Input "<< pos_it <<" file did not open\n";
      exit(1);
    }
    else{
      //Get Sim data from posit file
      ip.read(reinterpret_cast<char*>(&nchar), sizeof(int));
      std::string sid_str(nchar, ' ');
      ip.read(&sid_str[0], nchar);

      SID sid_posit = StringToSID(sid_str);
      posit_spec_[sid_posit];

      ip.read(reinterpret_cast<char*>(&(n_steps)), sizeof(int));
      ip.read(reinterpret_cast<char*>(&(posit_spec_[sid_posit])), sizeof(int));

      //Need n_steps to be the lowest
      if (params_->n_steps > n_steps) params_->n_steps = n_steps;
      //Will only graph as the largest n_posit 
      if (params_->n_graph < posit_spec_[sid_posit]) 
        params_->n_graph = posit_spec_[sid_posit];

      //Get the location where species data starts and close input for now
      std::ios::streampos beg = ip.tellg();
      ip.close();

      //Find species that correlates to SID and add it to specs
      for (auto spec_it : species_ )
        if (sid_posit == spec_it.first){
          spec_it.second->InitInputFile(pos_it, beg);
          specs[sid_posit] = spec_it.second;
          break;
        }
    }
  }
  //Replace species map
  species_ = specs; 
}

void OutputManager::ReadSpeciesPositions() {

  for( auto ps_it : posit_spec_) 
    if (*i_step_ % ps_it.second == 0 )
      species_[ps_it.first]->ReadPosits();

  if (*i_step_ % params_->n_graph == 0)
    GetGraphicsStructure();
}

void OutputManager::GetGraphicsStructure() {
  graph_array_->clear();
  for( auto spec : species_ )
    spec.second->Draw(graph_array_);
  //TODO Get uegine interactions to graph
}

void OutputManager::CalcOutputs() {
  if ( YAML::Node thermo = node_["thermo"] ){
    if (*i_step_ % thermo["steps"].as<int>() == 0 ) {
      YAML::Node thermo_spec = node_["thermo"]["species"];
      SID sid;
      for (yaml_it it=thermo_spec.begin(); it!=thermo_spec.end(); ++it){
        if (it->first.as<std::string>().compare("all") == 0){
          for(auto spec : species_)
            CalcSpeciesOutputs(spec.first, &(it->second));
        }
      }
    }
  }
}

void OutputManager::CalcSpeciesOutputs(SID sid, YAML::Node *data_list) {
  int n_dim = params_->n_dim;
  double const *spec_direct;
  
  for (yaml_it data_it=data_list->begin(); data_it!=data_list->end(); data_it++){
    if (data_it->as<std::string>().compare("virial") == 0){
      species_[sid]->GetVirial(tot_virial_);
    }
    else if (data_it->as<std::string>().compare("director") == 0){
      spec_direct = species_[sid]->GetDirector();
      for (int i=0; i<n_dim; i++)
        tot_direct_[i] += spec_direct[i];
      normalize_vector(tot_direct_, n_dim);
    }
    else if (data_it->as<std::string>().compare("polar_director") == 0){
      spec_direct = species_[sid]->GetPolarDirector();
      for (int i=0; i<n_dim; i++)
        tot_pol_direct_[i] += spec_direct[i];
    }
    else if (data_it->as<std::string>().compare("energy") == 0){
      tot_energy_ += species_[sid]->GetTotalEnergy();
    }
    //TODO Potential energy
    //TODO Kinetic energy
    //TODO Distributions
    //TODO Correlation functions
  }
}

void OutputManager::MakeHeaders() {
  if(node_["thermo"]){
    YAML::Node data_list = node_["thermo"]["species"];
    if(YAML::Node all = data_list["all"]){
      std::cout<<run_name_<<std::endl;
      thermo_file_.open(run_name_ + ".thermo", std::ios::out);
      thermo_file_<<"Time";
      for (yaml_it it=all.begin(); it!=all.end(); ++it){
        if (it->as<std::string>().compare("virial") == 0){
          for (int i=0; i<params_->n_dim; i++)
            for (int j=0; j<params_->n_dim; j++)
              thermo_file_<<" sig_"<<i<<j;
        }
        else if (it->as<std::string>().compare("director") == 0)
          for (int i=0; i<params_->n_dim; i++)
              thermo_file_<<" tot_n_"<<i;
        else if (it->as<std::string>().compare("polar_director") == 0)
          for (int i=0; i<params_->n_dim; i++)
              thermo_file_<<" tot_pn_"<<i;
        else
          thermo_file_<<" tot_"<<it->as<std::string>();
      }
      thermo_file_<<"\n";
    }
    //TODO Set up out files for other variables
  }
  //TODO Add other forms of output
}

void OutputManager::WriteTotalThermo() {
  if(YAML::Node thermo = node_["thermo"]){
    if (*i_step_ % thermo["steps"].as<int>() == 0 ) {
      YAML::Node data_list = thermo["species"];
      thermo_file_<<params_->delta * (*i_step_);
      for (yaml_it it=data_list["all"].begin(); it!=data_list["all"].end(); ++it){
        if (it->as<std::string>().compare("virial") == 0)
          for (int i=0; i<params_->n_dim; i++)
            for (int j=0; j<params_->n_dim; j++)
              thermo_file_<<" "<<tot_virial_[3*i+j];
        else if (it->as<std::string>().compare("director") == 0)
          for (int i=0; i<params_->n_dim; i++)
              thermo_file_<<" "<<tot_direct_[i];
        else if (it->as<std::string>().compare("polar_director") == 0)
          for (int i=0; i<params_->n_dim; i++)
              thermo_file_<<" "<<tot_pol_direct_[i];
        else if (it->as<std::string>() == "energy")
          thermo_file_<<" "<<tot_energy_;
      }
      thermo_file_<<"\n";
    }
  }
  Clear();
}

void OutputManager::WriteSpeciesThermo() {}

