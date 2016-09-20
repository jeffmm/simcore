#include "output_manager.h"

OutputManager::OutputManager(){}

void OutputManager::Init(system_parameters *params, 
    std::vector<graph_struct*> *graph_array, int *i_step) {
  params_ = params;
  i_step_ = i_step;
  graph_array_ = graph_array;
}

void OutputManager::WriteOutputs(){
  if ( (*i_step_ % params_->n_posit == 0) && params_->posit_flag)  {
    WriteSpeciesPosits();
  //Write observable output flags
    //TODO WriteVirial
  }
}

void OutputManager::WriteSpeciesPosits(){
  for (auto spec_it : species_)
    spec_it.second->WritePosits();
}

void OutputManager::AddSpecie(SpeciesBase *spec){
  SID sid = spec->GetSID();
  species_[sid] = spec;
  if (params_->posit_flag == 1){
    spec->InitOutputFile();
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
}

void OutputManager::Clear() {}

void OutputManager::CalcVirial() {}

void OutputManager::WriteVirial() {}

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








