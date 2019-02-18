# simcore

A modular, object-oriented program for coarse-grained physics simulations, using something I call **SIM**ple-**C**omposite **O**bject **RE**presentation.

[![DOI](https://zenodo.org/badge/78224592.svg)](https://zenodo.org/badge/latestdoi/78224592)

Updated: 11/30/2018.

## About simcore

simcore is written in c++ and designed for coarse-grained physics simulations with modularity and scalability in mind. All objects in the simulation are representable as a composite of what I call "simple" objects (points, spheres, rigid cylinders, and 2d polygon surfaces would all qualify). For short-range interactions, simcore uses cell and neighbor lists for improved performance and OpenMP for parallelization.

Although simcore is meant to be a generalized molecular/Brownian dynamics simulation engine, thanks to the narrow focus of my PhD research, it has up until now almost exclusively been used to model semiflexible filaments, and for that reason has come closer to resembling single-purpose software. It's still quite easy, for example, to use simcore for basic molecular dynamics simulations of interacting point-like particles. Modularity is still there in the basic design, so in the future I may add more object types, but as far as pre-written object types go, _it's all about the filaments_.

![A simulation using simcore](figs/simcore_snapshot.png "A simulation using simcore")

## Installation

For successful installation, make sure you have the following libraries/binaries installed:

  * your favorite c++ compiler
  * yaml-cpp (https://github.com/jbeder/yaml-cpp)
  * openGL
  * gsl
  * glew
  * glfw3
  * openmp
  * fftw

Once appropriate compiler and library paths are defined in the Makefile, the command 'make simcore' will build the simcore binary in the local directory in debug mode. You can add options such as CFG=release, THREADING=eomp, or NOGRAPH=true, to change optimization flags, toggle OpenMP for parallelization, or disable graphics. Note that the graphics are almost guaranteed to give some trouble, since almost all the libraries are at this point deprecated (especially since Apple recently declared OpenGL deprecated entirely).

## Running simcore

The simcore binary is run with

```
./simcore --flag1 --flag2 ... params_file 
```

The following flags are available:

* --help (-h)
  * show the help menu which gives short descriptions about each of the flags as well as binary usage
* --debug (-d)
  * run the simulation while setting the debug_trace flag, which is a global flag that can be used to output debug info or run unit tests 
* --run-name rname (-r rname)
  * overwrites the parameter "run_name" with rname which serves as a prefix for all outputs 
* --n-runs num (-n num)
  * overwrites the parameter "n_runs" with num, which tells the simulation how many times to run the given parameter set with different random number generator seeds.
* --movie (-m)
  * uses the parameters file params_file to load any output files that were generated from previous runs of the simulation to replay the graphics and record the frames as bitmaps into the directory specified with the "movie_directory" parameter.
* --posit (-p)
  * makes simcore use the position (.posit) files instead of species (.spec) files to run analyses or make movies. Useful for simple object types (like spheres) that don't require more than position, orientation, length, diameter, etc, or in case you don't need fine-grained visualization of more complicated objects.
* --analysis (-a)
  * loads posit/spec files into simulation for analysis in the same manner as the movie flag.
* --reduce reduce_factor (-R reduce_factor)
  * Reads in output files and writes new output files that are smaller by a factor of reduce_factor, effectively reducing time resolution of output data.
* --load (-l)
  * specifies to load any checkpoint files corresponding to the given parameter file, which can be used to continue a simulation that ended prematurely. New simulation will be given the name old_simulation_name_reload00n where n is the number of reloads performed on that simulation.
* --with-reloads (-w)
  * If running analyses or making movies, simcore will look for parameter files that have the same run name but with the reload00n addendum and attempt to open the corresponding output files whenever it reached EOF while reading an output file.
* --blank (-b)
  * generates all relevant parameter files using the SimulationManager without running the simulations. Useful for generating many parameter files from parameter sets (discussed below) and deploying simulations on different processors and/or machines.
* --single-frame (-M)
  * Still draws the graphics to the OpenGL window, but only records a bitmap of the last frame. Useful for generating final-state snapshots without creating a large amount of bitmaps used in making movies.
* --auto-graph (-G)
  * Does not wait for user input (usually the ESC key) when drawing the simulation output using OpenGL

## Parameters

There are three parameter types, but only two are necessary: global and species. Global parameters are seen by the entire system and species parameters are unique to the specified species. There is also an optional "global species" parameter type that affects every species.

What do I mean by species? simcore assumes that any given simulation will likely have many copies of one kind of thing, which I call a species, perhaps interacting with other species of other kinds. In a system of interacting spheres, the species is 'sphere.' In a system of interacting semiflexible filaments, the species is 'filament.' Simulations can have many types of species all interacting with each other with different species-species interaction potentials.
  
The parameter file must be in the YAML file format and is set up in the following way:

```
global_param_1: gp1_value
global_param_2: gp2_value
species: # This is the global species parameter
    global_species_param_1: gsp1_value
    global_species_param_2: gsp2_value
specific_species_name:
    species_param_1: sp1_value
    species_param_2: sp2_value
```

The "global species" parameter type houses parameters that are shared by every species, such as "num" which specifies how many items of each species to insert into the system, "insertion_type" which specifies how to insert them, "posit_flag" which instructs whether to output position files, etc. These do not have to be specified using the global species parameter, however, and may instead be specified by each specific species.  

If any parameter is not specified in the parameter file, any instance of that parameter in the simulation will assume its default value specified in the master_params.yaml file.

Some important global parameters to consider are:

* seed: simulation seed to use with random number generator 
* run_name: prefix for all output files
* n_runs: number of individual runs of each parameter type
* n_random: number of samples from a random parameter space (see more below)
* n_dim: number of dimensions of simulation
* n_periodic: number of periodic dimensions of simulation
* delta: simulation time step
* n_steps: total number of steps in each simulation
* system_radius: "box radius" of system
* graph_flag: run with graphics enabled
* n_graph: how many simulation steps to take between updating graphics
* movie_flag: whether to record the graphics frames into bitmaps
* movie_directory: local directory used to save the recorded bitmaps
* cell_length: linear size of cells used by cell lists (used by interactions)
* n_update_cells: how often to update the cell lists
* f_cutoff: maximum force to apply between pair interactions
* thermo_flag: whether to output thermodynamics outputs (stress tensors, etc)
* n_thermo: how often to output the thermodynamics outputs
* wca_eps/wca_sig/ss_a/ss_rs/ss_eps: interaction potential parameters

Some important species parameters to consider are:

* num: how many to insert into system
* insertion_type: how to insert object into system (e.g. random)
* overlap: whether species can overlap (should be 0 with interactions)
* draw_type: (orientation or flat) whether to color by orientation or not
* color: a double that specifies the flat rgb value of the object
* posit_flag: whether to output position files
* n_posit: how often to output position files
* spec_flag: whether to output species files
* n_spec: how often to output species files
* checkpoint_flag: whether to output checkpoint files
* n_checkpoint: how often to output checkpoint files

All parameters used in the simulation, along with their default values and data types, are specified in the master_params.yaml file in the src folder.

## Adding new parameters

simcore comes with it's own parameter initialization tool, simcore_config, which can be installed by following the above installation instructions for simcore and then doing 'make simcore_config'.  simcore_config makes it easy to add new parameters to the simulation without mucking around in the source code. Just add your new parameter to the master_params.yaml file using the following format: 

```
new_parameter_name: [default_parameter_value, parameter_type] 
```
 
Running simcore_config will look at all the parameters in the master_params.yaml file and add them seamlessly to the proper simcore files, and you can begin using them in your classes right away after recompiling simcore. NOTE: At the time of this writing (12/3/2017), yaml-cpp is inconsistent about its treatment of boolean values. For type-safety reasons, it's best to use integers instead of bools when adding flag parameters.  

### Parameter sets

Using parameter sets, it becomes easier to run many simulations over a given parameter space. There are two types of parameter sets possible with simcore: defined and random. Each parameter set type works the same with both global parameters and species parameters.

#### Defined parameter sets
  
Defined parameter sets are used in this way in the parameter file:

```
seed: 4916819461895
run_name: defined_set
n_runs: N
parameter_name1: param_value1
parameter_name2: [param_value2, param_value3]
parameter_name3: [param_value4, param_value5]
```

Parameters specified in this way (as lists of parameters) will be iterated over until every possible combination of parameters has been run. In this example, simcore will run N simulations each of the following 4 parameter sets:

```
seed: random_seed_1
run_name: defined_set_v000
n_runs: N
parameter_name1: param_value1
parameter_name2: param_value2
parameter_name3: param_value4

seed: random_seed_2
run_name: defined_set_v001
n_runs: N
parameter_name1: param_value1
parameter_name2: param_value2
parameter_name3: param_value5

seed: random_seed_3
run_name: defined_set_v002
n_runs: N
parameter_name1: param_value1
parameter_name2: param_value3
parameter_name3: param_value4

seed: random_seed_4
run_name: defined_set_v003
n_runs: N
parameter_name1: param_value1
parameter_name2: param_value3
parameter_name3: param_value5
```

#### Random parameter sets

Random parameter sets are designed specifically to be used with polynomial-chaos theory for n-dimensional parameter spaces for large n. Random sets are used in the following way:

```
seed: 2546954828254
n_runs: N
n_random: M
parameter_name1: param_value1
parameter_name2: [R, A, B] # sets to random real in range (A,B)
parameter_name3: [RINT, C, D] # sets to random int in range [C,D]
parameter_name4: [RLOG, F, G] # sets to 10^K for rand real K in range (F,G)
```

Given this parameter file, simcore will run N simulations each of M random parameter sets. The random parameter sets are generated in ranges specified in the lists that are prefixed by the R, RINT, RLOG options.

In this example, the sampled parameter space has dimensionality of n=3, since there are only three parameters we are sampling over. Each parameter set will have a random real number for parameter_name2 in the the range (A,B), a random integer in the range [C,D] for parameter_name3, and will set parameter_name4 to 10^K for random real number K in the range (F,G).  simcore will then run each parameter set N times each with a unique seed, and repeat this random process M times. It will therefore take N samples of M random points in the n-dimensional parameter space.  

## Interactions
  
The InteractionEngine in simcore was written with short-range interactions in mind. For this reason, interactions are treated by considering pair-wise interactions between neighboring interactor-elements that make up a composite object (e.g. small, rigid segments that compose a flexible filament). For this reason, interactions use cell lists to improve performance. Furthermore, simulating large objects in simcore requires representing the object as a composite of smaller, simple objects (thus, SIMple Composite Object REpresentation). An example of how a large object should be decomposed into simple objects is done in the Filament class.

### Potentials
  
simcore is designed to be able to use interchangable potentials for various objects. However, potentials need to be added manually as a subclass of PotentialBase, included in PotentialManager, and a corresponding potential_type added to definitions.h for lookup purposes (see the InitPotentials method in PotentialManager.h for examples).

## Outputs
  
simcore has four output types. Three are species specific (posit, spec, checkpoint), and the fourth is the statistical information file (thermo). All files are written in binary.

The posit file has the following header format:

```
int n_steps, int n_posit, double delta 
```

Followed by n_steps/n_posit lines of data with the format:

```
double position[3]
double scaled_position[3]
double orientation[3]
double diameter
double length
```

Where the scaled position is position mapped into the periodic coordinate space. The position itself gives the particle trajectory over time independent of periodicity.  

The spec file is a custom output file for each species, and can have the same information as the posit file or additional information if needed.

The checkpoint file is almost a copy of the spec file, except it also contains the random number generator information and is overwritten every n_checkpoint steps in the simulation. It can therefore be used to resume a simulation that ended prematurely.

The thermo file contains the following header information:

```
int n_steps, int n_thermo, double delta, int n_dim
```

followed by n_steps/n_thermo lines of data in the following format:

```
double unit_cell[9]
double pressure_tensor[9]
double pressure
double volume
```

Where the pressure is the isometric pressure, and the pressure tensor is calculated from the time-averaged stress tensor.

## Data analysis
  
If analysis operations of output files are already defined for your species, as is the case for the Filament species, analyzing outputs is a simple matter. First, make sure the desired analysis flag is set in the species parameters for that species.

For example, in the Filament species there is a persistence length analysis that produces .mse2e files that tracks the mean-square end-to-end distance of semiflexible filaments. This is triggered by a parameter lp_analysis=1, which can be set in the parameter file.

Anaylses are run by running simcore in the following way:
  
```
./simcore -a parameter_file.yaml.
```
  
NOTE: It is important to keep in mind that the parameter_file should be identical to the parameter file used to generate the outputs. There are a few exceptions that only affect post-processing, such as analysis flags, but this is true in general.

The way inputs and outputs are meant to work in simcore is such that during a simulation, output data are generated in the posit, spec, and checkpoint formats, and during analysis, the same output data are read back into the data structures in simcore for processing. The .posit files just contain bare-bones information that allow many types of simple analyses, but .spec files should in general contain all the necessary information to recreate the trajectory for a member of a species. 

For a new species analysis method, the analysis routines should be defined in the species container class, rather than the species member class, and called by the inherited RunAnalysis method of the SpeciesBase class (and likewise for analysis initialization and finalization, see examples for details).

For example, the RunSpiralAnalysis routine is called by the RunAnalysis method in FilamentSpecies, which uses the Filament .spec file as an input to do the necessary analysis, whose results are placed into a new file ending in filament.spiral. See Filament and FilamentSpecies for examples of how analyses can be initialized, processed, etc.

## Disclaimer

simcore was written for my personal academic use and in its current state is not intended to be used by the general public. If you are insane (and somehow also patient) and would like to run simcore for whatever reason, you can contact me for help and (if I have time) I will do what I can to offer assistance. In addition, the README provided here is in no way a complete documentation of the software. The simcore software is covered by the MIT license.



