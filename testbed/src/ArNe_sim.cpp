// Main Argon Neon simulation (2species)

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "Ar.h"
#include "ArSpecies.h"
#include "system.h"
#include "cell_list.h"

#include "lennard_jones_12_6.h"

#define BLEN 200

const int cellfreq = 4;

int main(int argc, char **argv) {
    // Get the input via command arguments
    if (argc != 2)
        return 1;
    std::ifstream input_file(argv[1]);
    std::string line;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN];
    FILE *fp, *traj, *erg;
    unsigned int natoms = 0, nsteps = 0;
    int nprint = 0;
    double mass, eps, sigma, rcut, box, dt;

    // Read in the test files
    if(std::getline(input_file, line)) {
        std::stringstream linestream(line);
        linestream >> natoms;
    }
    if(std::getline(input_file, line)) {
        std::stringstream linestream(line);
        linestream >> mass;
    }
    if(std::getline(input_file, line)) {
        std::stringstream linestream(line);
        linestream >> eps; 
    }
    if(std::getline(input_file, line)) {
        std::stringstream linestream(line);
        linestream >> sigma; 
    }
    if(std::getline(input_file, line)) {
        std::stringstream linestream(line);
        linestream >> rcut;
    }
    if(std::getline(input_file, line)) {
        std::stringstream linestream(line);
        linestream >> box;
    }
    if(std::getline(input_file, line)) {
        std::stringstream linestream(line);
        linestream >> restfile;
    }
    if(std::getline(input_file, line)) {
        std::stringstream linestream(line);
        linestream >> trajfile;
    }
    if(std::getline(input_file, line)) {
        std::stringstream linestream(line);
        linestream >> ergfile;
    }
    if(std::getline(input_file, line)) {
        std::stringstream linestream(line);
        linestream >> nsteps;
    }
    if(std::getline(input_file, line)) {
        std::stringstream linestream(line);
        linestream >> dt;
    }
    if(std::getline(input_file, line)) {
        std::stringstream linestream(line);
        linestream >> nprint;
    }

    // Creat the ArNe system
    std::cout << "Creating Argon/Neon System\n";
    auto ArNeSys = new SystemArch(natoms, box, dt);

    std::cout << "Setting Argon parameters\n";
    ArSpecies* ArSpec = new ArSpecies();
    ArSpec->setSpecies(mass, eps, sigma, rcut, 0, "Ar");
    ArNeSys->addSpecies(0, ArSpec);
    auto lj126_ar = potentialFactory<LJ126>(eps, sigma, rcut, box);

    /*std::cout << "Setting Neon parameters\n";
    ArSpecies* NeSpec = new ArSpecies();
    // Hard code in the neon parameters, because
    mass = 20.1797;
    eps = 0.069;
    sigma = 2.782;
    rcut = 2.5*sigma;
    NeSpec->setSpecies(mass, eps, sigma, rcut, 1, "Ne");
    ArNeSys->addSpecies(1, NeSpec);
    auto lj126_ne = potentialFactory<LJ126>(eps, sigma, rcut, box);*/

    // Read in the particles
    // Interplex Ar and Ne
    fp = fopen(restfile, "r");
    double rx, ry, rz, vx, vy, vz;
    if (fp) {
        for (int i = 0; i < (int)natoms; ++i) {
            fscanf(fp, "%lf%lf%lf", &rx, &ry, &rz);
            int idx;
            //if (i%2 == 0)
                idx = ArNeSys->addParticle(0);
            //else
            //    idx = ArNeSys->addParticle(1);
            assert(i == idx);
            auto arp = ArNeSys->getParticle(idx);
            arp->setXYZ(rx, ry, rz);
        }
        for (int i = 0; i < (int)natoms; ++i) {
            fscanf(fp, "%lf%lf%lf", &vx, &vy, &vz);
            auto arp = ArNeSys->getParticle(i);
            arp->setV(vx, vy, vz);
        }
        fclose(fp);
    } else {
        std::cerr << "Cannot read restart file\n";
    }

    // Initialize the engine
    ArNeSys->generateCellList();

    // Dump everything to check
    ArNeSys->dump();
    ArNeSys->checkConsistency();

    // Create the potentials
    // The order matters for Newton's 3rd law
    ArNeSys->addPotential(0, 0, lj126_ar);
    //ArNeSys->addPotential(1, 1, lj126_ne);
    // Write information
    ArNeSys->dumpPotentials();

    // Finish initialization
    ArNeSys->forceMP();
    ArNeSys->ukin();

    erg = fopen(ergfile, "w");
    traj = fopen(trajfile, "w");

    printf("Starting simulation with %d atoms for %d steps.\n", ArNeSys->nParticles(), nsteps);
    printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
    ArNeSys->output(erg, traj, 0);

    // main simulation loop
    for (unsigned int i_step = 1; i_step <= nsteps; ++i_step) {
        // output
        if ((i_step % nprint) == 0)
            ArNeSys->output(erg, traj, i_step);

        // do position update
        ArNeSys->velverlet();
        ArNeSys->ukin();

        // update the cell list
        if ((i_step % cellfreq) == 0)
            ArNeSys->updateCellList();
    }

    return 0;
}
