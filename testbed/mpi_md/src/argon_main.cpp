#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <string>
#include <cstdlib>
#include <cmath>
#include <memory>

// Boost yuckiness
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/back_inserter.hpp>

#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

#if defined(_OPENMP)
#define ENABLE_OPENMP
#endif

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

#include<mpi.h>

#define BLEN 200

const double kboltz = 0.0019872067;
const double mvsq2e = 2390.05736153349;

const double cellrat = 2.0;
const int cellfreq = 2;

// Cell list struct
struct _cell {
    int natoms;
    int owner;
    int *idxlist;
};
typedef struct _cell cell_t;

class SystemProperties {
  private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int /*version*/) {
        ar & natoms;
        ar & nfi;
        ar & nsteps;
        ar & nthreads;

        ar & dt;
        ar & mass;
        ar & eps;
        ar & sigma;
        ar & box;
        ar & rcut;
        ar & ekin;
        ar & epot;
        ar & temp;
    }
  public:
    SystemProperties() : natoms(0), nfi(0), nsteps(0), nthreads(1), dt(0.0), mass(0.0), eps(0.0), sigma(0.0), box(0.0), rcut(0.0),
                         ekin(0.0), epot(0.0), temp(0.0) {}

    void print() {
        int mpirank;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
        printf("\t[%d]sysprop(%p) ->\n", mpirank, this);
        printf("\t[%d]\t{natoms:%d}, {nfi:%d}, {nsteps:%d}, {dt:%2.2f}, {mass:%2.2f}\n",
                mpirank, natoms, nfi, nsteps, dt, mass);
        printf("\t[%d]\t{eps:%2.2f}, {sigma:%2.2f}, {box:%2.2f}, {rcut:%2.2f}\n",
                mpirank, eps, sigma, box, rcut);
    }

    int natoms, nfi, nsteps, nthreads;
    double dt, mass, eps, sigma, box, rcut;
    double ekin, epot, temp;
};

class MDSys {
  private:

  public:
    MDSys() {}
    ~MDSys() {
        delete(sys_prop);
    }

    void print() {
        printf("[%d]sys(%p) ->\n", mpirank, this);
        sys_prop->print();
    }

    int ngrid, ncell, npair, nidx;
    int* plist;

    double delta;
    double *pos, *vel, *frc, *buf;

    int mpirank, mpisize;
    MPI_Comm mpicomm;

    SystemProperties* sys_prop;
    cell_t* clist;
};

__attribute__((always_inline))
static void azzero(double *d, const int n) {
    for (int i = 0; i < n; ++i) {
        d[i] = 0.0;
    }
}

__attribute__((always_inline,pure))
static double pbc(double x, const double boxby2, const double box) {
    while (x >  boxby2) x -= box;
    while (x < -boxby2) x += box;
    return x;
}

// Hijac the hoomd version to serialize the damn code?
template<typename T>
void bcast(T& val, unsigned int root, const MPI_Comm mpi_comm) {
    int rank;
    MPI_Comm_rank(mpi_comm, &rank);

    char *buf = nullptr;
    int recv_count;
    if (rank == (int)root) {
        std::string str;
        boost::iostreams::back_insert_device <std::string> inserter(str);
        boost::iostreams::stream<boost::iostreams::back_insert_device<std::string>> s(inserter);
        boost::archive::binary_oarchive ar(s);

        // serialize object
        ar << val;

        // flush apparently
        s.flush();

        // copy and send
        recv_count = str.size();
        buf = new char[recv_count];
        str.copy(buf, recv_count);
    }

    // Broadcast that shit
    MPI_Bcast(&recv_count, 1, MPI_INT, root, mpi_comm);
    if (rank != (int)root) {
        buf = new char[recv_count];
    }
    MPI_Bcast(buf, recv_count, MPI_BYTE, root, mpi_comm);

    if (rank != (int)root) {
        // deserialize that shit
        std::string str(buf, recv_count);
        boost::iostreams::basic_array_source<char> dev(str.data(), str.size());
        boost::iostreams::stream<boost::iostreams::basic_array_source<char>> s(dev);
        boost::archive::binary_iarchive ar(s);

        ar >> val;
    }

    delete[] buf;
}

// Update the freaking cells
static void updatecells(MDSys* sys) {
    int ngrid, ncell, npair, midx, natoms;
    double delta, boxby2, boxoffs;
    boxby2 = 0.5 * sys->sys_prop->box;
    natoms = sys->sys_prop->natoms;

    if (sys->clist == nullptr) {
        int nidx;

        ngrid = floor(cellrat * sys->sys_prop->box / sys->sys_prop->rcut);
        ncell = ngrid*ngrid*ngrid;
        delta = sys->sys_prop->box / ngrid;
        boxoffs = boxby2 - 0.5*delta;

        sys->delta = delta;
        sys->ngrid = ngrid;
        sys->ncell = ncell;

        // allocate cell list storage
        sys->clist = (cell_t *) malloc(ncell*sizeof(cell_t));
        sys->plist = (int *) malloc(2*ncell*ncell*sizeof(int));

        // allocate index lists within cell.  cell density < 2x avg. density
        nidx = 2*natoms / ncell + 2;
        nidx = ((nidx/2) + 1) * 2;
        sys->nidx = nidx;
        for (int i=0; i<ncell; ++i) {
            sys->clist[i].idxlist = (int *) malloc(nidx*sizeof(int));
        }

        // bild the cell list
        npair = 0;
        for (int i=0; i < ncell-1; ++i) {
            int k;
            double x1,y1,z1;
            
            k  = i/ngrid/ngrid;
            x1 = k*delta - boxoffs;
            y1 = ((i-(k*ngrid*ngrid))/ngrid)*delta - boxoffs;
            z1 = (i % ngrid)*delta - boxoffs;

            for (int j=i+1; j<ncell; ++j) {
                double x2,y2,z2,rx,ry,rz;
                
                k  = j/ngrid/ngrid;
                x2 = k*delta - boxoffs;
                y2 = ((j-(k*ngrid*ngrid))/ngrid)*delta - boxoffs;
                z2 = (j % ngrid)*delta - boxoffs;

                rx=pbc(x1 - x2, boxby2, sys->sys_prop->box);
                ry=pbc(y1 - y2, boxby2, sys->sys_prop->box);
                rz=pbc(z1 - z2, boxby2, sys->sys_prop->box);

                // check for cells on a line
                if (fabs(rx) > sys->sys_prop->rcut+delta) continue;
                if (fabs(ry) > sys->sys_prop->rcut+delta) continue;
                if (fabs(rz) > sys->sys_prop->rcut+delta) continue;

                // Check or cells in a plane
                if (sqrt(rx*rx+ry*ry) > (sys->sys_prop->rcut+sqrt(2.0)*delta)) continue;
                if (sqrt(rx*rx+rz*rz) > (sys->sys_prop->rcut+sqrt(2.0)*delta)) continue;
                if (sqrt(ry*ry+rz*rz) > (sys->sys_prop->rcut+sqrt(2.0)*delta)) continue;

                // All other cell checks
                if (sqrt(rx*rx + ry*ry + rz*rz) > (sqrt(3.0)*delta+sys->sys_prop->rcut)) continue;
               
                // Cells are close enough, add
                sys->plist[2*npair  ] = i;
                sys->plist[2*npair+1] = j;
                ++npair;
            }
        }
        sys->npair = npair;
        if (sys->mpirank == 0) {
            printf("[%d]Cell list has %dx%dx%d=%d cells with %d/%d pairs and "
                   "%d atoms/celllist.\n", sys->mpirank, ngrid, ngrid, ngrid, sys->ncell, 
                   sys->npair, ncell*(ncell-1)/2, nidx);
        }
    }

    // reset cell list
    ncell=sys->ncell;
    delta=sys->delta;
    ngrid=sys->ngrid;
    
    for (int i=0; i < sys->ncell; ++i) {
        sys->clist[i].natoms=0;
    }

    boxoffs= boxby2 - 0.5*delta;
    midx=0;
    for (int i=0; i < natoms; ++i) {
        int idx,j,k,m,n;
        
        k=floor((pbc(sys->pos[i], boxby2, sys->sys_prop->box)+boxby2)/delta);
        m=floor((pbc(sys->pos[natoms + i], boxby2, sys->sys_prop->box)+boxby2)/delta);
        n=floor((pbc(sys->pos[2*natoms + i], boxby2, sys->sys_prop->box)+boxby2)/delta);
        j = ngrid*ngrid*k+ngrid*m+n;

        idx = sys->clist[j].natoms;
        sys->clist[j].idxlist[idx]=i;
        ++idx;
        sys->clist[j].natoms = idx;
        if (idx > midx) midx=idx;
    }
    if (midx > sys->nidx) {
        printf("[%d]overflow in cell list: %d/%d atoms/cells.\n", sys->mpirank, midx, sys->nidx);
        MPI_Abort(MPI_COMM_WORLD, 1);
        exit(1);
    }
    return;
}

static void free_cell_list(MDSys* sys) {
    if (sys->clist == nullptr)
        return;

    for (int i = 0; i < sys->ncell; ++i) {
        delete(sys->clist[i].idxlist);
    }

    delete(sys->clist);
    sys->clist = nullptr;
    sys->ncell = 0;
}

static void ekin(MDSys* sys) {
    sys->sys_prop->ekin = 0.0;
    for (int i = 0; i < 3*sys->sys_prop->natoms; ++i) {
        sys->sys_prop->ekin += sys->vel[i]*sys->vel[i];
    }
    sys->sys_prop->ekin *= 0.5*mvsq2e*sys->sys_prop->mass;
    sys->sys_prop->temp = 2.0*sys->sys_prop->ekin/(3.0*sys->sys_prop->natoms - 3.0)/kboltz;
}

static void force(MDSys* sys) {
#ifdef DEBUG
    printf("[%d of %d] Starting force calculation\n", sys->mpirank, sys->mpisize);
#endif
    double epot = 0.0;
    int natoms = sys->sys_prop->natoms;

    // Broadcast positions
    MPI_Bcast(sys->pos, 3*natoms, MPI_DOUBLE, 0, sys->mpicomm);

    // update the cell list if needed
    if ((sys->clist == nullptr) || (sys->sys_prop->nfi % cellfreq) == 0) {
        updatecells(sys);
    }

#ifdef ENABLE_OPENMP
#pragma omp parallel reduction(+:epot)
#endif
    {
        double c12, c6, boxby2, rcsq;
        double *fx, *fy, *fz;
        const double *rx, *ry, *rz;
        int tid;


        // Precompute constants
        c12 = 4.0*sys->sys_prop->eps*pow(sys->sys_prop->sigma, 12.0);
        c6  = 4.0*sys->sys_prop->eps*pow(sys->sys_prop->sigma,  6.0);
        rcsq = sys->sys_prop->rcut * sys->sys_prop->rcut;
        boxby2 = 0.5*sys->sys_prop->box;
        epot = 0.0;

        // get my force array
#ifdef ENABLE_OPENMP
        tid = omp_get_thread_num();
#else
        tid = 0;
#endif
#ifdef DEBUG
        printf("[\t{%d of %d} Getting force arrays\n", tid, sys->sys_prop->nthreads);
#endif
        fx = sys->buf + (3*tid*natoms);
        azzero(fx,3*natoms);
        fy = sys->buf + ((3*tid+1)*natoms);
        fz = sys->buf + ((3*tid+2)*natoms);
        rx = sys->pos;
        ry = sys->pos + natoms;
        rz = sys->pos + 2*natoms;

        // Our incrementer, based on mpirank and threads
        int incr = sys->mpisize * sys->sys_prop->nthreads;

        // Interactions within our cell
        for (int n = 0; n < sys->ncell; n += incr) {
            const cell_t *c1;

            int i = n + sys->mpirank*sys->sys_prop->nthreads + tid;
            if (i >= sys->ncell) break;
            c1 = sys->clist + i;

            for (int j = 0; j < c1->natoms - 1; ++j) {

                int ii = c1->idxlist[j];
                double rx1 = rx[ii];
                double ry1 = ry[ii];
                double rz1 = rz[ii];

                for (int k = j + 1; k < c1->natoms; ++k) {
                    int jj = c1->idxlist[k];

                    // Calc distance
                    double rx2 = pbc(rx1 - rx[jj], boxby2, sys->sys_prop->box);
                    double ry2 = pbc(ry1 - ry[jj], boxby2, sys->sys_prop->box);
                    double rz2 = pbc(rz1 - rz[jj], boxby2, sys->sys_prop->box);
                    double rsq = rx2*rx2 + ry2*ry2 + rz2*rz2;

                    // If within cutoff, do this
                    if (rsq < rcsq) {
                        double r6,rinv,ffac;

                        rinv = 1.0/rsq;
                        r6 = rinv*rinv*rinv;

                        ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
                        epot += r6*(c12*r6 - c6);

                        fx[ii] += rx2*ffac;
                        fy[ii] += ry2*ffac;
                        fz[ii] += rz2*ffac;
                        fx[jj] -= rx2*ffac;
                        fy[jj] -= ry2*ffac;
                        fz[jj] -= rz2*ffac;
                    }
                } // particle 2
            } // particle 1
        } // end of interactions in this cell

        // Interactions with neighbor cells
        for (int n = 0; n < sys->npair; n += incr) {

            int i = n + sys->mpirank*sys->sys_prop->nthreads + tid;
            if (i >= sys->npair) break;
            cell_t* c1 = sys->clist + sys->plist[2*i];
            cell_t* c2 = sys->clist + sys->plist[2*i+1];

            for (int j = 0; j < c1->natoms; ++j) {

                int ii = c1->idxlist[j];
                double rx1 = rx[ii];
                double ry1 = ry[ii];
                double rz1 = rz[ii];

                for (int k = 0; k < c2->natoms; ++k) {
                
                    int jj = c2->idxlist[k];
                    double rx2 = pbc(rx1 - rx[jj], boxby2, sys->sys_prop->box);
                    double ry2 = pbc(ry1 - ry[jj], boxby2, sys->sys_prop->box);
                    double rz2 = pbc(rz1 - rz[jj], boxby2, sys->sys_prop->box);
                    double rsq = rx2*rx2 + ry2*ry2 + rz2*rz2;

                    // If within cutoff, do this
                    if (rsq < rcsq) {
                        double r6,rinv,ffac;

                        rinv = 1.0/rsq;
                        r6 = rinv*rinv*rinv;

                        ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
                        epot += r6*(c12*r6 - c6);

                        fx[ii] += rx2*ffac;
                        fy[ii] += ry2*ffac;
                        fz[ii] += rz2*ffac;
                        fx[jj] -= rx2*ffac;
                        fy[jj] -= ry2*ffac;
                        fz[jj] -= rz2*ffac;
                    }
                } // atom2
            } // atom1
        } // Done with neighbor interactions

        // Before reducing hte forces, do a local reduce in openmp
#ifdef ENABLE_OPENMP
#pragma omp barrier
#endif
        int n = 1 + (3*natoms / sys->sys_prop->nthreads);
        int fromidx = tid * n;
        int toidx = fromidx + n;
        if (toidx > 3*natoms) toidx = 3*natoms;

        // Reduce ze forces
        for (int nn = 1; nn < sys->sys_prop->nthreads; ++nn) {
            int offs, j;

            offs = 3*nn*natoms;
            
            for (j = fromidx; j < toidx; ++j) {
                sys->buf[j] += sys->buf[offs+j];
            }
        }
    } // omp parallel reduction (+:epot)

    // Reduce
    MPI_Reduce(sys->buf, sys->frc, 3*natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
    MPI_Reduce(&epot, &sys->sys_prop->epot, 1, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
}

//vel verlet
static void velverlet(MDSys* sys) {
    double dtmf = 0.5*sys->sys_prop->dt / mvsq2e / sys->sys_prop->mass;

    for (int i = 0; i < 3*sys->sys_prop->natoms; ++i) {
        sys->vel[i] += dtmf * sys->frc[i];
        sys->pos[i] += sys->sys_prop->dt*sys->vel[i];
    }

    force(sys);

    for (int i = 0; i < 3*sys->sys_prop->natoms; ++i) {
        sys->vel[i] += dtmf * sys->frc[i];
    }
}

// output
static void output(MDSys *sys, FILE *erg, FILE *traj)
{
    int i,natoms;
    natoms=sys->sys_prop->natoms;
    
    printf("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->sys_prop->nfi, sys->sys_prop->temp, sys->sys_prop->ekin, sys->sys_prop->epot, 
            sys->sys_prop->ekin+sys->sys_prop->epot);
    fprintf(erg,"% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->sys_prop->nfi, sys->sys_prop->temp, sys->sys_prop->ekin, sys->sys_prop->epot, 
            sys->sys_prop->ekin+sys->sys_prop->epot);
    fprintf(traj,"%d\n nfi=%d etot=%20.8f\n", sys->sys_prop->natoms, sys->sys_prop->nfi, sys->sys_prop->ekin+sys->sys_prop->epot);
    for (i=0; i<natoms; ++i) {
        fprintf(traj, "Ar  %20.8f %20.8f %20.8f\n", sys->pos[i], sys->pos[natoms+i], sys->pos[2*natoms+i]);
    }
}

int main(int argc, char** argv) {

    MDSys sys;
    MPI_Comm mpicomm = MPI_COMM_WORLD;
    int mpirank, mpisize;

    // Variables that will get assigned to structure
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN];
    int natoms, nsteps, nprint;
    double mass, eps, sigma, rcut, box, dt;

    FILE *fp, *erg, *traj;

    // Init MPI and get the rank, size
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(mpicomm, &mpirank);
    MPI_Comm_size(mpicomm, &mpisize);
    sys.mpirank = mpirank;
    sys.mpisize = mpisize;
    sys.mpicomm = mpicomm;
    sys.sys_prop = new SystemProperties();

#ifdef ENABLE_OPENMP
#pragma omp parallel
    {
        if (omp_get_thread_num() == 0) {
            sys.sys_prop->nthreads = omp_get_num_threads();
        }
    }
#else
    sys.sys_prop->nthreads = 1;
#endif


    if (mpirank == 0) {
        printf("Running in hybrid parallel with %d MPI tasks using %d threads each\n", sys.mpisize, sys.sys_prop->nthreads);
        std::ifstream input_file(argv[1]);
        std::string line;
        if (std::getline(input_file, line)) {
            std::stringstream linestream(line);
            linestream >> natoms;
            sys.sys_prop->natoms = natoms;
        }
        if(std::getline(input_file, line)) {
            std::stringstream linestream(line);
            linestream >> mass;
            sys.sys_prop->mass = mass;
        }
        if(std::getline(input_file, line)) {
            std::stringstream linestream(line);
            linestream >> eps;
            sys.sys_prop->eps = eps;
        }
        if(std::getline(input_file, line)) {
            std::stringstream linestream(line);
            linestream >> sigma; 
            sys.sys_prop->sigma = sigma;
        }
        if(std::getline(input_file, line)) {
            std::stringstream linestream(line);
            linestream >> rcut;
            sys.sys_prop->rcut = rcut;
        }
        if(std::getline(input_file, line)) {
            std::stringstream linestream(line);
            linestream >> box;
            sys.sys_prop->box = box;
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
            sys.sys_prop->nsteps = nsteps;
        }
        if(std::getline(input_file, line)) {
            std::stringstream linestream(line);
            linestream >> dt;
            sys.sys_prop->dt = dt;
        }
        if(std::getline(input_file, line)) {
            std::stringstream linestream(line);
            linestream >> nprint;
        }
    }
    // Broadcast, oh joy of joys
    bcast(sys.sys_prop, 0, sys.mpicomm);
    //sys.print();

    // Now we're getting dangerous, since we're initializing mallocs
    sys.pos = (double *)malloc(3*sys.sys_prop->natoms*sizeof(double));
    sys.vel = (double *)malloc(3*sys.sys_prop->natoms*sizeof(double));
    sys.frc = (double *)malloc(3*sys.sys_prop->natoms*sizeof(double));

    sys.buf = (double *)malloc(3*sys.sys_prop->natoms*sys.sys_prop->nthreads*sizeof(double));

    if (sys.mpirank == 0) {
        fp = fopen(restfile,"r");
        if (fp) {
            natoms = sys.sys_prop->natoms;
            for (int i = 0; i < natoms; ++i) {
                fscanf(fp, "%lf%lf%lf", sys.pos+i, sys.pos+natoms+i, sys.pos+2*natoms+i);
            }
            for (int i = 0; i < natoms; ++i) {
                fscanf(fp, "%lf%lf%lf", sys.vel+i, sys.vel+natoms+i, sys.vel+2*natoms+i);
            }
            fclose(fp);
            azzero(sys.frc, 3*sys.sys_prop->natoms);
        } else {
            perror("cannot read restart file");
            MPI_Abort(sys.mpicomm, 3);
            return 3;
        }
    } else {
        azzero(sys.vel, 3*sys.sys_prop->natoms);
        azzero(sys.frc, 3*sys.sys_prop->natoms);
    }

    // Setup cell list
    sys.clist = nullptr;
    sys.plist = nullptr;

    sys.sys_prop->nfi = 0;
    force(&sys);
    ekin(&sys);

    if (sys.mpirank == 0) {
        erg = fopen(ergfile,"w");
        traj=fopen(trajfile,"w");

        printf("Starting simulation with %d atoms for %d steps.\n",sys.sys_prop->natoms, sys.sys_prop->nsteps);
        printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
        output(&sys, erg, traj);
    }

    // Loopy loop
    auto starttime = MPI_Wtime();
    for (sys.sys_prop->nfi = 1; sys.sys_prop->nfi <= sys.sys_prop->nsteps; ++sys.sys_prop->nfi) {
        if (sys.mpirank == 0) {
            if ((sys.sys_prop->nfi % nprint) == 0) {
                output(&sys, erg, traj);
            }
        }

        velverlet(&sys);
        ekin(&sys);
    }
    auto endtime = MPI_Wtime();

    // Clean up after ourselves
    if (sys.mpirank == 0) {
        printf("Simulation Done. Loop time: %.2f seconds\n", endtime-starttime);
        fclose(erg);
        fclose(traj);
    }

    delete(sys.pos);
    delete(sys.vel);
    delete(sys.frc);
    delete(sys.buf);
    free_cell_list(&sys);

    MPI_Barrier(sys.mpicomm);
    MPI_Finalize();
    return 0;
}
