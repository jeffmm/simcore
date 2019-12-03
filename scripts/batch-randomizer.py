"""Module writes num_cores summit batch files with approx num_jobs/num_cores jobs to run
on each core. Assumes parameter files exist with format <run_name>_v<run_number>.yaml.
"""

import random
import sys
import numpy as np


def writeSummitJob(run_name, run_num, time=4):
    """ Function inputs: run_name (str), run_num (int), time in hours (int)
    Function behavior: writes slurm batch file for simcore job on summit with
    title derived from run_name + run_num
    """
    assert isinstance(
        run_name, str
    ), "writeSummitJob did not receive a run_name as string type"

    assert isinstance(time, int), "writeSummitJob did not receive a time as int type"
    assert isinstance(
        run_num, int
    ), "writeSummitJob did not receive a run_num as int type"
    job_name = run_name + "_v" + str(run_num).zfill(3)
    f = open("sjob_" + job_name + ".sh", "w")
    f.write("#!/bin/bash\n")
    f.write("#SBATCH --job-name=" + job_name + "\n")
    f.write("#SBATCH --time " + str(time).zfill(2) + ":00:00\n")
    f.write("#SBATCH --nodes 1\n")
    f.write("#SBATCH --ntasks 1\n")
    f.write("#SBATCH --ntasks-per-node 1\n")
    f.write("#SBATCH --cpus-per-task 1\n")
    f.write("#SBATCH --output " + job_name + ".out\n")
    f.write("#SBATCH --error " + job_name + ".err\n")
    f.write("#SBATCH --account ucb-summit-smr\n")
    f.write("#SBATCH --qos=condo\n")
    f.write("#SBATCH --partition=shas\n")
    f.write("#SBATCH --mail-type=ALL\n")
    f.write("#SBATCH --mail-user=jemo9179@colorado.edu\n\n")

    f.write("module purge\n")
    f.write("module load singularity/3.3.0\n\n")
    f.write(
        'echo "Executing simcore job ${SLURM_JOB_ID} on'
        ' ${SLURM_CPUS_PER_TASK} threads"\n'
    )
    f.write("export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}\n\n")
    f.write(
        "singularity run simcore_latest.sif ./simcore.exe "
        "" + job_name + "_params.yaml\n"
    )
    f.close()


def writeMultiJob(
    run_name, num_list, num_cores, job_num=0, time=24, email="none", n_run=-1
):
    """ Function inputs: run_name (str), num_list (list of ints), core_num
    (int), time in hours (int)

        Function behavior: writes slurm batch file for many simcore jobs on
        summit based on numbers in num_list on summit with jobs derived from
        run_name + number
    """
    assert isinstance(
        run_name, str
    ), "writeMultiJob did not receive a run_name as string type"
    assert isinstance(time, int), "writeMultiJob did not receive a time as int type"
    assert isinstance(
        job_num, int
    ), "writeMultiJob did not receive a core_num as int type"
    assert isinstance(
        num_cores, int
    ), "writeMultiJob did not receive a num_cores as int type"
    assert isinstance(
        num_list, list
    ), "writeMultiJob did not receive a num_list as list type"
    for num in num_list:
        assert isinstance(
            num, int
        ), "writeMultiJob did not receive a list of ints for num_list"
    job_name_core = run_name + "_j" + str(job_num).zfill(2)
    f = open("sjob_" + job_name_core + ".sh", "w")
    f.write("#!/bin/bash\n")
    f.write("#SBATCH --job-name=" + job_name_core + "\n")
    f.write("#SBATCH --time " + str(time).zfill(2) + ":00:00\n")
    f.write("#SBATCH --nodes 1\n")
    f.write("#SBATCH --ntasks 1\n")
    f.write("#SBATCH --ntasks-per-node 1\n")
    f.write("#SBATCH --cpus-per-task 1\n")
    f.write("#SBATCH --output " + job_name_core + ".out\n")
    f.write("#SBATCH --error " + job_name_core + ".err\n")
    f.write("#SBATCH --account ucb-summit-smr\n")
    f.write("#SBATCH --qos=condo\n")
    f.write("#SBATCH --partition=shas\n")
    if email != "none":
        f.write("#SBATCH --mail-type=ALL\n")
        f.write("#SBATCH --mail-user=" + email + "\n\n")

    f.write("module purge\n")
    f.write("module load singularity/3.3.0\n\n")
    f.write(
        'echo "Executing simcore job ${SLURM_JOB_ID} on'
        ' ${SLURM_CPUS_PER_TASK} threads"\n'
    )
    f.write("export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}\n\n")

    for num in num_list:
        job_name = run_name + "_v" + str(num).zfill(3)
        if n_run > -1:
            job_name += "_r" + str(n_run).zfill(3)
        f.write('echo "starting simcore job ' + job_name + '"\n')
        f.write(
            "singularity run simcore_latest.sif ./simcore.exe "
            "" + job_name + "_params.yaml\n"
        )
    f.close()


def randomizeSummitJobs(
    run_name, num_jobs, num_cores, time=24, email="none", n_runs=1, jobs_per_core=1
):
    """ Function inputs: run_name (str), num_jobs (int), num_cores (int), time
    in hours (int), email (str)

        Function behavior: writes num_cores summit batch files with approx
        num_jobs/num_cores jobs to run on each core. Assumes parameter files
        exist with format  <run_name>_v<run_number>.yaml.
    """
    for run in range(n_runs):
        num_list = list(range(num_jobs))
        random.shuffle(num_list)
        # jobs_per_core = int(np.ceil(num_jobs / num_cores))
        n_job_files = int(np.ceil(num_jobs / jobs_per_core))
        for job in range(n_job_files):
            job_nums = num_list[jobs_per_core * job : jobs_per_core * (job + 1)]
            if n_runs == 1:
                writeMultiJob(run_name, job_nums, num_cores, job, time, email)
            else:
                job_num = run * n_job_files + job
                writeMultiJob(run_name, job_nums, num_cores, job_num, time, email, run)


def printHelp():
    """ Prints expected arguments for summit_job_randomizer.py """
    print(
        "=========================================================\n"
        "batch-randomizer.py: generate slurm job files for simcore\n"
        "=========================================================\n"
        "The expected format is:\n"
    )
    print(
        "  python batch-randomizer.py run_name n_sims n_cores time "
        "n_runs email n_sims_per_job\n"
    )
    print(
        "with arguments:\n"
        "  run_name (str): simulation name prefix\n"
        "  n_sims (int): number of sims with unique params\n"
        "  n_cores (int): number of cpus to use per sim\n"
        "  time (int, optional): time in hours per job, default 24\n"
        "  n_runs (int, optional): number of unique seeds per sim, default 1\n"
        "  email (str, optional): user email, default none\n"
        "  n_sims_per_job (int, optional): number of sims to launch per slurm job,"
        " default 1\n"
    )


if __name__ == "__main__":
    num_args = len(sys.argv)
    if num_args == 1:
        printHelp()
        sys.exit(1)
    if num_args < 4 or num_args > 8:
        print("Did not receive correct number of arguments.\n")
        printHelp()
        sys.exit(1)
    try:
        run_name = sys.argv[1]
        num_jobs = int(sys.argv[2])
        num_cores = int(sys.argv[3])
        time = 24
        if num_args > 4:
            time = int(sys.argv[4])
        n_runs = 1
        if num_args > 5:
            n_runs = int(sys.argv[5])
        email = "none"
        if num_args > 6:
            email = str(sys.argv[6])
        jobs_per_core = 1
        if num_args > 7:
            jobs_per_core = int(sys.argv[7])
    except ValueError:
        print("Arguments received are not in the correct format.\n")
        printHelp()
        sys.exit(1)

    print(
        "randomizing job with parameters:\n"
        "  run_name: " + str(run_name) + "\n"
        "  n_sims: " + str(num_jobs) + "\n"
        "  n_cores: " + str(num_cores) + "\n"
        "  time: " + str(time) + "\n"
        "  n_runs: " + str(n_runs) + "\n"
        "  email: " + email + "\n"
        "  n_sims_per_job: " + str(jobs_per_core) + "\n"
    )
    randomizeSummitJobs(
        run_name, num_jobs, num_cores, time, email, n_runs, jobs_per_core
    )
    sys.exit(0)
