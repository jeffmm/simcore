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
    f.write('echo "starting simcore job ' + job_name + '"\n\n')
    f.write("time ./simcore " + job_name + "_params.yaml\n\n")
    f.close()


def writeMultiJob(run_name, num_list, core_num=0, time=4):
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
        core_num, int
    ), "writeMultiJob did not receive a core_num as int type"
    assert isinstance(
        num_list, list
    ), "writeMultiJob did not receive a num_list as list type"
    for num in num_list:
        assert isinstance(
            num, int
        ), "writeMultiJob did not receive a list of nits for num_list"
    job_name_core = run_name + "_j" + str(core_num).zfill(2)
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
    for num in num_list:
        job_name = run_name + "_v" + str(num).zfill(3)
        f.write('echo "starting simcore job ' + job_name + '"\n')
        f.write("time ./simcore " + job_name + "_params.yaml\n")
    f.write("\n")
    f.close()


def randomizeSummitJobs(run_name, num_jobs, num_cores, time=4):
    """ Function inputs: run_name (str), num_jobs (int), num_cores (int), time
    in hours (int)

        Function behavior: writes num_cores summit batch files with approx
        num_jobs/num_cores jobs to run on each core. Assumes parameter files
        exist with format  <run_name>_v<run_number>.yaml.
    """
    num_list = list(range(num_jobs))
    random.shuffle(num_list)
    jobs_per_core = int(np.ceil(num_jobs / num_cores))
    for core in range(num_cores):
        job_nums = num_list[jobs_per_core * core : jobs_per_core * (core + 1)]
        writeMultiJob(run_name, job_nums, core, time)


def printHelp():
    """ Prints expected arguments for summit_job_randomizer.py """
    print("Expected format:\n")
    print("  python summit_job_randomizer.py run_name num_jobs num_cores (time)\n")
    print(
        "  with arguments: run_name (str), num_jobs (int), num_cores (int), time in "
        "hours (int, optional)\n"
    )


if __name__ == "__main__":
    num_args = len(sys.argv)
    if num_args < 4 or num_args > 5:
        print("Did not receive correct number of arguments.\n")
        printHelp()
        sys.exit(1)
    try:
        run_name = sys.argv[1]
        num_jobs = int(sys.argv[2])
        num_cores = int(sys.argv[3])
        time = 4
        if num_args == 5:
            time = int(sys.argv[4])
    except TypeError:
        print("Arguments received are not in the correct format.\n")
        printHelp()
        sys.exit(1)

    randomizeSummitJobs(run_name, num_jobs, num_cores, time)
    sys.exit(0)
