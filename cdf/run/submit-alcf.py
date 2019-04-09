#!/usr/bin/env python
#COBALT -A datascience
###COBALT -q debug-flat-quad
###COBALT -n 8
###COBALT -t 5

# Load python modules
import subprocess
import argparse
import os
import sys
print("Using Python version: "+str(sys.version_info[0]))

# Machine Specific
machine   = "theta"
srcroot   = "/home/zamora/hdf5_root_dir/xgitlabrepos/pnetcdf_vol"
fsroot    = "/projects/datascience/rzamora"
ppn_in    = 32
nranksmin = 1024

#machine   = "other"
#srcroot   = "/Users/rzamora/remote-fs/theta_pnetcdf_vol"
#fsroot    = "/Users/rzamora/remote-fs"
#ppn       = 4
#nranksmin = 2

# Important Benchmarking Settings
benchname = "hdf5-cdfvl-20190904"
runroot   = fsroot + "/benchscratch"
benchroot = srcroot + "/cdf"
outroot   = benchroot + "/run/results/" + benchname
execname  = benchroot + "/bin/vol_test.ex"
lfs_count = 52    # Number of Stripes in LUSTRE
lfs_size  = 8     # Size of Stripes in LUSTRE (MB Units)
nodes_in  = 1     # Ignored if machine is "vesta" or "theta"
nocheck   = False # Can turn off data validation in benchmark

# Parse command line inputs
parser = argparse.ArgumentParser()
parser.add_argument("--machine", dest="machine", default=machine,
                    help="system name -- Available: theta, vesta, other [default="+machine+"]")
parser.add_argument("--execname", dest="execname", default=execname,
                    help="Path to Exerciser executable [default="+execname+"]")
parser.add_argument("--ppn", dest="ppn", type=int, default=ppn,
                    help="Processes to use per node [default="+str(ppn_in)+"]")
parser.add_argument("--lfs_count", dest="lfs_count", type=int, default=lfs_count,
                    help="Lustre Stripe Count [default="+str(lfs_count)+"]")
parser.add_argument("--lfs_size", dest="lfs_size", type=int, default=lfs_size,
                    help="Lustre Stripe Size (Units=MB) [default="+str(lfs_size)+"]")
#parser.add_argument("--ccio", dest="ccio", action="store_true", default=ccio,
#                    help="Using the CCIO version of HDF5 (Set CCIO Env Vars) [default="+str(ccio)+"]")
args = parser.parse_args()
machine   = args.machine
execname  = args.execname
lfs_count = args.lfs_count
lfs_size  = args.lfs_size
ppn       = args.ppn

# Check that machine is supported
if not machine in ["theta", "vesta", "other"]:
    print("Error - machine "+machine+" not recongnized")
    sys.exit(1)

# Define Env vars that we wont change here
envs_const = [ ]
if machine in ["theta", "vesta"]:
    nodes_in = int(os.environ['COBALT_JOBSIZE'])
if machine == "theta":
    # Allow module load/swap/list etc:
    execfile(os.environ['MODULESHOME']+'/init/python.py')
    #os.environ['MPICH_MPIIO_HINTS'] = '*:cray_cb_write_lock_mode=1'
    os.environ['MPICH_NEMESIS_ASYNC_PROGRESS'] = 'ML'
    os.environ['MPICH_MAX_THREAD_SAFETY'] = 'multiple'
    module('load','cray-parallel-netcdf')
elif machine == "vesta":
    envs_const.append("BGLOCKLESSMPIO_F_TYPE=0x47504653")

# Env Var Helper Funciton
def export_envs( envs_dyn ):

    for env in envs_dyn:
        env_split  = env.split("=")
        env_name   = env_split[0]
        env_value  = env_split[1]
        subprocess.call(["echo","Setting "+env_name+" to "+env_value+"."], stdout=outf)
        os.environ[ env_name ] = env_value

# Run-command Helper Funciton
def get_runjob_cmd( envs_dyn, nranks, ppn ):

    if machine == "vesta":
        cmd = ["runjob"]
        cmd.append("--np");    cmd.append(str(nranks))
        cmd.append("-p");      cmd.append(str(ppn))
        cmd.append("--block"); cmd.append(os.environ['COBALT_PARTNAME'])

        # Environment variables added here
        for env in envs_const:
            cmd.append("--envs");  cmd.append(env)
        for env in envs_dyn:
            cmd.append("--envs");  cmd.append(env)

        # Benchmark
        cmd.append(":"); cmd.append(execname)

    elif machine == "theta":
        export_envs( envs )
        # Define env and part of aprun command that wont change
        #os.environ["MPICH_MPIIO_TIMERS"]="1"
        #os.environ["MPICH_MPIIO_STATS"]="1"
        #os.environ["MPICH_MPIIO_AGGREGATOR_PLACEMENT_DISPLAY"]="1"
        #os.environ["PMI_LABEL_ERROUT"]="1"
        #os.environ["MPICH_MPIIO_CB_ALIGN"]="2"
        #os.environ["MPICH_MPIIO_HINTS"]="*:romio_ds_write=disable"
        subprocess.Popen('ulimit -c unlimited', shell=True)
        cmd = ["aprun"]
        cmd.append("-n"); cmd.append(str(nranks)); cmd.append("-N"); cmd.append(str(ppn))
        cmd.append("-d"); cmd.append("1"); cmd.append("-j"); cmd.append("1")
        cmd.append("-cc"); cmd.append("depth"); cmd.append(execname)

    else:
        export_envs( envs )
        cmd = ["mpirun"]
        cmd.append("-n"); cmd.append(str(nranks))
        cmd.append(execname)

    return cmd

def print_cmd( cmd ):
    sys.stdout.write('Using Run Command:\n')
    for val in cmd:
        sys.stdout.write('%s ' % (val))
    sys.stdout.write('\n')

# Create top "scratch" directory (runroot):
if not os.path.isdir(runroot): subprocess.call(["mkdir",runroot])

# Create "scratch" directory for current benchmark campaign (runbench):
runbench = runroot + "/" + benchname
if not os.path.isdir(runbench): subprocess.call(["mkdir",runbench])

# Create directory for actual I/O operations (rundir):
rundir = runbench+"/count."+str(lfs_count)+".size."+str(lfs_size)+".nodes."+str(nodes_in)+".ppn."+str(ppn)
if not os.path.isdir(rundir): subprocess.call(["mkdir",rundir])

# Create "output" directory for current benchmark campaign (outroot):
if not os.path.isdir(outroot): subprocess.call(["mkdir",outroot])

# Create directory for actual results (outdir):
outdir = outroot+"/count."+str(lfs_count)+".size."+str(lfs_size)+".nodes."+str(nodes_in)+".ppn."+str(ppn)
if not os.path.isdir(outdir): subprocess.call(["mkdir",outdir])

# Determine the "jobid"
if machine == "other":
    jobid_i = 0
    while os.path.exists( outdir+"/results."+str(jobid_i) ):
        jobid_i += 1
    jobid = str(jobid_i)
else: jobid = os.environ['COBALT_JOBID']

# Go to run directory
os.chdir(rundir)

# Now, run the job (with output written to "outdir/results.<jobid>")
with open(outdir+"/results."+jobid, "a") as outf:

    if machine == "theta":
        # Set lustre stripe properties
        subprocess.call(["lfs","setstripe","-c",str(lfs_count),"-S",str(lfs_size)+"m","."])

    ntrials = 10
    run_all = False
    dimlens = 64
    nranks = ppn * nodes_in

    while nranks >= nranksmin:

        # Independent I/O
        subprocess.call(["echo",""], stdout=outf)
        subprocess.call(["echo","[Independent]: "+str(nranks)+" procs"], stdout=outf)
        for itrial in range(ntrials):
            envs = [ ]
            cmd = list( get_runjob_cmd( envs, nranks, ppn ) );
            cmd.append("--dimlen"); cmd.append(str(dimlens));
            if run_all: cmd.append("--all");
            if nocheck: cmd.append("--nocheck");
            print_cmd(cmd)
            subprocess.call(cmd, stdout=outf)

        # Collective I/O
        subprocess.call(["echo",""], stdout=outf)
        subprocess.call(["echo","[Collective]: "+str(nranks)+" procs"], stdout=outf)
        for itrial in range(ntrials):
            envs = [ ]
            cmd = list( get_runjob_cmd( envs, nranks, ppn ) );
            cmd.append("--dimlen"); cmd.append(str(dimlens));
            cmd.append("--col");
            if run_all: cmd.append("--all");
            if nocheck: cmd.append("--nocheck");
            print_cmd(cmd)
            subprocess.call(cmd, stdout=outf)

        # Divide system in half for next run
        nranks = int(nranks / 2)
        if (nranks % ppn): break

# ---------------------------------------------------------------------------- #
#  Done.
# ---------------------------------------------------------------------------- #

cmd = ["echo","done"]
subprocess.call(cmd)
