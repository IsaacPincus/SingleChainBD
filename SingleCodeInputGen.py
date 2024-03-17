import math

# Function to set up subfile, actual meat of the code is at the bottom.
def create_subfile(RunName, walltime, mem, nproc, cluster, RunTime=1000, TimesRestart=10):
    """ Creates a submission script for Gadi
    RunName should be a string
    walltime is in minutes
    mem is in GB
    nproc is number of processors to be used
    """

    hours = math.floor(walltime/60)
    mins = walltime%60

    if cluster=="Gadi":
        if mem >= nproc*4:
            print("more mem per processor than allowed, you will be overcharged")

        with open("submit_script.sh", "w") as file:
            file.write("#!/bin/bash\n")
            file.write("#PBS -N " + RunName + "\n")
            file.write("#PBS -P g16\n")
            file.write("#PBS -q normal\n")
            file.write("#PBS -l walltime=" + str(hours) + ":" + str(mins) + ":00\n")
            file.write("#PBS -l mem=" + str(mem) + "GB\n")
            file.write("#PBS -l ncpus=" + str(nproc) + "\n")
            file.write("#PBS -l wd\n")
            file.write("#PBS -M isaac.pincus@monash.edu\n")
            file.write("#PBS -m abe\n")
            file.write("\n")
            file.write("source gadi_modules.sh\n")
            #file.write("echo $HOSTNAME \n")
            file.write("\n")
            file.write("mpirun ./sens >> output.dat \n")
            file.write("\n")
            file.write("times_restarted=0\n")
            file.write("while [ $times_restarted -lt "+str(TimesRestart)+" ]\n")
            file.write("do\n")

            file.write("runtime="+str(RunTime)+"\n")
            file.write("""x=$(sed -rn 's/final time is:\ *(.*)/\\1/p' output.dat)
            y=$(echo $x | awk '{print $NF}')
            new_time=$(python3 -c \"print($runtime+$y)\")
            echo new time is $new_time
            echo old time is $y
            sed -i \"s/Tpr.*/Tpr $new_time/\" inputc.dat

            x=$(sed -rn 's/Restart should be:\ *(.*)/\\1/p' output.dat)
            y=$(echo $x | awk '{print $NF}')
            echo New restart is $y
            sed -i \"s/Restart.*/Restart $y/\" inputc.dat

            times_restarted=$(python3 -c \"print($times_restarted+1)\")
            echo $times_restarted

            mpirun ./sens >> output.dat
            done
            \n""")

    elif cluster=="Massive":
        with open("sub_script.pbs", "w") as file:
            file.write("#!/bin/bash\n")

            file.write("#SBATCH --job-name="  + RunName + "\n")
            file.write("#SBATCH --account=oy14\n")
            file.write("#SBATCH --time=" + str(hours) + ":" + str(mins) + ":00\n")
            file.write("#SBATCH --mem=" + str(mem) + "G\n")
            file.write("#SBATCH --ntasks=" + str(nproc) + "\n")
            file.write("#SBATCH --ntasks-per-node=" + str(nproc) + "\n")
            file.write("#SBATCH --cpus-per-task=1\n")
            file.write("#SBATCH --output=output.dat\n")
            file.write("#SBATCH --error=errset.err\n")
            file.write("#SBATCH --mail-type=ALL\n")
            file.write("#SBATCH --mail-user=isaac.pincus@monash.edu\n")
            file.write("\n")
            file.write("source massive_modules.sh\n")
            file.write("echo $HOSTNAME\n")
            file.write("\n")
            file.write("mpirun sens\n")
            file.write("times_restarted=0\n")
            file.write("while [ $times_restarted -lt "+str(TimesRestart)+" ]\n")
            file.write("do\n")

            file.write("runtime="+str(RunTime)+"\n")
            file.write("""x=$(sed -rn 's/final time is:\ *(.*)/\\1/p' output.dat)
            y=$(echo $x | awk '{print $NF}')
            new_time=$(python -c \"print($runtime+$y)\")
            echo new time is $new_time
            echo old time is $y
            sed -i \"s/Tpr.*/Tpr $new_time/\" inputc.dat

            x=$(sed -rn 's/Restart should be:\ *(.*)/\\1/p' output.dat)
            y=$(echo $x | awk '{print $NF}')
            echo New restart is $y
            sed -i \"s/Restart.*/Restart $y/\" inputc.dat

            times_restarted=$(python -c \"print($times_restarted+1)\")
            echo $times_restarted

            mpirun sens
            done
            \n""")
            file.write('sacct --units=G --format="JobID,JobName,CPUTime,MaxRSS,NCPUs"\n')

    elif cluster=="Monarch":
        with open("sub_script.pbs", "w") as file:
            file.write("#!/bin/bash\n")
            file.write("#SBATCH --job-name="  + RunName + "\n")
            file.write("#SBATCH --time=" + str(hours) + ":" + str(mins) + ":00\n")
            if hours<24:
                file.write("#SBATCH --partition=short,comp\n")
            else:
                file.write("#SBATCH --partition=comp\n")
            file.write("#SBATCH --mem=" + str(mem) + "G\n")
            file.write("#SBATCH --ntasks=" + str(nproc) + "\n")
            file.write("#SBATCH --ntasks-per-node=" + str(nproc) + "\n")
            file.write("#SBATCH --cpus-per-task=1\n")
            file.write("#SBATCH --output=output.dat\n")
            file.write("#SBATCH --error=errset.err\n")
            file.write("#SBATCH --mail-type=ALL\n")
            file.write("#SBATCH --mail-user=isaac.pincus@monash.edu\n")
            file.write("\n")
            file.write("source monarch_modules.sh\n")
            file.write("echo $HOSTNAME\n")
            file.write("\n")
            file.write("mpirun sens\n")
            file.write("times_restarted=0\n")
            file.write("while [ $times_restarted -lt "+str(TimesRestart)+" ]\n")
            file.write("do\n")

            file.write("runtime="+str(RunTime)+"\n")
            file.write("""x=$(sed -rn 's/final time is:\ *(.*)/\\1/p' output.dat)
            y=$(echo $x | awk '{print $NF}')
            new_time=$(python -c \"print($runtime+$y)\")
            echo new time is $new_time
            echo old time is $y
            sed -i \"s/Tpr.*/Tpr $new_time/\" inputc.dat

            x=$(sed -rn 's/Restart should be:\ *(.*)/\\1/p' output.dat)
            y=$(echo $x | awk '{print $NF}')
            echo New restart is $y
            sed -i \"s/Restart.*/Restart $y/\" inputc.dat

            times_restarted=$(python -c \"print($times_restarted+1)\")
            echo $times_restarted

            mpirun sens
            done
            \n""")
            file.write('sacct --units=G --format="JobID,JobName,CPUTime,MaxRSS,NCPUs"\n')
            file.write('sacct --units=G --format="JobID,JobName,CPUTime,MaxRSS,NCPUs"\n')

    elif cluster=='Local':
        with open('run_script.sh', 'w') as file:
            file.write('#!/bin/bash\n')
            file.write('\n')
            file.write('time /usr/bin/mpirun -n ' + str(nproc) + ' ./sens > terminal_out.dat 2>error.dat \n')


# Function to create input script
def create_inputc(SPtype, FlowType, Nbeads, hstar, zstar, dstar, dQ, sigma, gdot,
                  Nsamples, Teq, Tpr,
                  ndelts, dtseq, dtsne, nblock, ntot, tol,
                  EV="noEV",Bens=False, NetCDF=True, Restart=0,
                  VR=False,  phiFile=False, phiSDK=False, Ntrajdone=0,
                  InitialConfiguration="RandomSpherical",
                  BendingPotentialType="NoBendingPotential",
                  NaturalAngleFromFile=False,
                  BendingStiffness=0, NaturalAngleScalarInput=0,
                  COMUpdateOn=True,
                  TrapOneInitialPosition=0, TrapTwoInitialPosition=0,
                  TrapOneStrength=0, TrapTwoStrength=0,
                  TrapTwoVelocity=0,
                  LookupTableOpt=False, LookupTableTol=0, max_gamma=0,
                  min_EV_cutoff=0.7, max_EV_cutoff=1.5, contour_dist_for_EV=1,
                  delSCalcMethod=0, EigsCalcMethod=0, nchebMultiplier=1.0,
                  fdErrMax=2.5e-3, ChebUpdateMethod=0):

    SpringTypeInt = SPtype
    FlowTypeInt = FlowType
    EVOptionInt = EV
    InitialConfigInt = InitialConfiguration
    BendingPotentialInt = BendingPotentialType

    with open("inputc.dat", "w") as file:
        file.write("SpType\t" + str(SpringTypeInt) + "\n")
        file.write("LookupTableOpt\t" + str(int(LookupTableOpt==True)) + "\n")
        file.write("LookupTableTolerance\t" + str(LookupTableTol) + "\n")
        file.write("MaxGamma\t" + str(max_gamma) + "\n")
        file.write("FlowType\t" + str(FlowTypeInt) + "\n")
        file.write("NBeads\t" + str(Nbeads) + "\n")
        file.write("hstar\t" + str(hstar) + "\n")
        file.write("EV\t" + str(EVOptionInt) + "\n")
        file.write("zstar\t" + str(zstar) + "\n")
        file.write("dstar\t" + str(dstar) + "\n")
        file.write("contour_dist_for_EV\t" + str(contour_dist_for_EV) + "\n")
        file.write("min_EV_cutoff\t" + str(min_EV_cutoff) + "\n")
        file.write("max_EV_cutoff\t" + str(max_EV_cutoff) + "\n")
        file.write("BendingPotentialType\t" + str(BendingPotentialInt) + "\n")
        file.write("BendingStiffness\t" + str(BendingStiffness) + "\n")
        file.write("NaturalAngleFromFile\t" + str(int(NaturalAngleFromFile==True)) + "\n")
        file.write("natural_angle_scalar_input\t" + str(NaturalAngleScalarInput) + "\n")
        file.write("sqrtb\t" + str(dQ) + "\n")
        file.write("sigma\t" + str(sigma) + "\n")
        file.write("gdots\t" + str(gdot) + "\n")
        file.write("bens\t" + str(int(Bens==True)) + "\n")
        file.write("NetCDF\t" + str(int(NetCDF==True)) + "\n")
        file.write("Restart\t" + str(Restart) + "\n")
        file.write("variance_reduction\t" + str(int(VR==True)) + "\n")
        file.write("COM_update_on\t" + str(int(COMUpdateOn==True)) + "\n")
        file.write("InitialConfiguration\t" + str(InitialConfigInt) + "\n")
        file.write("delSCalcMethod\t" + str(delSCalcMethod) + "\n")
        file.write("EigsCalcMethod\t" + str(EigsCalcMethod) + "\n")
        file.write("ChebUpdateMethod\t" + str(ChebUpdateMethod) + "\n")
        file.write("nchebMultiplier\t" + str(nchebMultiplier) + "\n")
        file.write("fd_err_max\t" + str(fdErrMax) + "\n")
        file.write("TrapOneInitialPosition\t" + str(TrapOneInitialPosition) + "\n")
        file.write("TrapTwoInitialPosition\t" + str(TrapTwoInitialPosition) + "\n")
        file.write("TrapOneStrength\t" + str(TrapOneStrength) + "\n")
        file.write("TrapTwoStrength\t" + str(TrapTwoStrength) + "\n")
        file.write("TrapTwoVelocity\t" + str(TrapTwoVelocity) + "\n")
        #file.write("emax\t" + str(emax) + "\n")
        file.write("Nsamples\t" + str(Nsamples) + "\n")
        file.write("phiFromFile\t" + str(int(phiFile==True)) + "\n")
        file.write("phi_data\t" + str(int(phiSDK==True)) + "\n")
        file.write("Teq\t" + str(Teq) + "\n")
        file.write("Tpr\t" + str(Tpr) + "\n")
        file.write("Ntrajdone\t" + str(Ntrajdone) + "\n")

        file.write("ndelts\t" + str(ndelts) + "\n")
        for count,ele in enumerate(dtseq):
            file.write("dtseq\t" + str(ele) + "\n" +
                       "dtsne\t" + str(dtsne[count]) + "\n" +
                       "nblock\t" + str(nblock) + "\n" +
                       "ntot\t" + str(ntot) + "\n" +
                       "tol\t" + str(tol[count]) + "\n")
