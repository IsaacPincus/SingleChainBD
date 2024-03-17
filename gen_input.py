# This file generates the inputc file for FENE-Fraenkel spring simulations
# Values can be input in either Hookean or Rodlike unit systems

import math
import os
import shutil
import numpy as np
import SingleCodeInputGen as gen

# Generic stuff
SimulationNamePrefix = ""
SimulationNameSuffix = ""
#Options are Hookean, FENE, InverseLangevin, Wormlike, Fraenkel, FENEFraenkel
SPTypeDict = {"Hookean": 1, "FENE": 2, "InverseLangevin": 3,
              "Wormlike": 4, "Fraenkel": 5, "FENEFraenkel": 6,
              "WLC_bounded":7}
SpringType = "Hookean"
# Options are SH, UA, PL, UR, PU, PP, EQ
FTDict = {"Equilibrium no EV or HI": 0, "Planar Shear": 1, "Uniaxial Elongational": 2,
          "Planar Elongational": 3, "Uniaxial Extension Relaxation": 4,
          "Periodic Uniaxial Extension": 5, "Periodic Planar Extension": 6,
          "Equilibrium": 7, "BeadPulling": 8}
FlowTypeInput = "Planar Shear"
BendTypeDict = {"NoBendingPotential": 0, "OneMinusCosTheta": 1}
BendingTypeInput = "NoBendingPotential"
BendingStiffness = 0
NaturalAngleFromFile = False
NaturalAngles = np.array([0])
NaturalAngleScalarInput = 0
ICDict = {"Random Spherical":1, "x-axis Aligned":2}
InitialConfigurationOption = "Random Spherical"
NumberOfBeads = 10
COMUpdateOn = True
LookupTableOpt = False
LookupTableTol = 0
MaxGamma = 0
#Options are noEV, Gaussian, LJ, SDK
EVOptionDict = {"noEV": 0, "Gaussian": 1, "LJ": 2, "SDK": 3}
ExcludedVolumeOption = "noEV"
BlockEnsemble = False
NetCDF = True
Restart = 0
VarianceReduction = False
NumberOfSamplesPerTrajectory = 3
phiFromFile = False
phiFromSDK = False
DelSCalcDict = {"Chebyshev": 0, "Cholesky": 1, "ExactSqrt": 2}
delSCalcMethod = DelSCalcDict['Chebyshev']
EigsCalcDict = {"EigsFixman": 0, "EigsExact": 1}
EigsCalcMethod = EigsCalcDict['EigsFixman']
ChebUpdateDict = {"UpdateChebNew": 0, "UpdateChebAddOne": 1}
ChebUpdateMethod = ChebUpdateDict['UpdateChebNew']
nchebMultiplier = 1
fdErrMax = 2.5e-3
# Trajectories already completed is PER PROCESSOR!
NumberOfTrajectoriesAlreadyCompleted = 0
NumberOfTimestepWidths = 1
NoProcessors = 24
TotalTrajectories = 5
RodlikeUnits = False
#mem in GB
gbOfMemory=10
#walltime in minutes
minsWalltime=500
cluster="Local"
if cluster=="Gadi":
    codeLocation = "/scratch/g16/ip9701/Single_Chain/code/Single-Chain-Development/"
    sensLocation = os.path.join(codeLocation,'sens')
    modulesLocation = os.path.join(codeLocation,'gadi_modules.sh')
elif cluster=="Massive":
    codeLocation = "/home/ipincus/oy14/ipincus/code/single-chain-development/"
    sensLocation = os.path.join(codeLocation,'sens')
    modulesLocation = os.path.join(codeLocation,'massive_modules.sh')
elif cluster=="Monarch":
    codeLocation = "/home/ipincus/wm73/ipincus/Single-Chain/Code/single-chain-development/"
    sensLocation = os.path.join(codeLocation, 'sens')
    modulesLocation = os.path.join(codeLocation, 'monarch_modules.sh')
elif cluster=="Local":
    codeLocation = "/home/ipincus/single_chain_code/SingleChainBD/"
    sensLocation = os.path.join(codeLocation, 'sens')

# Parameters in Hookean units
hstar = 0.5
zstar = 0
dstar = 0
# dQ (aka sqrtb) shouldn't be 0 even if you're using a Hookean or Fraenkel spring, else
# you will get divide by 0 errors. This is stupid and I hope I have time to fix
# it properly in the code.
dQ = 1
sigma = 1
gammaDot = np.array([0])
EquilibrationTime = 12
ProductionTime = 20
TrapOneInitialPosition = 0
TrapTwoInitialPosition = 0
TrapOneStrength = 0
TrapTwoStrength = 0
TrapTwoVelocity = 0
contour_dist_for_EV = 1

dtEquilibration = np.array([0.01])
dtProduction = np.array([0.01])
tolerance = np.array([0.00001])

Wi_vals = [0.1]

# Parameters in Rodike Units, comment out if not needed
if RodlikeUnits:
    H = 800
    s = 0.05
    tau = 1
    #tau = 0.0142*(NumberOfBeads+1)**2    #If using Wi
    hstarRod = 0.5
    zstarRod = 0
    dstarRod = 0
    Wi_vals = np.logspace(-1, 3, 9)
    gammaDotRod = Wi_vals/tau
    EqTimeRod = 10
    ProdTimeRod = 30
    TrapOneInitialPositionRod = 0
    TrapTwoInitialPositionRod = 0
    TrapTwoVelocityRod = 0

    #dtEquilibrationRod = np.array([0.001])
    #dtProductionRod = np.array([0.001])
    toleranceRod = np.array([0.001, 0.001])

    #calculations
    sigma = math.sqrt(H)
    dQ = s*sigma
    hstar = (4*sigma*hstarRod)/(3*math.sqrt(math.pi))
    gammaDot = gammaDotRod/(4*H)
    EquilibrationTime = EqTimeRod*4*H
    ProductionTime = ProdTimeRod*4*H
    #dtEquilibration = dtEquilibrationRod*4*H
    #dtProduction = dtProductionRod*4*H
    tolerance = toleranceRod*sigma
    zstar = zstarRod #Energies are both scaled by k_BT
    dstar = dstarRod*sigma
    TrapOneInitialPosition = TrapOneInitialPositionRod*sigma
    TrapTwoInitialPosition = TrapTwoInitialPositionRod*sigma
    TrapTwoVelocity = TrapTwoVelocityRod/(4*sigma)

if ExcludedVolumeOption == "LJ":
    max_EV_cutoff = 1.5*dstar
elif ExcludedVolumeOption == "SDK":
    max_EV_cutoff = 1.82*dstar
else:
    max_EV_cutoff = 1.5*dstar

min_EV_cutoff = 0.7*dstar

# processors etc
TrajectoriesPerProcessor = math.ceil(TotalTrajectories/NoProcessors)
TrueTotalTrajectories = TrajectoriesPerProcessor*NoProcessors

if ((len(dtEquilibration)!=NumberOfTimestepWidths) or
(len(dtProduction)!=NumberOfTimestepWidths) or
(len(tolerance)!=NumberOfTimestepWidths)):
    print("Number of timestep widths doesn't match dt array lengths")


# Actual writing and everything
topDir = os.getcwd()
for item, Wi in enumerate(Wi_vals):
    WiString = "Wi{0:g}".format(Wi)
    WiDir = os.path.join(topDir,WiString)
    try:
        os.mkdir(WiDir)
    except FileExistsError:
        pass

    os.chdir(WiDir)
    shutil.copy2(sensLocation,os.path.join(WiDir, 'sens'))
    #shutil.copy2(os.path.join(codeLocation, 'Matlab_scripts/Calculate_averages_from_trajectories.m'),os.path.join(WiDir, 'sens'))
    if cluster=="Gadi":
        shutil.copy2(modulesLocation,os.path.join(WiDir, 'gadi_modules.sh'))
    elif cluster=="Massive":
        shutil.copy2(modulesLocation,os.path.join(WiDir, 'massive_modules.sh'))
    elif cluster=="Monarch":
        shutil.copy2(modulesLocation,os.path.join(WiDir, 'monarch_modules.sh'))
    elif cluster=="Local":
        print("no modules needed!")

    if Wi<1:
        VarianceReduction=True
    else:
        VarianceReduction=False

    gen.create_subfile(SimulationNamePrefix+
                   WiString+
                   SimulationNameSuffix,
                   minsWalltime,
                   gbOfMemory,
                   NoProcessors,
                   cluster=cluster)

    gen.create_inputc(SPtype=SPTypeDict[SpringType], FlowType=FTDict[FlowTypeInput],
                  Nbeads=NumberOfBeads, hstar=hstar, zstar=zstar, dstar=dstar,
                  dQ=dQ, sigma=sigma, gdot=gammaDot[item], Nsamples=NumberOfSamplesPerTrajectory,
                  Teq=EquilibrationTime, Tpr=ProductionTime,
                  ndelts=NumberOfTimestepWidths, dtseq=dtEquilibration, dtsne=dtProduction,
                  nblock=TrajectoriesPerProcessor, ntot=TrueTotalTrajectories, tol=tolerance,
                  VR=VarianceReduction, EV=EVOptionDict[ExcludedVolumeOption],
                  Bens=BlockEnsemble, NetCDF=NetCDF, Restart=Restart,
                  InitialConfiguration=ICDict[InitialConfigurationOption],
                  BendingPotentialType=BendTypeDict[BendingTypeInput],
                  NaturalAngleFromFile=NaturalAngleFromFile,
                  BendingStiffness=BendingStiffness,
                  NaturalAngleScalarInput=NaturalAngleScalarInput,
                  COMUpdateOn=COMUpdateOn,
                  TrapOneInitialPosition=TrapOneInitialPosition,
                  TrapTwoInitialPosition=TrapTwoInitialPosition,
                  TrapOneStrength=TrapOneStrength,
                  TrapTwoStrength=TrapTwoStrength,
                  TrapTwoVelocity=TrapTwoVelocity,
                  LookupTableOpt=LookupTableOpt,
                  LookupTableTol=LookupTableTol,
                  max_gamma=MaxGamma,
                  max_EV_cutoff=max_EV_cutoff,
                  min_EV_cutoff=min_EV_cutoff,
                  contour_dist_for_EV=contour_dist_for_EV,
                  delSCalcMethod=delSCalcMethod,
                  EigsCalcMethod=EigsCalcMethod,
                  nchebMultiplier=nchebMultiplier,
                  fdErrMax=fdErrMax,
                  ChebUpdateMethod=ChebUpdateMethod)


    if NaturalAngleFromFile:
        with open("natural_angles.txt", "w") as file:
            file.write(str(NaturalAngles))

    os.chdir(topDir)
