# This is an example to run MD simulations with the CHARMM36m force field
# with either CHARMM tip3p water or an example modified alternative water 
# model (see Huang et al, Nat. Methods, 14, 71, 2017)
# CHARMM psf and crd files are read, and dcd files are written
# This input can be run with "python run.py 0 5"
# The environment variable OPENMM_INCLUDE_PATH, OPENMM_LIB_PATH,
# OPENMM_PLUGIN_DIR, CUDA_HOME, OPENMM_CUDA_COMPILER might need to be set.
# Jing Huang (jing.huang.ff@gmail.com)

from simtk import openmm as mm
from simtk.openmm import app
from simtk.unit import *
from sys import stdout, exit, stderr
from subprocess import call

start=int(sys.argv[1]) # the number of dcd's to run
stop=int(sys.argv[2]) # the number of dcd's to stop
jobname='csptm.c36m'

prefix='/tmp/' # temporary directory to write dcd

# Read the CHARMM PSF file
psf = app.CharmmPsfFile('csptm.psf')

boxsize=6.9 # boxsize in nanometer
psf.setBox(boxsize*nanometer,boxsize*nanometer,boxsize*nanometer)

# Get the coordinates from the CHARMM crd file
crd = app.CharmmCrdFile('csptm.crd') 
# pdb = app.PDBFile('csptm.pdb')  # in case a pdb file instead of crd is used

# Load the parameter set.
lib = 'toppar/' # path of force field topology and parameter files
params = app.CharmmParameterSet(lib+'top_all36_prot.rtf', lib+'par_all36m_prot.prm', lib+'toppar_water_ions.str')
# params = app.CharmmParameterSet(lib+'top_all36_prot.rtf', lib+'par_all36m_prot.prm', lib+'toppar_n010water_ions.str') # in case alternative water model is used to sample more extended states of IDPs

platform = mm.Platform.getPlatformByName('CUDA') # make sure the GPU is used

system = psf.createSystem(params, nonbondedMethod=app.PME, nonbondedCutoff=1.2*nanometer, switchDistance=1.0*nanometer, ewaldErrorTolerance = 0.0001, constraints=app.HBonds)

system.addForce(mm.AndersenThermostat(300*kelvin, 1/picosecond)) # thermostat
system.addForce(mm.MonteCarloBarostat(1*bar, 300*kelvin))        # barostat
integrator = mm.VerletIntegrator(0.002*picoseconds)        # 2 fs time step

simulation = app.Simulation(psf.topology, system, integrator, platform)

simulation.context.setPositions(crd.positions) 
# simulation.context.setPositions(pdb.getPositions()) # in case the pdb file is used

# simulation.minimizeEnergy(maxIterations = 100) # minimization

if start > 0: # restart from the positions and velocities
    restart=str(start-1)
    #with open(jobname+'.'+restart+'.rst', 'r') as f:  # the .rst file should only be used if .chk restart file cannot be loaded
    #    simulation.context.setState(mm.XmlSerializer.deserialize(f.read()))
    with open(jobname+'.'+restart+'.chk', 'rb') as f:
        simulation.context.loadCheckpoint(f.read())

nsavcrd=50000 # save coordinates every 100 ps
nstep=50000000 # write dcd files every 100 ns
nprint=5000000 # report every 10 ns

for ii in range(start,stop+1):
    dcd=app.DCDReporter(prefix+jobname+'.'+str(ii)+'.dcd', nsavcrd)
    firstdcdstep = ii*nstep + nsavcrd
    while (firstdcdstep > 2000000000): # reset frame number to avoid charmm overfloat
        firstdcdstep -= 2000000000
    dcd._dcd = app.DCDFile(dcd._out, simulation.topology, simulation.integrator.getStepSize(), firstdcdstep, nsavcrd) # charmm doesn't like first step to be 0
    simulation.reporters.append(dcd)
    simulation.reporters.append(app.StateDataReporter(jobname+'.'+str(ii)+'.out', nprint, step=True, kineticEnergy=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, speed=True))
    simulation.step(nstep)
    simulation.reporters.pop()
    simulation.reporters.pop()
    dcd._out.close() # close the dcd file to make sure the buffers are written.
   
    # write restart file
    state = simulation.context.getState( getPositions=True, getVelocities=True )
    with open(jobname+'.'+str(ii)+'.rst', 'w') as f:
        f.write(mm.XmlSerializer.serialize(state))
    with open(jobname+'.'+str(ii)+'.chk', 'wb') as f:
        f.write(simulation.context.createCheckpoint())

    # if dcd is written in a tmp direcotry, move it back to home
    cmdstr='mv '+prefix+jobname+'.'+str(ii)+'.dcd .'
    call(cmdstr, shell=True)


