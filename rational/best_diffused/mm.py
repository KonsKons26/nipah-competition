import argparse
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pdbfixer import PDBFixer
from sys import stdout

# ==============================================================================
# USER OPTIONS
# ==============================================================================

# Input file
argparser = argparse.ArgumentParser()
argparser.add_argument(
    "-n",
    "--name",
    type=str,
    required=True,
    help="Base name for input/output files."
)

args = argparser.parse_args()
name = args.name

if name.endswith(".pdb"):
    name = name[:-4]

# Force field
ffs = [
    "amber14-all.xml",
    "amber14/tip3pfb.xml"
]

# Model the system
pH = 7.0
padding = 1.5 * unit.nanometers
ionic_strength = 0.1 * unit.molar
add_extra_particles = False

# Create the system
nonbonded_method = app.PME
switch_distance = 0.9 * unit.nanometers
constraints = app.HBonds

# Integrator setup
temperature = 300 * unit.kelvin
friction = 1 / unit.picosecond
timestep = 2 * unit.femtoseconds
integrator = mm.LangevinMiddleIntegrator(
    temperature,
    friction,
    timestep,
)

# Energy minimization
energy_minimization_time = None  # in picoseconds, None to converge
save_minimized_structure = True

# Equilibration
type_of_ensemble = "BOTH"  # "NVT", "NPT", "BOTH"
nvt_time = 100 * unit.picoseconds
npt_time = 100 * unit.picoseconds
barostat_time = 0.05 * unit.picoseconds
track_equilibration_progress = True
save_equilibration_progress = True
save_equilibrated_structure = True

# Production run
production_time = 1 * unit.nanoseconds
report_time = 10 * unit.picoseconds
trajectory_time = 10 * unit.picoseconds
output_name = name + "_production"


# ==============================================================================
# INITIAL SETUP
# ==============================================================================

print("1. Initializing Setup...")

# # Platform selection
# try:
#     platform = mm.Platform.getPlatformByName('CUDA')
#     properties = {'CudaPrecision': 'mixed'}
#     print("   > Using CUDA platform")
# except:
#     try:
#         platform = mm.Platform.getPlatformByName('OpenCL')
#         properties = {}
#         print("   > Using OpenCL platform")
#     except:
#         platform = mm.Platform.getPlatformByName('CPU')
#         properties = {}
#         print("   > Using CPU platform")

platform = mm.Platform.getPlatformByName('OpenCL')
properties = {}

# Force field
forcefield = app.ForceField(*ffs)

# Input file
file_name = name + ".pdb"

print(f"   > Fixing topology for {file_name} using PDBFixer...")
fixer = PDBFixer(filename=file_name)

# Fix Topology
# 1. Detect missing residues (necessary for standard terminals)
fixer.findMissingResidues()
# 2. Detect missing atoms (heavy atoms + OXT)
fixer.findMissingAtoms()
# 3. Add the missing atoms
fixer.addMissingAtoms()

# Model the system
print("   > Modelling the system (adding hydrogens and solvent)...")
# Initialize Modeller with the fixed topology and positions
modeller = app.Modeller(
    fixer.topology, fixer.positions
)
modeller.deleteWater()
modeller.addHydrogens(
    forcefield,
    pH=pH,
)
modeller.addSolvent(
    forcefield,
    padding=padding,
    ionicStrength=ionic_strength,
)
if add_extra_particles:
    modeller.addExtraParticles(forcefield)

# Create system
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=nonbonded_method,
    switchDistance=switch_distance,
    constraints=constraints,
    rigidWater=True,
)

# Simulation setup
simulation = app.Simulation(
    modeller.topology,
    system,
    integrator,
    platform,
    properties,
)
simulation.context.setPositions(modeller.positions)


# ==============================================================================
# ENERGY MINIMIZATION
# ==============================================================================

print("\n2. Minimizing Energy...")

# Convert energy minimization times to number of steps
if energy_minimization_time is None:
    energy_minimization_steps = 0
    print(f"   > Minimizing until convergence...")
else:
    energy_minimization_steps = int(energy_minimization_time / timestep)
    print(f"   > Minimizing for {energy_minimization_steps} steps...")

simulation.minimizeEnergy(maxIterations=energy_minimization_steps)

if save_minimized_structure:
    print(f"   > Saving minimized structure to {name + '_minimized.pdb'}...")
    positions = simulation.context.getState(getPositions=True).getPositions()
    with open(name + "_minimized.pdb", "w") as f:
        app.PDBFile.writeFile(
            simulation.topology,
            positions,
            f,
        )

print("   > Energy minimization complete.")

# ==============================================================================
# EQUILIBRATION
# ==============================================================================
print("\n3. Running Equilibration...")

# Convert equilibration times to number of steps
nvt_steps = int(nvt_time / timestep)
npt_steps = int(npt_time / timestep)

# Set initial velocities
simulation.context.setVelocitiesToTemperature(temperature)

# NVT Equilibration
if type_of_ensemble.upper() in ["NVT", "BOTH"]:
    print(f"   > Running NVT equilibration for {nvt_time}...")

    # Reporters for tracking equilibration progress
    if track_equilibration_progress:
        simulation.reporters.append(
            app.StateDataReporter(
                stdout,
                1000,
                step=True,
                temperature=True,
                potentialEnergy=True,
                kineticEnergy=True,
                totalEnergy=True,
                volume=True,
                density=True,
                progress=True,
                remainingTime=True,
                speed=True,
                totalSteps=nvt_steps,
                separator="\t",
            )
        )

    # Reporter for saving equilibration progress
    if save_equilibration_progress:
        simulation.reporters.append(
            app.StateDataReporter(
                f"{name}_equilibration_nvt.log",
                1000,
                step=True,
                time=True,
                temperature=True,
                potentialEnergy=True,
                volume=True,
                density=True,
                separator="\t",
            )
        )

    simulation.step(nvt_steps)
    simulation.reporters.clear()

# NPT Equilibration
if type_of_ensemble.upper() in ["NPT", "BOTH"]:
    barostat_steps = int(barostat_time / timestep)
    print(f"   > Adding barostat and running NPT equilibration for {npt_time}...")
    barostat = mm.MonteCarloBarostat(
        1.0 * unit.atmosphere,
        temperature,
        barostat_steps,
    )
    system.addForce(barostat)
    simulation.context.reinitialize(preserveState=True)

    # Reporters for tracking equilibration progress
    if track_equilibration_progress:
        simulation.reporters.append(
            app.StateDataReporter(
                stdout,
                1000,
                step=True,
                temperature=True,
                potentialEnergy=True,
                kineticEnergy=True,
                totalEnergy=True,
                volume=True,
                density=True,
                progress=True,
                remainingTime=True,
                speed=True,
                totalSteps=npt_steps,
                separator="\t",
            )
        )

    # Reporter for saving equilibration progress
    if save_equilibration_progress:
        simulation.reporters.append(
            app.StateDataReporter(
                f"{name}_equilibration_npt.log",
                1000,
                step=True,
                time=True,
                temperature=True,
                potentialEnergy=True,
                volume=True,
                density=True,
                separator="\t",
            )
        )

    simulation.step(npt_steps)
    simulation.reporters.clear()

# Save equilibrated structure
if save_equilibrated_structure:
    print(f"   > Saving equilibrated structure to {name + '_equilibrated.pdb'}...")
    positions = simulation.context.getState(getPositions=True).getPositions()
    with open(name + "_equilibrated.pdb", "w") as f:
        app.PDBFile.writeFile(
            simulation.topology,
            positions,
            f,
        )

print("   > Equilibration complete.")


# ==============================================================================
# PRODUCTION RUN
# ==============================================================================
print("\n4. Starting Production Run...")

# Convert production times to number of steps
total_steps = int(production_time / timestep)
report_steps = int(report_time / timestep)
trajectory_steps = int(trajectory_time / timestep)

# Reporter for the log file
simulation.reporters.append(
    app.StateDataReporter(
        f"{output_name}.log",
        report_steps,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=True,
        density=True,
        speed=True,
        progress=True,
        remainingTime=True,
        separator="\t",
        totalSteps=total_steps,
    )
)
# Reporter for the trajectory file
simulation.reporters.append(
    app.DCDReporter(
        output_name + ".dcd",
        trajectory_steps,
    )
)
# Reporter for the terminal (stdout) (more concise)
simulation.reporters.append(app.StateDataReporter(
    stdout,
    report_steps,
    progress=True,
    step=True,
    potentialEnergy=True,
    temperature=True,
    speed=True,
    totalSteps=total_steps,
    separator="\t",
))

print(f"   > Running production for {production_time}...")
print(f"   > Output will be saved to {output_name}.log and {output_name}.dcd")

simulation.step(total_steps)

print("\n5. Simulation Finished Successfully!")

# Save final structure
positions = simulation.context.getState(getPositions=True).getPositions()
with open(output_name + "_final.pdb", "w") as f:
    app.PDBFile.writeFile(
        simulation.topology,
        positions,
        f,
    )

print(f"   > Final structure saved to {output_name}_final.pdb")