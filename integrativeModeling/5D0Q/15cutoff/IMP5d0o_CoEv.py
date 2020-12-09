'''
This is Caitie's adaptation for 5D0Q.
'''
# Imports
from __future__ import print_function
import IMP
import RMF
import IMP.rmf
import IMP.pmi.dof
import IMP.pmi.io.crosslink
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.bayesianem
import IMP.bayesianem.restraint
import IMP.container
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.restraints.basic
import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros
import IMP.pmi.topology
import IMP.pmi.analysis
import IMP.em

import os
import sys


# Identify data files
gmm_data = "/home/cns-mccafferty/Ecoli_CoEv_modeling/CoEv_IMP/5d0q/5d0q.50gmm.txt"

# Topology File
topology_file = "/home/cns-mccafferty/Ecoli_CoEv_modeling/CoEv_IMP/5d0q/topology.txt"
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Here is where the work begins
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# All IMP systems start out with a Model
mdl = IMP.Model()
# Read the topology file for a given state
t = IMP.pmi.topology.TopologyReader(topology_file,pdb_dir='/home/cns-mccafferty/Ecoli_CoEv_modeling/CoEv_IMP/5d0q/', fasta_dir='/home/cns-mccafferty/Ecoli_CoEv_modeling/CoEv_IMP/5d0q/', gmm_dir='/home/cns-mccafferty/Ecoli_CoEv_modeling/CoEv_IMP/5d0q/')
# Create a BuildSystem macro to and add a state from a topology file
bs = IMP.pmi.macros.BuildSystem(mdl)
bs.add_state(t)
# executing the macro will return the root hierarchy and degrees of freedom (dof) objects
root_hier, dof = bs.execute_macro(max_rb_trans=4.0,
max_rb_rot=1.0,
max_bead_trans=4.0,
max_srb_trans=1.0,
max_srb_rot=0.1)
## output to RMF
fname = 'test.rmf'
rh = RMF.create_rmf_file(fname)
IMP.rmf.add_hierarchy(rh, root_hier)
IMP.rmf.save_frame(rh)
#####################################################
##################### RESTRAINTS ####################
#####################################################
# Restraints define functions that score the model based on
# input information.  
#
# Restraint objects are first created in the definition.
# To be evaluated, the restraint object must be add_to_model().
#
# In some cases, sampled parameters for restraints must be added to the DOF
# object
# The output_objects list is used to collect all restraints
# where we want to log the output in the STAT file.
# Each restraint should be appended to this list.
output_objects = []
display_restraints = [] # display as springs in RMF
# -----------------------------
# %%%%% CONNECTIVITY RESTRAINT
#
# Restrains residues/particles that are collected in sequence
# This should be used for any system without an atomic force field (e.g. CHARMM)
# We apply the restraint to each molecule
mols = IMP.pmi.tools.get_molecules(root_hier)
for mol in mols:
    molname=mol.get_name()
    IMP.pmi.tools.display_bonds(mol)
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
    cr.add_to_model()
    cr.set_label(molname)
    output_objects.append(cr)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(mdl))
print("Conn", sf.evaluate(False))
# Excluded Volume Restraint
#  To speed up this expensive restraint, we operate it at resolution 20
ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                         included_objects=root_hier,
                                         resolution=1000)
ev.add_to_model()
output_objects.append(ev)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(mdl))
print("EV", sf.evaluate(True))
# -------------------------
# %%%%% COEVOLUTIONARY RESTRAINT
#
# Incorporating a coevolutioary restraint as a basic restraint

protein1 = [480, 372, 351, 363, 361, 38, 373, 351, 363, 362, 366, 652, 177, 351, 745, 365, 358, 652] #test data
protein2 = [187, 190, 178, 187, 134, 99, 197, 181, 184, 181, 185, 191, 114, 204, 205, 181, 139, 140]
evscore = [1.73, 1.4, 1.37, 1.34, 1.3, 0.96, 0.96, 0.81, 0.75, 0.75, 0.65, 0.65, 0.62, 0.62, 0.62, 0.61, 0.58, 0.57]
for i in range(len(protein1)): #loops through a list of coev pairs and their scores
    coevr2 = IMP.pmi.restraints.basic.DistanceRestraint(representation=None,
                     tuple_selection1=(protein1[i],15,"P0A940",0),
                     tuple_selection2=(protein2[i],25,"P0AC02",0),
                     distancemin=0,
                     distancemax=20.4,
                     resolution=1.0,
                     kappa=1.0,
                     root_hier=root_hier,
                     weight=evscore[i])

    coevr2.add_to_model()
    output_objects.append(coevr2)

protein1 = [375, 481, 376, 621] #test data
protein2 = [33, 29, 33, 71]
evscore = [0.6, 0.57, 0.51, 0.51]
for i in range(len(protein1)): #loops through a list of coev pairs and their scores
    coevr3 = IMP.pmi.restraints.basic.DistanceRestraint(representation=None,
                     tuple_selection1=(protein1[i],15,"P0A940",0),
                     tuple_selection2=(protein2[i],25,"P0A937",0),
                     distancemin=0,
                     distancemax=25.3,
                     resolution=1.0,
                     kappa=1.0,
                     root_hier=root_hier,
                     weight=evscore[i])

    coevr3.add_to_model()
    output_objects.append(coevr3)


protein1 = [109, 95, 44] #test data
protein2 = [184, 82, 172]
evscore = [0.38, 0.38, 0.38]
for i in range(len(protein1)): #loops through a list of coev pairs and their scores
    coevr1 = IMP.pmi.restraints.basic.DistanceRestraint(representation=None,
                     tuple_selection1=(protein1[i],15,"P0A903",0),
                     tuple_selection2=(protein2[i],25,"P0AC02",0),
                     distancemin=0,
                     distancemax=29.7,
                     resolution=1.0,
                     kappa=1.0,
                     root_hier=root_hier,
                     weight=evscore[i])

    coevr1.add_to_model()
    output_objects.append(coevr1)

protein1 = [135] #test data
protein2 = [44]
evscore = [0.42]
for i in range(len(protein1)): #loops through a list of coev pairs and their scores
    coevr1 = IMP.pmi.restraints.basic.DistanceRestraint(representation=None,
                     tuple_selection1=(protein1[i],15,"P0A903",0),
                     tuple_selection2=(protein2[i],25,"P0A937",0),
                     distancemin=0,
                     distancemax=29,
                     resolution=1.0,
                     kappa=1.0,
                     root_hier=root_hier,
                     weight=evscore[i])

    coevr1.add_to_model()
    output_objects.append(coevr1)

protein1 = [237, 230, 230, 156, 143, 199, 184, 208, 61, 208, 188, 137, 124, 154, 143, 203, 146, 199, 192, 192, 84] #test data
protein2 = [68, 66, 74, 97, 59, 67, 32, 66, 36, 58, 32, 77, 61, 40, 39, 67, 27, 39, 63, 83, 84]
evscore = [2.11, 1.4, 1.37, 1.23, 1.22, 1.21, 1.2, 1.19, 1.15, 1.04, 1.04, 1.03, 1.01, 1.01, 0.99, 0.95, 0.94, 0.94, 0.92, 0.83, 0.83, 0.83, 0.83, 0.82, 0.82, 0.82, 0.81, 0.8, 0.79, 0.77, 0.77, 0.77, 0.75, 0.75, 0.74, 0.74, 0.74, 0.74, 0.74, 0.74]
for i in range(len(protein1)): #loops through a list of coev pairs and their scores
    coevr1 = IMP.pmi.restraints.basic.DistanceRestraint(representation=None,
                     tuple_selection1=(protein1[i],15,"P0AC02",0),
                     tuple_selection2=(protein2[i],25,"P0A937",0),
                     distancemin=0,
                     distancemax=21.2,
                     resolution=1.0,
                     kappa=1.0,
                     root_hier=root_hier,
                     weight=evscore[i])

    coevr1.add_to_model()
    output_objects.append(coevr1)
# -------------------------
# %%%%% EM RESTRAINT
#
# Scores a model based on its cross-correlation to an EM density.
# Since cross-correlation is very expensive, we approximate both
# the EM map and model as a set of 3D Gaussians (done in Representation).
#  
# First, collect all density particles from the model. 
em_components = IMP.pmi.tools.get_densities(root_hier)
gemt = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(em_components,gmm_data,scale_target_to_mass=True,slope=0.0000001,weight=80)
# gemt = IMP.pmi.restraints.em.GaussianEMRestraint(em_components,
#                                                 gmm_data,
#                                                 scale_target_to_mass=True,
#                                                 slope=0.0001,
#                                                 weight=100.0)

gemt.add_to_model()
output_objects.append(gemt)
sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(mdl))
print(sf.evaluate(True))
#####################################################
###################### SAMPLING #####################
#####################################################
# With our representation and scoring functions determined, we can now sample
# the configurations of our model with respect to the information.
# First shuffle all particles to randomize the starting point of the 
# system. For larger systems, you may want to increase max_translation
IMP.pmi.tools.shuffle_configuration(root_hier,max_translation=50)
# Shuffling randomizes the bead positions. It's good to 
# allow these to optimize first to relax large connectivity
# restraint scores.  100-500 steps is generally sufficient.
dof.optimize_flexible_beads(500)


# Now, add all of the other restraints to the scoring function to start sampling


num_frames = 10000
if '--test' in sys.argv: num_frames=100
num_mc_steps = 20

# This object defines all components to be sampled as well as the sampling protocol
mc1=IMP.pmi.macros.ReplicaExchange0(mdl,
                                    root_hier=root_hier,
                                    monte_carlo_sample_objects=dof.get_movers(),
                                    output_objects=output_objects,
                                    monte_carlo_temperature=1.0,
                                    simulated_annealing=True,
                                    simulated_annealing_minimum_temperature=1.0,
                                    simulated_annealing_maximum_temperature=5.0,
                                    simulated_annealing_minimum_temperature_nframes=200,
                                    simulated_annealing_maximum_temperature_nframes=20,
                                    replica_exchange_minimum_temperature=1.0,
                                    replica_exchange_maximum_temperature=2.5,
                                    number_of_best_scoring_models=10,
                                    monte_carlo_steps=num_mc_steps,
                                    number_of_frames=num_frames,
                                    global_output_directory="/home/cns-mccafferty/Ecoli_CoEv_modeling/CoEv_IMP/5d0o/15cutoff/run1/output/")

# Start Sampling
mc1.execute_macro()
