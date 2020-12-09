'''
This is Caitie's adaptation for 1L0V.
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
gmm_data = "/home/cns-mccafferty/Ecoli_CoEv_modeling/CoEv_IMP/1l0v/1l0v.50gmm.txt"

# Topology File
topology_file = "/home/cns-mccafferty/Ecoli_CoEv_modeling/CoEv_IMP/1l0v/topology.txt"
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Here is where the work begins
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# All IMP systems start out with a Model
mdl = IMP.Model()
# Read the topology file for a given state
t = IMP.pmi.topology.TopologyReader(topology_file,pdb_dir='/home/cns-mccafferty/Ecoli_CoEv_modeling/CoEv_IMP/1l0v/', fasta_dir='/home/cns-mccafferty/Ecoli_CoEv_modeling/CoEv_IMP/1l0v/', gmm_dir='/home/cns-mccafferty/Ecoli_CoEv_modeling/CoEv_IMP/1l0v/')
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
protein1 = [226, 163, 92, 203, 12, 100, 165, 226, 96, 138, 54, 161] #test data
protein2 = [94, 12, 10, 29, 10, 5, 32, 87, 23, 75, 116, 26]
evscore = [11.2, 10.2, 9.6, 8.8, 8.1, 7.8, 7.7, 7.6, 7.4, 7.4, 7.1, 6.7]
for i in range(len(protein1)): #loops through a list of coev pairs and their scores
    coevr1 = IMP.pmi.restraints.basic.DistanceRestraint(representation=None,
                     tuple_selection1=(protein1[i],15,"P0AC47",0),
                     tuple_selection2=(protein2[i],25,"P0A8Q0",0),
                     distancemin=0,
                     distancemax=29,
                     resolution=1.0,
                     kappa=1.0,
                     root_hier=root_hier,
                     weight=evscore[i])

    coevr1.add_to_model()
    output_objects.append(coevr1)

protein1 = [37, 116, 33, 7, 120, 77, 87, 37, 47, 34, 91, 33, 27, 121, 31, 5, 27, 50, 113, 117, 31, 31, 38, 124, 14, 113, 73, 75, 7, 64, 117, 68, 41, 91, 27, 31, 108, 10, 17, 54, 78, 57, 44, 47, 15, 37, 70, 88, 122, 15, 57, 54, 94, 116, 59, 62, 125, 27, 109, 43, 28, 23, 27, 13] #test data
protein2 = [79, 27, 83, 77, 27, 107, 23, 83, 57, 83, 23, 86, 86, 31, 85, 77, 90, 52, 23, 27, 86, 90, 73, 31, 66, 27, 37, 71, 13, 9, 84, 50, 110, 14, 85, 92, 23, 88, 9, 53, 4, 98, 25, 107, 86, 80, 104, 89, 62, 110, 50, 52, 16, 37, 110, 38, 31, 88, 79, 101, 37, 90, 31, 66]
evscore = [15.31, 13.85, 13.11, 11.39, 11.33, 10.92, 10.66, 10.6, 10.53, 10.28, 10.19, 10.15, 10.11, 10.04, 9.95, 9.76, 9.69, 9.46, 9.4, 9.28, 9.19, 9.17, 8.85, 8.82, 8.81, 8.78, 8.67, 8.59, 8.57, 8.31, 8.25, 8.22, 8.1, 8.01, 8.0, 7.88, 7.83, 7.76, 7.71, 7.7, 7.66, 7.61, 7.61, 7.58, 7.56, 7.47, 7.41, 7.36, 7.32, 7.31, 7.17, 7.16, 7.16, 7.15, 7.09, 7.06, 6.97, 6.96, 6.95, 6.92, 6.92, 6.91, 6.87, 6.87]
for i in range(len(protein1)): #loops through a list of coev pairs and their scores
    coevr2 = IMP.pmi.restraints.basic.DistanceRestraint(representation=None,
                     tuple_selection1=(protein1[i],15,"P0A8Q0",0),
                     tuple_selection2=(protein2[i],25,"P0A8Q3",0),
                     distancemin=0,
                     distancemax=29.7,
                     resolution=1.0,
                     kappa=1.0,
                     root_hier=root_hier,
                     weight=evscore[i])

    coevr2.add_to_model()
    output_objects.append(coevr2)

protein1 = [180, 53, 561, 510, 406, 424, 510] #test data
protein2 = [24, 50, 34, 107, 81, 81, 55]
evscore = [0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045]
for i in range(len(protein1)): #loops through a list of coev pairs and their scores
    coevr3 = IMP.pmi.restraints.basic.DistanceRestraint(representation=None,
                     tuple_selection1=(protein1[i],15,"P00363",0),
                     tuple_selection2=(protein2[i],25,"P0A8Q0",0),
                     distancemin=0,
                     distancemax=33.9,
                     resolution=1.0,
                     kappa=1.0,
                     root_hier=root_hier,
                     weight=evscore[i])

    coevr3.add_to_model()
    output_objects.append(coevr3)

protein1 = [432, 491, 424, 180, 180, 180, 180, 180, 180] #test data
protein2 = [39, 10, 8, 106, 115, 32, 60, 68, 83]
evscore = [0.079, 0.063, 0.063, 0.062, 0.062, 0.062, 0.062, 0.062, 0.062]
for i in range(len(protein1)): #loops through a list of coev pairs and their scores
    coevr4 = IMP.pmi.restraints.basic.DistanceRestraint(representation=None,
                     tuple_selection1=(protein1[i],15,"P00363",0),
                     tuple_selection2=(protein2[i],25,"P0A8Q3",0),
                     distancemin=0,
                     distancemax=27.4,
                     resolution=1.0,
                     kappa=1.0,
                     root_hier=root_hier,
                     weight=evscore[i])

    coevr4.add_to_model()
    output_objects.append(coevr4)

protein1 = [233] #test data
protein2 = [16]
evscore = [8.7]
for i in range(len(protein1)): #loops through a list of coev pairs and their scores
    coevr5 = IMP.pmi.restraints.basic.DistanceRestraint(representation=None,
                     tuple_selection1=(protein1[i],15,"P0AC47",0),
                     tuple_selection2=(protein2[i],25,"P0A8Q3",0),
                     distancemin=0,
                     distancemax=28.7,
                     resolution=1.0,
                     kappa=1.0,
                     root_hier=root_hier,
                     weight=evscore[i])

    coevr5.add_to_model()
    output_objects.append(coevr5)
# -------------------------
# %%%%% EM RESTRAINT
#
# Scores a model based on its cross-correlation to an EM density.
# Since cross-correlation is very expensive, we approximate both
# the EM map and model as a set of 3D Gaussians (done in Representation).
#  
# First, collect all density particles from the model. 
em_components = IMP.pmi.tools.get_densities(root_hier)
gemt = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(em_components,gmm_data,scale_target_to_mass=True,slope=0.0001,weight=100)
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
                                    global_output_directory="/home/cns-mccafferty/Ecoli_CoEv_modeling/CoEv_IMP/1l0v/20cutoff/run1/output/")

# Start Sampling
mc1.execute_macro()
