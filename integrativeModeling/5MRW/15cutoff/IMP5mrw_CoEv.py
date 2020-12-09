'''
This is Caitie's adaptation for 5MRW.
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
gmm_data = "/home/cns-mccafferty/Ecoli_CoEv_modeling/CoEv_IMP/5mrw/5mrw.50gmm.txt"

# Topology File
topology_file = "/home/cns-mccafferty/Ecoli_CoEv_modeling/CoEv_IMP/5mrw/topology.txt"
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Here is where the work begins
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# All IMP systems start out with a Model
mdl = IMP.Model()
# Read the topology file for a given state
t = IMP.pmi.topology.TopologyReader(topology_file,pdb_dir='/home/cns-mccafferty/Ecoli_CoEv_modeling/CoEv_IMP/5mrw/', fasta_dir='/home/cns-mccafferty/Ecoli_CoEv_modeling/CoEv_IMP/5mrw/', gmm_dir='/home/cns-mccafferty/Ecoli_CoEv_modeling/CoEv_IMP/5mrw/')
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
protein1 = [426, 453, 418, 546, 537, 393, 542, 411, 526, 433, 556, 548, 521] #test data
protein2 = [234, 242, 227, 242, 584, 580, 235, 219, 651, 238, 604, 591, 647]
evscore = [1.85, 1.6, 1.57, 1.36, 1.36, 1.24, 1.06, 1.0, 0.88, 0.88, 0.83, 0.81, 0.79]
for i in range(len(protein1)): #loops through a list of coev pairs and their scores
    coevr1 = IMP.pmi.restraints.basic.DistanceRestraint(representation=None,
                     tuple_selection1=(protein1[i],15,"P03959",0),
                     tuple_selection2=(protein2[i],25,"P03960",0),
                     distancemin=0,
                     distancemax=18.3,
                     resolution=1.0,
                     kappa=1.0,
                     root_hier=root_hier,
                     weight=evscore[i])

    coevr1.add_to_model()
    output_objects.append(coevr1)

protein1 = [186, 183, 448, 204, 249, 11, 249, 171, 191, 183, 187, 198, 325, 76, 183, 241, 325, 46, 244, 174, 325, 198, 50, 446, 190, 308, 72, 207, 50, 353, 199, 447, 248, 214, 130, 131, 309, 250, 46, 215, 49, 202, 139, 119, 320, 204, 8, 171, 208, 354, 325, 196, 308, 238, 208, 215, 207, 325, 242, 448, 191, 352, 170, 49, 202, 455, 214, 122, 356, 72, 205, 76, 472, 312, 50, 554, 250, 139, 250, 170, 319, 174, 314, 122, 309, 446, 215, 187, 447, 4, 208, 349, 355, 135, 190, 314, 135, 72, 130, 237, 175, 179, 167, 204, 213, 553, 315, 352, 206, 312, 201, 253, 93, 8, 537, 221, 134, 310, 244, 246, 350, 72, 135, 306, 443, 142] #test data
protein2 = [26, 26, 91, 49, 175, 170, 167, 14, 34, 25, 29, 24, 175, 11, 29, 125, 119, 13, 55, 18, 131, 27, 10, 87, 26, 95, 11, 179, 9, 85, 39, 91, 167, 55, 23, 15, 98, 167, 17, 80, 13, 49, 11, 81, 118, 42, 166, 10, 60, 125, 134, 27, 90, 82, 59, 79, 183, 133, 128, 92, 37, 88, 14, 9, 42, 88, 42, 80, 90, 8, 48, 15, 125, 95, 13, 90, 171, 10, 168, 10, 118, 14, 102, 81, 95, 93, 82, 27, 88, 166, 64, 125, 89, 14, 38, 121, 15, 4, 24, 82, 17, 18, 7, 50, 80, 90, 118, 91, 58, 99, 80, 171, 131, 165, 92, 26, 23, 89, 53, 125, 125, 7, 11, 94, 93, 18]
evscore = [6.09, 5.97, 4.89, 4.57, 4.34, 3.68, 3.18, 3.1, 3.08, 3.01, 2.98, 2.78, 2.76, 2.73, 2.61, 2.59, 2.42, 2.37, 2.36, 2.31, 2.26, 2.26, 2.14, 2.1, 2.09, 2.04, 2.03, 2.01, 2.01, 1.91, 1.9, 1.88, 1.85, 1.83, 1.78, 1.75, 1.73, 1.72, 1.69, 1.63, 1.63, 1.63, 1.62, 1.6, 1.52, 1.52, 1.51, 1.48, 1.46, 1.45, 1.42, 1.4, 1.4, 1.38, 1.36, 1.34, 1.31, 1.3, 1.3, 1.29, 1.28, 1.28, 1.28, 1.27, 1.27, 1.27, 1.26, 1.25, 1.24, 1.23, 1.23, 1.23, 1.22, 1.21, 1.15, 1.15, 1.14, 1.13, 1.12, 1.09, 1.08, 1.08, 1.07, 1.06, 1.06, 1.05, 1.04, 1.04, 1.04, 1.04, 1.04, 1.03, 1.03, 1.03, 1.01, 1.0, 1.0, 0.99, 0.99, 0.99, 0.98, 0.97, 0.96, 0.96, 0.94, 0.94, 0.93, 0.93, 0.91, 0.91, 0.9, 0.9, 0.89, 0.89, 0.88, 0.87, 0.86, 0.86, 0.86, 0.86, 0.83, 0.82, 0.81, 0.81, 0.81, 0.8]
for i in range(len(protein1)): #loops through a list of coev pairs and their scores
    coevr2 = IMP.pmi.restraints.basic.DistanceRestraint(representation=None,
                     tuple_selection1=(protein1[i],15,"P03959",0),
                     tuple_selection2=(protein2[i],25,"P03961",0),
                     distancemin=0,
                     distancemax=19.2,
                     resolution=1.0,
                     kappa=1.0,
                     root_hier=root_hier,
                     weight=evscore[i])

    coevr2.add_to_model()
    output_objects.append(coevr2)

protein1 = [609] #test data
protein2 = [131]
evscore = [1.38]
for i in range(len(protein1)): #loops through a list of coev pairs and their scores
    coevr3 = IMP.pmi.restraints.basic.DistanceRestraint(representation=None,
                     tuple_selection1=(protein1[i],15,"P03960",0),
                     tuple_selection2=(protein2[i],25,"P03961",0),
                     distancemin=0,
                     distancemax=19.3,
                     resolution=1.0,
                     kappa=1.0,
                     root_hier=root_hier,
                     weight=evscore[i])

    coevr3.add_to_model()
    output_objects.append(coevr3)

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
                                    global_output_directory="/home/cns-mccafferty/Ecoli_CoEv_modeling/CoEv_IMP/5mrw/15cutoff/run1/output/")

# Start Sampling
mc1.execute_macro()
