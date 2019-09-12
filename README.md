# Pharmmaker

This is a suite of programs for making pharmacophore models from druggability simulations. See http://prody.csb.pitt.edu/tutorials/pharmmaker/ for more details.


Druggability simulations are molecular dynamics simulations of target proteins in solutions containing drug-like probe molecules. DruGUI is a toolkit for setting up and analyzing druggability simulations. See http://prody.csb.pitt.edu/drugui/ for more details.


Pharmmaker is to analyze the trajectory to construct pharmacophore models (PMs) to use for virtual screening of libraries of small molecules. It is using results of druggability simulations and analyses, and composed of a suite of steps as below:


1) Identify high affinity residues on the target protein for each probe molecule type from druggability simulations
2) Preselect high affinity residues near hot spots in a druggable site (hot spots are outputs of DruGUI analysis)
3) Analyze interactions between high affinity residues and probe at the hot spots for each probe type and Rank the interaction pairs between residue and probe.
4) Select snapshots with the top ranking interaction pairs from the druggability simulations
     

The selected snapshots have target conformations and poses of probes, and these are are used for the construction of Pharmacophore models, and the Pharmacophore models are then used as filters for identifying hits in structure-based virtual screening. 


A strong aspect of the method is that Pharmmaker uses multiple target conformations dependent on the binding poses of probes where they interact during druggability simulations. Therefore, the binding score in virtual screening can be more evaluated in a more realistic manner. Also, we can have multiple Pharmacophore models with different target conformations and probe poses, which can be analyzed statistically.
