Script to quickly calculate alphafold3 ligand docking scores.

Usage

$ python FastScore.py <input_dir> (-r)

-r: Operates recursively through the top level of a directory (eg directory with many AF3 outputs)

For single directories returns a csv with all ligands docked with a score for rigid and flexible docking using Vina conda implementation.
Requires open babel:

$sudo apt-get install obabel

and MGL tools:
https://ccsb.scripps.edu/mgltools/downloads/

you may need to reassert usr/bin/ to $PATH as MGL tools has a version of openBabel
