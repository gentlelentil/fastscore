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
note the prepare_ligand4.py from MGL tools (/path/to/MGLxxx/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand.py) is bugged, the existing one in this directory is a template for the correct one (uncomment line 105, comment out 106)

added deprotonate function using schrodinger ligprep to generate AF3 structures at pH 7.4