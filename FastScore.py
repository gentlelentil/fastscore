import os
import subprocess
import sys
from vina import Vina
import csv
import argparse

def extract_ligands_and_protein(input_pdb, output_dir):
    """
    Extract ligands and protein from a PDB file.
    Ligands are saved separately, while protein atoms are saved in a protein-only PDB.
    """
    ligands = {}
    protein = []
    
    with open(input_pdb, 'r') as file:
        for line in file:
            if line.startswith("HETATM"):  # Ligands are HETATM lines
                ligand_chain = line[21:22]  # Chain ID for ligand
                if ligand_chain not in ligands:
                    ligands[ligand_chain] = []
                ligands[ligand_chain].append(line)
            elif line.startswith("ATOM"):  # Protein atoms are ATOM lines
                protein.append(line)

    # Save protein PDB
    protein_pdb = os.path.join(output_dir, "protein_no_ligand.pdb")
    with open(protein_pdb, 'w') as protein_file:
        protein_file.writelines(protein)
    print(f"Protein without ligands saved to {protein_pdb}")

    # Save each ligand as a separate PDB file
    for chain, atoms in ligands.items():
        ligand_pdb = os.path.join(output_dir, f"ligand_{chain}.pdb")
        with open(ligand_pdb, 'w') as ligand_file:
            ligand_file.writelines(atoms)
        print(f"Ligand {chain} extracted and saved to {ligand_pdb}")
    
    return protein_pdb, [os.path.join(output_dir, f"ligand_{chain}.pdb") for chain in ligands]

def fix_ligand_and_trim_pdb(input_pdb, output_pdb):
    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        lines = infile.readlines()

        # Find the last ATOM/HETATM line
        last_atom_index = max(
            i for i, line in enumerate(lines)
            if line.startswith("ATOM") or line.startswith("HETATM")
        )

        # Write up to the last ATOM/HETATM line (inclusive)
        for i, line in enumerate(lines):
            if i <= last_atom_index:
                # Convert ATOM to HETATM for ligands
                if line.startswith("ATOM") and "LIG" in line[17:20]:
                    line = line.replace("ATOM  ", "HETATM", 1)
                outfile.write(line)

        # Add the END line
        outfile.write("END\n")

def convert_cif_to_pdb(input_cif, output_pdb):
    command = [
        "/usr/bin/obabel",  # Assuming Open Babel is installed at this location
        input_cif,
        "-O",
        output_pdb,
        "-h"
    ]
    subprocess.run(command, check=True)
    print(f"Converted {input_cif} to {output_pdb}")

def prepare_receptor(receptor_pdb, output_pdbqt, mgltools_path):
    prepare_receptor_script = os.path.join(mgltools_path, "MGLToolsPckgs", "AutoDockTools", "Utilities24", "prepare_receptor4.py")
    command = f"pythonsh {prepare_receptor_script} -r {receptor_pdb} -o {output_pdbqt}"
    subprocess.run(command, shell=True, check=True)
    print(f"Receptor prepared: {output_pdbqt}")

def prepare_ligand(ligand_pdb, output_pdbqt, mgltools_path):
    prepare_ligand_script = os.path.join(mgltools_path, "MGLToolsPckgs", "AutoDockTools", "Utilities24", "prepare_ligand4.py")
    command = f"pythonsh {prepare_ligand_script} -l {ligand_pdb} -o {output_pdbqt}"
    subprocess.run(command, shell=True, check=True)
    print(f"Ligand prepared: {output_pdbqt}")

def calculate_centroid(pdb_file):
    x, y, z, count = 0, 0, 0, 0
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("HETATM"):
                coords = line[30:54].split()
                x += float(coords[0])
                y += float(coords[1])
                z += float(coords[2])
                count += 1
    if count == 0:
        raise ValueError(f"No atoms found in {pdb_file}")
    return [x / count, y / count, z / count]

def run_docking(receptor_pdbqt, ligand_pdbqt, center, output_dir, box_size=(20, 20, 20), ligand_name="ligand", protein_name="protein"):

    error = False

    # Vina docking
    box_size_str = f"--size_x {box_size[0]} --size_y {box_size[1]} --size_z {box_size[2]}"
    center_str = f"--center_x {center[0]} --center_y {center[1]} --center_z {center[2]}"
    
    # log_file = os.path.join(output_dir, "vina_output.log")
    output_pdbqt = os.path.join(output_dir, f"{ligand_name}_docked_out.pdbqt")

    v = Vina(sf_name='vina')
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)
    v.compute_vina_maps(center=center, box_size=box_size)


    try:
        # In-place binding (rigid)
        rigid_binding = v.score()

        # Localized docking (slight flexibility)
        v.dock(exhaustiveness=8, n_poses=1)

        # Retrieve docked poses
        poses = v.poses(n_poses=1)
        flex_energy = v.score()

        # Write the best docked pose
        v.write_poses(output_pdbqt, n_poses=1)

        # Prepare result for CSV
        result = [protein_name, ligand_name, rigid_binding[0], flex_energy[0]]

            # Append results to CSV
        output_csv = os.path.join(output_dir, "docking_results.csv")
        write_docking_results_to_csv([result], output_csv)

    except RuntimeError as e:
        error = True
        result = [protein_name, ligand_name, f"Error: {e}", ""]
        

        # Append results to CSV
        output_csv = os.path.join(output_dir, "docking_results.csv")
        write_docking_results_to_csv([result], output_csv)

    return result, error

def process_cif_file(cif_file, output_dir):
    pdb_file = os.path.join(output_dir, os.path.basename(cif_file).replace(".cif", "_bound.pdb"))
    convert_cif_to_pdb(cif_file, pdb_file)

    fixed_pdb_file = os.path.join(output_dir, os.path.basename(pdb_file).replace("_bound.pdb", "_fixed.pdb"))
    fix_ligand_and_trim_pdb(pdb_file, fixed_pdb_file)

    # Add hydrogens using Open Babel
    hydrogenated_pdb_file = os.path.join(output_dir, os.path.basename(fixed_pdb_file).replace("_fixed.pdb", "_withH.pdb"))
    subprocess.run(['obabel', fixed_pdb_file, '-O', hydrogenated_pdb_file, '-h'])

    return hydrogenated_pdb_file


def write_docking_results_to_csv(results, output_csv):
    """Write the docking results (protein, ligand, energies) to a CSV file."""
    with open(output_csv, mode='a', newline='') as file:
        writer = csv.writer(file)
        # Write header if the file is empty
        if file.tell() == 0:
            writer.writerow(['Protein', 'Ligand', 'Rigid Binding Energy (Before Minimization)', 'Flexible Binding Energy (After Minimization)'])
        
        # Write results
        for result in results:
            writer.writerow(result)




def fastdock(input_dir):
    output_dir = os.path.join(input_dir, "output")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    #add path to mgltools
    mgltools_path = ""
    if mgltools_path == "":
        print("Please provide the path to MGLTools")
        sys.exit(1)

    receptor_pdbqt = os.path.join(output_dir, "protein_no_ligand.pdbqt")

    cif_files = [
        os.path.join(input_dir, filename)
        for filename in os.listdir(input_dir) if filename.endswith(".cif")
    ]
    print(f"Found CIF files: {cif_files}")

    jobrun_results = []

    for cif_file in cif_files:

        protein_name = os.path.splitext(os.path.basename(cif_file))[0]

        ligand_receptor_pdb = process_cif_file(cif_file, output_dir)

        

        protein_pdb, ligand_pdbs = extract_ligands_and_protein(ligand_receptor_pdb, output_dir)
        
        prepare_receptor(protein_pdb, receptor_pdbqt, mgltools_path)

        for ligand_pdb in ligand_pdbs:

            ligand_name = os.path.splitext(os.path.basename(ligand_pdb))[0]

            ligand_pdbqt = os.path.join(output_dir, os.path.basename(ligand_pdb).replace(".pdb", ".pdbqt"))
            prepare_ligand(ligand_pdb, ligand_pdbqt, mgltools_path)

            centroid = calculate_centroid(ligand_pdb)
            # docking_score = run_docking(receptor_pdbqt, ligand_pdbqt, centroid, output_dir)
            results, error = run_docking(receptor_pdbqt, ligand_pdbqt, centroid, output_dir, ligand_name=ligand_name, protein_name=protein_name)
            # print(f"Docking score for {ligand_pdbqt}: {docking_score} kcal/mol")

            if error:
                print(results[-2])
                
            else:
                print('Score before minimization: %.3f (kcal/mol)' % results[-2])
                print('Score after minimization: %.3f (kcal/mol)' % results[-1])
            
            jobrun_results.append(results)
    
    return jobrun_results

        

def main(input_dir, recursive=False):
    # If recursive flag is set, find subdirectories
    if recursive:
        try:
            subdirectories = [os.path.join(input_dir, subdir) for subdir in os.listdir(input_dir)
                            if os.path.isdir(os.path.join(input_dir, subdir))]
            print(f"Found subdirectories: {subdirectories}")

            full_results = []
            
            # Execute fast_dock for each subdirectory
            for subdir in subdirectories:
                print(f"Running fast_dock for {subdir}")
                
                subdir_results = fastdock(subdir)
                full_results.extend(subdir_results)


            output_csv = os.path.join(input_dir, f"{input_dir.rstrip(os.sep)}_docking_results.csv")
            

            write_docking_results_to_csv(full_results, output_csv)
        
        except Exception as e:
            print(f"An error occurred: {e}")


    else:
        # Execute fast_dock normally for the provided directory
        print(f"Running fast_dock for {input_dir}")
        fastdock(input_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Docking script with optional recursion")
    parser.add_argument("input_dir", type=str, help="Input directory containing CIF files")
    parser.add_argument("--r", action="store_true", help="Process subdirectories recursively")
    args = parser.parse_args()

    main(args.input_dir, recursive=args.r)