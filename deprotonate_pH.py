import json
import os
import subprocess
import argparse

def run_ligprep(smiles_input, output_smi_file='schroout.smi'):
    # Define temporary input and output file paths
    input_smi_file = 'input.smi'

    # Write the input SMILES string to a temporary file
    with open(input_smi_file, 'w') as f:
        f.write(smiles_input + '\n')

    # Construct the LigPrep command
    ligprep_command = [
        os.path.join(os.getenv('SCHRODINGER'), 'ligprep'),
        '-s1',
        '-i', '2',
        '-W', 'i,-ph,7.4,-pht,0.0',
        '-t1',
        '-ismi', input_smi_file,
        '-osmi', output_smi_file,
        '-WAIT'
    ]

    try:
        # Run the LigPrep command
        subprocess.run(ligprep_command, check=True)

        # Read the output SMILES file and extract the SMILES string
        with open(output_smi_file, 'r') as f:
            result_smiles = f.readline().strip()

            # Remove anything after the first space, if present (to clean up the SMILES string)
            result_smiles = result_smiles.split()[0]

        # Return the SMILES string from the output file
        return result_smiles

    except subprocess.CalledProcessError as e:
        print(f"Error while running LigPrep: {e}")
        return None

    finally:
        # Clean up temporary input and output files
        if os.path.exists(input_smi_file):
            os.remove(input_smi_file)
        if os.path.exists(output_smi_file):
            os.remove(output_smi_file)


def process_json(file_path, output_path):
    """
    Parse a JSON file, update all ligand SMILES with protonation adjustments, and save the result.

    Args:
        file_path (str): Path to the input JSON file.
        output_path (str): Path to save the updated JSON file.
    """
    # Load the JSON file
    with open(file_path, "r") as f:
        data = json.load(f)
    
    # Iterate through sequences to find and process all ligand blocks
    for entry in data.get("sequences", []):
        if "ligand" in entry:
            ligand = entry["ligand"]
            smiles = ligand.get("smiles")

            if smiles:
                try:
                    # Update SMILES with protonation adjustment
                    updated_smiles = run_ligprep(smiles)

                    if updated_smiles:
                        ligand["smiles"] = updated_smiles
                        print(f"Updated SMILES for ligand {ligand['id']}: {updated_smiles}")
                except Exception as e:
                    print(f"Failed to update ligand {ligand['id']}: {e}")
    
    # Save the updated JSON
    with open(output_path, "w") as f:
        json.dump(data, f, indent=4)
    print(f"Updated JSON saved to {output_path}")



def process_directory(directory_path):
    """
    Process all JSON files in the specified directory, updating the ligand SMILES for each,
    and deleting the schroout.log file from the current working directory.
    
    Args:
        directory_path (str): Path to the directory containing JSON files.
    """
    # Check if the directory exists
    if not os.path.isdir(directory_path):
        print(f"Directory '{directory_path}' does not exist.")
        return
    
    for file_name in os.listdir(directory_path):
        if file_name.endswith(".json"):
            file_path = os.path.join(directory_path, file_name)
            output_path = os.path.join(directory_path, f"{file_name}")

            print(f"Processing {file_path}...")
            process_json(file_path, output_path)

    # Delete the schroout.log file from the current execution directory (not the provided directory)
    log_file_path = os.path.join(os.getcwd(), "schroout.log")
    if os.path.exists(log_file_path):
        os.remove(log_file_path)
        print(f"Deleted log file: {log_file_path}")

if __name__ == '__main__':
    # Set up argument parser to take directory input from the command line
    parser = argparse.ArgumentParser(description="Process JSON files and update ligand SMILES with protonation adjustments.")
    parser.add_argument('directory', type=str, help="Directory containing the JSON files to process.")
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Process the directory
    process_directory(args.directory)
