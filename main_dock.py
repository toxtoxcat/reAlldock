import os, sys, subprocess
from argparse import ArgumentParser
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

os.environ['CUDA_VISIBLE_DEVICES'] = '1'

def get_parser():

    parser = ArgumentParser()
    parser.add_argument('--ligand_path', type=str, default=None)
    parser.add_argument('-s','--smiles', type=str, default=None)
    parser.add_argument('--ligand_name', type=str, default=None)
    parser.add_argument('-o','--output_path', type=str, default=None, required=True)
    parser.add_argument('--vina_executable_path', type=str, default=None, required=True)

    # Vina-GPU option
    parser.add_argument('--thread', type=str, default=1000)
    parser.add_argument('--size_x', type=str, default=20)
    parser.add_argument('--size_y', type=str, default=20)
    parser.add_argument('--size_z', type=str, default=20)

    return parser

_args = get_parser().parse_args()
input_ligand_path = _args.ligand_path
smiles = _args.smiles
if _args.ligand_name is not None:
    ligand_name = _args.ligand_name
elif smiles is not None:
    ligand_name = "ligand_from_smiles"
elif input_ligand_path is not None:
    ligand_name = os.path.basename(os.path.splitext(input_ligand_path)[0])

output_dir = _args.output_path
vina_executable = _args.vina_executable_path

size_x = _args.size_x
size_y = _args.size_y
size_z = _args.size_z
num_threads = _args.thread

if input_ligand_path is None and smiles is None:
    print("error: option -l/--ligand_path or -s/--smiles is required")
    sys.exit()

if input_ligand_path is not None and smiles is not None:
    print("error: You cannot specify both -l/--ligand_path and -s/--smiles")
    sys.exit()

os.makedirs(f'{output_dir}/dock_result/{ligand_name}', exist_ok=True)


#------------------------------------------------------------------------------------------------------
# convert smiles to sdf
if smiles is not None:
    mol = Chem.MolFromSmiles(smiles)
    hmol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(hmol,AllChem.ETKDG())
    writer = Chem.SDWriter(f'{output_dir}/dock_result/{ligand_name}/{ligand_name}.sdf')
    writer.write(hmol)
    writer.close()
    input_ligand_path = f'{output_dir}/dock_result/{ligand_name}/{ligand_name}.sdf'


#------------------------------------------------------------------------------------------------------
# convert ligand file to pdbqt
def convert_to_pdbqt(input_ligand_path, output_pdbqt_ligand):
    # Command to convert PDB to PDBQT with Open Babel
    command = f'obabel {input_ligand_path} -O {output_pdbqt_ligand} -h -p 7.4 --partialcharge gasteiger'

    # Run the command using subprocess
    subprocess.run(command, shell=True)

# Input and output file path
output_pdbqt_ligand = f'{output_dir}/dock_result/{ligand_name}/{ligand_name}.pdbqt'

# Call the function to convert and add charges
convert_to_pdbqt(input_ligand_path, output_pdbqt_ligand)


#------------------------------------------------------------------------------------------------------
# Run Vina-GPU
# Setting up necessary directories and parameters
config_folder = f'{output_dir}/prepare_result/config_docking'
output_folder = f'{output_dir}/dock_result/{ligand_name}/{ligand_name}_vina_results'
os.makedirs(output_folder, exist_ok=True)

# Navigate to the Vina-GPU-CUDA directory
vina_executable_dir = os.path.dirname(vina_executable)
os.chdir(f"{vina_executable_dir}")
print(os.getcwd())
# Add executable permission to the binary
os.system(f"chmod +x {vina_executable}")

# Get all config files in the config folder
config_files = [os.path.join(config_folder, f) for f in os.listdir(config_folder) if f.endswith('_config.txt')]

# Perform docking for each config file
for config_file in config_files:
    result_log = os.path.basename(config_file).replace("_config.txt","_out.txt")
    result_pdbqt = os.path.basename(config_file).replace("_config.txt","_out.pdbqt")
    if os.path.exists(f"{output_folder}/{result_pdbqt}") == False:
        # Construct the Vina-GPU execution command
        vina_cmd = f'{vina_executable} \
            --ligand {output_pdbqt_ligand} \
            --size_x {size_x} \
            --size_y {size_y} \
            --size_z {size_z} \
            --thread {num_threads} \
            --log {output_folder}/{result_log} \
            --out {output_folder}/{result_pdbqt} \
            --config {config_file}'

        # Execute the command
        os.system(vina_cmd)


#------------------------------------------------------------------------------------------------------
# Aggregate docking scores into a table

# Set the path to the folder containing receptor files
receptor_folder = f'{output_dir}/prepare_result/protein_rawpdb'

# Get the list of receptor files
receptor_names = [rn for rn in os.listdir(receptor_folder) if rn.endswith(".pdb")]

# Initialize a DataFrame
df = pd.DataFrame()

# Process docking results for each receptor
for rn in receptor_names:
    # Get the file name without extension from the receptor file name
    rn = rn[:-4]
    # Specify the path to the docking result file
    result_file = os.path.join(output_folder, f'{rn}_out.pdbqt')

    # If the docking result file exists
    if os.path.exists(result_file):
        with open(result_file) as file:
            lines = file.readlines()
            s_lines = [line.strip() for line in lines]
            s_lines_vina = [s_line for s_line in s_lines if 'REMARK VINA RESULT' in s_line]
            if len(s_lines_vina) != 0:
                score = float(s_lines_vina[0][23:29])
                df.loc[rn, "vina_score"] = score
            else:
                df.loc[rn, "vina_score"] = np.nan
    else:
        df.loc[rn, "vina_score"] = np.nan

# Creating a CSV file that consolidates docking scores.
# Set the path for the CSV file
csv_file = f'{output_dir}/dock_result/{ligand_name}/{ligand_name}_vina_results.csv'

#  Save the DataFrame as a CSV file
df.to_csv(csv_file)
fpocket_centroid_path = f'{output_dir}/prepare_result/fpocket-centroids_pLDDT.csv'
df_fpocket_plddt_new = pd.read_csv(fpocket_centroid_path)
ncbi_id = df_fpocket_plddt_new.columns.values[0].split("_")[1]
df_new = pd.read_csv(csv_file)
df_new.columns = [f"ncbi_{ncbi_id}_{ligand_name}"] + list(df.columns)
df_new.to_csv(csv_file, index=False)


#------------------------------------------------------------------------------------------------------
# Create a histogram of docking scores
# Display the histogram
df.hist(bins=40, color = "blue", grid =True)
plt.xlim(-15, 0)
plt.ylabel('Frequency')
plt.xlabel('Score (kcal/mol)')
plt.legend()
plt.title('Docking Scores Histogram')

# Save the image
output_image_path = f'{output_dir}/dock_result/{ligand_name}/{ligand_name}_docking_scores_histogram.png'
plt.savefig(output_image_path)
plt.close()