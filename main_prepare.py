from argparse import ArgumentParser
import os, sys, subprocess, tarfile, gzip, shutil, math
import numpy as np
import pandas as pd

def get_parser():

    parser = ArgumentParser()
    parser.add_argument('-o','--output_path', type=str, default=None, required=True)

    return parser

_args = get_parser().parse_args()
output_dir = _args.output_path

os.makedirs(f'{output_dir}/prepare_result', exist_ok=True)
output_dir = f'{output_dir}/prepare_result'

AF_ver = "v4"

#------------------------------------------------------------------------------------------------------
# Dictionary of species names and corresponding URLs
download_urls = {
    'Arabidopsis': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000006548_3702_ARATH_{AF_ver}.tar',
    'Nematode worm': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000001940_6239_CAEEL_{AF_ver}.tar',
    'C. albicans': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000000559_237561_CANAL{AF_ver}.tar',
    'Zebrafish': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000000437_7955_DANRE_{AF_ver}.tar',
    'Dictyostelium': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000002195_44689_DICDI_{AF_ver}.tar',
    'Fruit fly': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000000803_7227_DROME_{AF_ver}.tar',
    'E. coli': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000000625_83333_ECOLI_{AF_ver}.tar',
    'Soybean': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000008827_3847_SOYBN_{AF_ver}.tar',
    'Human': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000005640_9606_HUMAN_{AF_ver}.tar',
    'M. jannaschii': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000000805_243232_METJA_{AF_ver}.tar',
    'Mouse': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000000589_10090_MOUSE_{AF_ver}.tar',
    'Asian rice': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000059680_39947_ORYSJ_{AF_ver}.tar',
    'Rat': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000002494_10116_RAT_{AF_ver}.tar',
    'Budding yeast': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000002311_559292_YEAST_{AF_ver}.tar',
    'Fission yeast': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000002485_284812_SCHPO_{AF_ver}.tar',
    'Maize': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000007305_4577_MAIZE_{AF_ver}.tar',
    'Ajellomyces capsulatus': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000001631_447093_AJECG_{AF_ver}.tar',
    'Brugia malayi': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000006672_6279_BRUMA_{AF_ver}.tar',
    'C. jejuni': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000000799_192222_CAMJE_{AF_ver}.tar',
    'Cladophialophora carrionii': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000094526_86049_9EURO1_{AF_ver}.tar',
    'Dracunculus medinensis': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000274756_318479_DRAME_{AF_ver}.tar',
    'Enterococcus faecium': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000325664_1352_ENTFC_{AF_ver}.tar',
    'Fonsecaea pedrosoi': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000053029_1442368_9EURO2_{AF_ver}.tar',
    'H. influenzae': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000000579_71421_HAEIN_{AF_ver}.tar',
    'H. pylori': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000000429_85962_HELPY_{AF_ver}.tar',
    'K. pneumoniae': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000007841_1125630_KLEPH_{AF_ver}.tar',
    'L. infantum': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000008153_5671_LEIIN_{AF_ver}.tar',
    'Madurella mycetomatis': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000078237_100816_9PEZI1_{AF_ver}.tar',
    'Mycobacterium leprae': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000000806_272631_MYCLE_{AF_ver}.tar',
    'M. tuberculosis': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000001584_83332_MYCTU_{AF_ver}.tar',
    'Mycobacterium ulcerans': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000020681_1299332_MYCUL_{AF_ver}.tar',
    'N. gonorrhoeae': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000000535_242231_NEIG1_{AF_ver}.tar',
    'Nocardia brasiliensis': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000006304_1133849_9NOCA1_{AF_ver}.tar',
    'Onchocerca volvulus': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000024404_6282_ONCVO_{AF_ver}.tar',
    'Paracoccidioides lutzii': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000002059_502779_PARBA_{AF_ver}.tar',
    'P. falciparum': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000001450_36329_PLAF7_{AF_ver}.tar',
    'P. aeruginosa': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000002438_208964_PSEAE_{AF_ver}.tar',
    'S. typhimurium': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000001014_99287_SALTY_{AF_ver}.tar',
    'Schistosoma mansoni': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000008854_6183_SCHMA_{AF_ver}.tar',
    'S. dysenteriae': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000002716_300267_SHIDS_{AF_ver}.tar',
    'Sporothrix schenckii': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000018087_1391915_SPOS1_{AF_ver}.tar',
    'S. aureus': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000008816_93061_STAA8_{AF_ver}.tar',
    'S. pneumoniae': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000000586_171101_STRR6_{AF_ver}.tar',
    'Strongyloides stercoralis': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000035681_6248_STRER_{AF_ver}.tar',
    'Trichuris trichiura': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000030665_36087_TRITR_{AF_ver}.tar',
    'Trypanosoma brucei': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000008524_185431_TRYB2_{AF_ver}.tar',
    'T. cruzi': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000002296_353153_TRYCC_{AF_ver}.tar',
    'Wuchereria bancrofti': f'https://ftp.ebi.ac.uk/pub/databases/alphafold/{AF_ver}/UP000270924_6293_WUCBA_{AF_ver}.tar',
}

file_list = list(download_urls.keys())
print("available species:")
for index, filename in enumerate(file_list, start=1):
    print(f"{index}. {filename}")

try:
    user_choice = int(input("Please select the number of the species you download. (Enter 0 to exit.): "))
    if user_choice == 0:
        print("Exit program.")
        sys.exit()
    elif 1 <= user_choice <= len(file_list):
        selected_species = file_list[user_choice - 1]
        selected_path = download_urls[selected_species]
        os.system(f"wget {selected_path} -P {output_dir} -c")

    else:
        print("Invalid selection. Please select the correct number.")
        sys.exit()
except ValueError:
    print("Invalid selection. Please enter the number.")
    sys.exit()

## extract download data
def extract_tar_file(tar_path, extract_to):

    # Create the extraction destination folder (do nothing if it already exists)
    os.makedirs(extract_to, exist_ok=True)

    with tarfile.open(tar_path, "r:") as tar:
        tar.extractall(path=extract_to)
        print(f"File {tar_file} extracted to {extract_to}.")

tar_file = os.path.basename(selected_path)
tar_path = f'{output_dir}/{tar_file}'
extract_to = f'{output_dir}/protein_rawpdb'
extract_tar_file(tar_path, extract_to)

# Get all .gz files in the directory
gz_files = [f for f in os.listdir(extract_to) if f.endswith('.gz')]

# Extract each .gz file
for gz_file in gz_files:
    with gzip.open(f'{extract_to}/{gz_file}', 'rb') as f_in:
        with open(f'{extract_to}/{gz_file[:-3]}', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

# Calculate mean pLDDT
rawpdb_files = [f for f in os.listdir(f'{output_dir}/protein_rawpdb') if f.endswith('.pdb')]
rawpdb_files_plddt = {}

for rawpdb_file in rawpdb_files:
    with open(f'{output_dir}/protein_rawpdb/{rawpdb_file}') as file:
        lines = file.readlines()
        s_lines = [line.strip() for line in lines]
        s_lines_ATOM = [s_line for s_line in s_lines if s_line[:4]=="ATOM"]
        res_plddt = {}
        for s_line_ATOM in s_lines_ATOM:
            res_plddt[s_line_ATOM[22:26]] = float(s_line_ATOM[61:66])
        mean_plddt = round(sum(res_plddt.values())/len(res_plddt), 2)
        rawpdb_files_plddt[rawpdb_file] = mean_plddt


#------------------------------------------------------------------------------------------------------
# Create pdbqt files
# Add hydrogen and charge to proteins using OpenBabel
# Input folder containing PDB files
if not os.path.exists(f'{output_dir}/protein_pdbqt'):

    input_folder = f"{output_dir}/protein_rawpdb"  # Specify the path to the input folder containing PDB files

    # Output folder for converted PDBQT files
    output_folder = f"{output_dir}/added_pdbqt_folder"  # Specify the path to the output folder for converted PDBQT files

    # Create the output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Get a list of PDB files in the input folder
    pdb_files = [file for file in os.listdir(input_folder) if file.endswith('.pdb')]

    # Loop through each PDB file and convert to PDBQT format with added hydrogens and charges using Open Babel
    for pdb_file in pdb_files:
        input_pdb_path = os.path.join(input_folder, pdb_file)
        output_pdbqt_path = os.path.join(output_folder, os.path.splitext(pdb_file)[0] + '.pdbqt')

        # Command to add hydrogens, compute Gasteiger charges, and convert to PDBQT using Open Babel
        command = f'obabel {input_pdb_path} -xr --addpolarh --partialcharge gasteiger -O {output_pdbqt_path}'

        # Run the command using os.system
        try:
            if os.path.exists(output_pdbqt_path) == False :
                os.system(command)
                print(f'Converted file: {pdb_file}')
            elif os.path.getsize(output_pdbqt_path) == 0:
                os.system(command)
                print(f'Converted file: {pdb_file}')
        except Exception as e:
            print(f'Error converting file: {pdb_file}')
            print(e)

    # Remove everything except the "ATOM" lines from the PDBQT file and shorten the file name
    def extract_atom_lines_from_pdbqt(input_folder, output_folder):
        # Create the output folder if it doesn't exist
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # Get a list of PDBQT files in the input folder
        pdbqt_files = [file for file in os.listdir(input_folder) if file.endswith('.pdbqt')]

        # Loop through each PDBQT file and extract 'ATOM' lines
        for pdbqt_file in pdbqt_files:
            input_pdbqt_path = os.path.join(input_folder, pdbqt_file)
            output_pdbqt_path = os.path.join(output_folder, pdbqt_file)

            with open(input_pdbqt_path, 'r') as f_in:
                lines = f_in.readlines()

            atom_lines = [line for line in lines if line.startswith('ATOM')]

            with open(output_pdbqt_path, 'w') as f_out:
                f_out.writelines(atom_lines)

    # Input and output folder paths
    input_folder = f'{output_dir}/added_pdbqt_folder'  # Specify the path to the input folder containing PDBQT files
    output_folder = f'{output_dir}/atom_pdbpt_folder'  # Specify the path to the output folder for extracted PDBQT files

    # Call the function to extract 'ATOM' lines from PDBQT files
    extract_atom_lines_from_pdbqt(input_folder, output_folder)

    os.rename(f'{output_dir}/atom_pdbpt_folder',f'{output_dir}/protein_pdbqt')
    shutil.rmtree(f'{output_dir}/added_pdbqt_folder')


#------------------------------------------------------------------------------------------------------
# Run fpocket
os.makedirs(f'{output_dir}/fpocketout', exist_ok=True)

for i, pdb_file in enumerate(rawpdb_files):
    pdb_file_without_ext = os.path.splitext(os.path.basename(pdb_file))[0]
    if os.path.exists(f'{output_dir}/fpocketout/{pdb_file_without_ext}_out') == False:
        subprocess.run(['fpocket', '-f', f'{output_dir}/protein_rawpdb/{pdb_file}'])
        shutil.move(f'{output_dir}/protein_rawpdb/{pdb_file_without_ext}_out', f'{output_dir}/fpocketout')

# Calculate centroid of fpocket
def calculate_centroid_of_cluster1(pdb_file):
    """Calculate the centroid coordinates of cluster 1 from the _out.pdb file"""
    coordinates = []
    with open(pdb_file, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith("HETATM"):
                cluster_number = line[25]
                if cluster_number == "1":
                    x, y, z = float(line[30:37]), float(line[38:45]), float(line[46:53])
                    coordinates.append([x, y, z])

    if coordinates:
        centroid = np.mean(coordinates, axis=0)
        return centroid
    else:
        return None

directory = f'{output_dir}/fpocketout'
fpocket_results = {}

for folder in os.listdir(directory):
    folder_path = os.path.join(directory, folder)
    if os.path.isdir(folder_path):
        pdb_file = os.path.join(folder_path, f"{folder}.pdb")
        if os.path.isfile(pdb_file):
            centroid = calculate_centroid_of_cluster1(pdb_file)
            if centroid is not None:
                fpocket_results[folder] = centroid

# Write fpocket_centroids and pLDDT CSV
df_fpocket_plddt = pd.DataFrame(np.nan, index=rawpdb_files, columns=["center_x","center_y","center_z","pLDDT"])
for rawpdb_file in rawpdb_files:
    if f"{rawpdb_file[:-4]}_out" in fpocket_results.keys() :
        df_fpocket_plddt.loc[rawpdb_file, 'center_x'] = fpocket_results[f"{rawpdb_file[:-4]}_out"][0]
        df_fpocket_plddt.loc[rawpdb_file, 'center_y'] = fpocket_results[f"{rawpdb_file[:-4]}_out"][1]
        df_fpocket_plddt.loc[rawpdb_file, 'center_z'] = fpocket_results[f"{rawpdb_file[:-4]}_out"][2]
    else:
        df_fpocket_plddt.loc[rawpdb_file, 'center_x'] = np.nan
        df_fpocket_plddt.loc[rawpdb_file, 'center_y'] = np.nan
        df_fpocket_plddt.loc[rawpdb_file, 'center_z'] = np.nan
    df_fpocket_plddt.loc[rawpdb_file, 'pLDDT'] = rawpdb_files_plddt[rawpdb_file]
df_fpocket_plddt.to_csv(f'{output_dir}/fpocket-centroids_pLDDT.csv')


#------------------------------------------------------------------------------------------------------
# Create config files to run vina-GPU
# Define paths and parameters
receptor_folder = f'{output_dir}/protein_pdbqt'
config_folder = f'{output_dir}/config_docking'

# Create directories if they don't exist
os.makedirs(config_folder, exist_ok=True)

proteins = df_fpocket_plddt.index.values.tolist()

for protein in proteins:
    protein = protein[:-4]
    if math.isnan(df_fpocket_plddt.at[f"{protein}.pdb","center_x"]) == False:
        # Create the config file path
        config_file = os.path.join(config_folder, f'{protein}_config.txt')

        # Write to the config file
        with open(config_file, 'w') as config:
            config.write(f'receptor = {os.path.join(receptor_folder, protein + ".pdbqt")}\n')
            config.write(f'center_x = {df_fpocket_plddt.at[f"{protein}.pdb","center_x"]}\n')
            config.write(f'center_y = {df_fpocket_plddt.at[f"{protein}.pdb","center_y"]}\n')
            config.write(f'center_z = {df_fpocket_plddt.at[f"{protein}.pdb","center_z"]}\n')