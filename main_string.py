import os, math, requests
from argparse import ArgumentParser
from time import sleep
from PIL import Image
import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO
import seaborn as sns

def get_parser():

    parser = ArgumentParser()
    parser.add_argument('-o','--output_path', type=str, default=None, required=True)
    parser.add_argument('--vina_csv_path', type=str, default=None, required=True)
    parser.add_argument('--top_n', type=float, default = 100 )
    parser.add_argument('--thre_score', type=float, default = -8 )

    return parser

_args = get_parser().parse_args()
output_dir = _args.output_path
vina_csv_path = _args.vina_csv_path
ligand_name = os.path.basename(vina_csv_path).replace('_vina_results.csv', '')
top_n = _args.top_n
thre_score = _args.thre_score

os.makedirs(f'{output_dir}/string_result', exist_ok=True)
os.makedirs(f'{output_dir}/string_result/{ligand_name}', exist_ok=True)


#------------------------------------------------------------------------------------------------------
# Read vina_result_csv and extract uniprot list
df_vina_result = pd.read_csv(vina_csv_path, index_col=0)
AFs = df_vina_result.index.values.tolist()
UFs = []
uniprots = []
for AF in AFs:
    UFs.append(AF.replace('AF-','').replace('-model_v4',''))
for UF in UFs:
    uniprots.append(UF.split("-")[0])
uniprots = list(set(uniprots))


#------------------------------------------------------------------------------------------------------
# Dictionary mapping species names to their NCBI IDs
species_ncbi_ids = {
    'Arabidopsis': 3702,
    'Nematode worm': 6239,
    'C. albicans': 237561,
    'Zebrafish': 7955,
    'Dictyostelium': 44689,
    'Fruit fly': 7227,
    'E. coli': 83333,
    'Soybean': 3847,
    'Human': 9606,
    'M. jannaschii': 243232,
    'Mouse': 10090,
    'Asian rice': 39947,
    'Rat': 10116,
    'Budding yeast': 559292,
    'Fission yeast': 284812,
    'Maize': 4577,
    'Ajellomyces capsulatus': 447093,
    'Brugia malayi': 6279,
    'C. jejuni': 192222,
    'Cladophialophora carrionii': 86049,
    'Dracunculus medinensis': 318479,
    'Enterococcus faecium': 1352,
    'Fonsecaea pedrosoi': 1442368,
    'H. influenzae': 71421,
    'H. pylori': 85962,
    'K. pneumoniae': 1125630,
    'L. infantum': 5671,
    'Madurella mycetomatis': 100816,
    'Mycobacterium leprae': 272631,
    'M. tuberculosis': 83332,
    'Mycobacterium ulcerans': 1299332,
    'N. gonorrhoeae': 242231,
    'Nocardia brasiliensis': 1133849,
    'Onchocerca volvulus': 6282,
    'Paracoccidioides lutzii': 502779,
    'P. falciparum': 36329,
    'P. aeruginosa': 208964,
    'S. typhimurium': 99287,
    'Schistosoma mansoni': 6183,
    'S. dysenteriae': 300267,
    'Sporothrix schenckii': 1391915,
    'S. aureus': 93061,
    'S. pneumoniae': 171101,
    'Strongyloides stercoralis': 6248,
    'Trichuris trichiura': 36087,
    'Trypanosoma brucei': 185431,
    'T. cruzi': 353153,
    'Wuchereria bancrofti': 6293,
}

# Prompt user to select a species
print("Select a species:")
for index, species in enumerate(species_ncbi_ids.keys(), 1):
    print(f"{index}. {species}")

selection = input("Enter the number corresponding to the desired species: ")

# Validate user input
try:
    selection_index = int(selection)
    if 1 <= selection_index <= len(species_ncbi_ids):
        selected_species = list(species_ncbi_ids.keys())[selection_index - 1]
        ncbi_id = species_ncbi_ids[selected_species]
    else:
        raise ValueError
except ValueError:
    print("Invalid input. Please enter a valid number.")
    exit()

print(f"NCBI ID for {selected_species}: {ncbi_id}")


#------------------------------------------------------------------------------------------------------
# Using the UniProt IDs, retrieve the corresponding STRING IDs through the STRING database API
# Save the results to a text file
string_api_url = "https://version-11-5.string-db.org/api"
output_format = "tsv-no-header"
method = "get_string_ids"

uniprots_split_list = [uniprots[i:i + 2000] for i in range(0, len(uniprots), 2000)]

output_lines = []
uniprot_string = {}

for uniprots_split in uniprots_split_list:
    params = {
        "identifiers": "\r".join(uniprots_split),  # Assuming 'UniprotID'
        "species": species_ncbi_ids[selected_species],  # Use the NCBI ID corresponding to the selected species
        "limit": 1,  # only one (best) identifier per input protein
        "echo_query": 1,  # see your input identifiers in the output
        "caller_identity": "www.awesome_app.org"  # your app name
    }

    ##
    ## Construct URL
    ##
    request_url = "/".join([string_api_url, output_format, method])

    ##
    ## Call STRING
    ##
    results = requests.post(request_url, data=params)

    ##
    ## Read and parse the results
    ##
    for line in results.text.strip().split("\n"):
        l = line.split("\t")
        input_identifier, string_identifier = l[0], l[2]
        uniprot_string[input_identifier] = string_identifier
        output_lines.append(f"Input: {input_identifier}\tSTRING: {string_identifier}")

    sleep(1)

##
## Write results to text file
##
string_id_file_path = f'{output_dir}/string_result/string_api_mapping.txt'
with open(string_id_file_path, "w") as output_file:
    output_file.write("\n".join(output_lines))

print(f"Mapping results saved to {string_id_file_path}")


#------------------------------------------------------------------------------------------------------
# Read STRING IDs from a text file, generate a network image, and compress the output
output_format = "image"
method = "network"

# Set the path to the output folder
output_folder = f"{output_dir}/string_result/network_results"

# Create the output folder if it does not exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

##
## Construct URL
##
request_url = "/".join([string_api_url, output_format, method])

## For each protein call STRING
for uniprot in uniprot_string.keys():
    string = uniprot_string[uniprot]
    file_name = os.path.join(output_folder, f"{uniprot}_{string}_network.png")

    if os.path.exists(file_name) == False:
        ##
        ## Set parameters
        ##
        params = {
            "identifiers": string,  # your protein
            "species": species_ncbi_ids[selected_species],  # species NCBI identifier
            "add_white_nodes": 15,  # add 15 white nodes to my protein
            "network_flavor": "confidence",  # show confidence links
            "caller_identity": "www.awesome_app.org"  # your app name
        }

        ##
        ## Call STRING
        ##
        response = requests.post(request_url, data=params)

        if response.status_code == 200:
            ##
            ## Save the network to file
            ##
            print("Saving interaction network to %s" % file_name)

            with open(file_name, 'wb') as fh:
                fh.write(response.content)

            img = Image.open(file_name)
            new_img = Image.new("RGB", img.size, "WHITE")
            new_img.paste(img, (0, 0), img)
            new_img.save(file_name)

        sleep(2)


#------------------------------------------------------------------------------------------------------
# Extract TOP_n and under thre_score proteins
uniprot_score = {}

for uniprot in uniprots:
    s_UFs = [UF for UF in UFs if UF.startswith(f"{uniprot}-")]
    min_score_in_s_UFs = 0
    for s_UF in s_UFs:
        if math.isnan(df_vina_result.at[f"AF-{s_UF}-model_v4","vina_score"]) == False:
            s_UF_score = float(df_vina_result.at[f"AF-{s_UF}-model_v4","vina_score"])
        else:
            s_UF_score = 0
        if s_UF_score < min_score_in_s_UFs:
            min_score_in_s_UFs = s_UF_score
    uniprot_score[uniprot] = min_score_in_s_UFs

uniprot_score = {uniprot: score for uniprot, score in uniprot_score.items() if score <= thre_score } # Extract under thre_score
uniprot_score = sorted(uniprot_score.items(), key=lambda x:x[1]) # Sort by ascending score
uniprot_score_top_n = uniprot_score[:top_n]

uniprot_string_top_n = {}
for uniprot_score in uniprot_score_top_n:
    uniprot = uniprot_score[0]
    if uniprot in uniprot_string.keys():
        uniprot_string_top_n[uniprot] = uniprot_string[uniprot]

identifiers_top_n = "%0d".join(uniprot_string_top_n.values())


#------------------------------------------------------------------------------------------------------
# Output the network image of the high docking score proteins
# Set up STRING API
string_api_url = "https://string-db.org/api"
output_format = "image"
method_network = "network"

# API endpoint
network_url = f"{string_api_url}/{output_format}/{method_network}"

# Parameters for STRING API
params_network = {
    "identifiers": identifiers_top_n,
    "species": species_ncbi_ids[selected_species],
    "required_score": 400,
    "network_type": "functional"
}

# Fetch network image
network_response = requests.get(network_url, params=params_network)

if network_response.status_code == 200:
    output_image_path = f'{output_dir}/string_result/{ligand_name}/network_image_of_TOP{top_n}_proteins_for_{ligand_name}.png'
    with open(output_image_path, 'wb') as fh:
        fh.write(network_response.content)
        fh.close
        img = Image.open(output_image_path)
        new_img = Image.new("RGB", img.size, "WHITE")
        new_img.paste(img, (0, 0), img)
        new_img.save(output_image_path)
else:
    print("Network retrieval failed:", network_response.status_code, network_response.text)


#------------------------------------------------------------------------------------------------------
# Retrieve functional annotation such as GO_term of the top n docking score proteins
# # Set up your STRING API credentials
string_api_url = "https://version-11-5.string-db.org/api"
output_format = "tsv-no-header"
method = "functional_annotation"

# Function to get Functional Enrichments data for the entire protein set
def get_functional_annotation(identifiers):
    request_url = f"{string_api_url}/{output_format}/{method}"
    params = {
        "identifiers": identifiers,
        "species": species_ncbi_ids[selected_species], 
        "caller_identity": "www.awesome_app.org"  # Your app name
    }
    response = requests.post(request_url, data=params)
    if response.status_code == 200:
        return response.text
    else:
        print(f"Failed to retrieve data for the protein set, status code: {response.status_code}, response: {response.text}")
        return None
    
# Get and parse Functional annotation data
functional_annotation_data = get_functional_annotation(identifiers_top_n)
if functional_annotation_data:
    df = pd.read_csv(StringIO(functional_annotation_data), sep='\t', header=None)
    # Manually inspect and set the appropriate column names
    column_names = [
        'category', 'term', 'number_of_genes', 'ratio_in_set',
        'ncbiTaxonId', 'inputGenes', 'preferredNames', 'description'
    ]

    if len(df.columns) >= len(column_names):
        df.columns = column_names + [f"extra_{i}" for i in range(len(df.columns) - len(column_names))]

    # Select relevant columns including description
    grouped_results = df[['category', 'term', 'description', 'number_of_genes', 'inputGenes']]
    grouped_results.columns = ['category', 'term', 'description', 'Count_in_Set', 'inputGenes']

    # Save the results to a CSV file (optional)
    grouped_results.to_csv(f'{output_dir}/string_result/{ligand_name}/functional_annotation_of_TOP{top_n}_proteins_for_{ligand_name}.csv', index=False)
else:
    print("Failed to retrieve functional annotation data.")


#------------------------------------------------------------------------------------------------------
# Output the results in descending order of data quantity
# Load the CSV file
csv_file = f'{output_dir}/string_result/{ligand_name}/functional_annotation_of_TOP{top_n}_proteins_for_{ligand_name}.csv'
df = pd.read_csv(csv_file)

# Sort the DataFrame by 'Count_in_Set' in descending order
df = df.sort_values(by='Count_in_Set', ascending=False)
df = df.head(30)
df = df.sort_values(['category','Count_in_Set'], ascending=[False, False])

# Draw bar graph
plt.figure(figsize=(25, 20))
plt.rcParams["font.size"] = 20
barplot = sns.barplot(x='Count_in_Set', y='description', hue='category',palette="viridis", data=df, dodge=False)

plt.xlabel('Number of Proteins')
plt.ylabel('Function')
plt.title('Number of Proteins for Each Function')
plt.legend(title='category')

plt.tight_layout()

# Save the image
output_image_path = f'{output_dir}/string_result/{ligand_name}/Number_of_proteins_for_Each_Function_in_TOP{top_n}_Proteins_for_{ligand_name}.png'
plt.savefig(output_image_path)
plt.close()


#------------------------------------------------------------------------------------------------------
# Perform enrichment analysis on the top n docking score proteins
# Obtain 'term', 'description', 'Count_in_Set', 'Count_in_Background', and 'FDR'
# # Set up your STRING API credentials
string_api_url = "https://version-11-5.string-db.org/api"
output_format = "tsv-no-header"
method = "enrichment"

# Function to get Functional Enrichments data for the entire protein set
def get_functional_enrichments(identifiers):
    request_url = f"{string_api_url}/{output_format}/{method}"
    params = {
        "identifiers": identifiers,
        "species": species_ncbi_ids[selected_species], 
        "caller_identity": "www.awesome_app.org"  # Your app name
    }
    response = requests.post(request_url, data=params)
    if response.status_code == 200:
        return response.text
    else:
        print(f"Failed to retrieve data for the protein set, status code: {response.status_code}, response: {response.text}")
        return None

# Get and parse Functional Enrichments data
enrichment_data = get_functional_enrichments(identifiers_top_n)
if enrichment_data:
    df = pd.read_csv(StringIO(enrichment_data), sep='\t', header=None)
    # Manually inspect and set the appropriate column names
    column_names = [
        'category', 'term', 'number_of_genes', 'number_of_genes_in_background',
        'species', 'genes_in_network', 'genes_in_background', 'p_value', 'fdr', 'description'
    ]
    if len(df.columns) >= len(column_names):
        df.columns = column_names + [f"extra_{i}" for i in range(len(df.columns) - len(column_names))]

    # Select relevant columns including description
    grouped_results = df[['term', 'description', 'number_of_genes', 'number_of_genes_in_background', 'fdr']]
    grouped_results.columns = ['term', 'description', 'Count_in_Set', 'Count_in_Background', 'FDR']

    # Save the results to a CSV file (optional)
    grouped_results.to_csv(f'{output_dir}/string_result/{ligand_name}/functional_enrichment_results_of_TOP{top_n}_proteins_for_{ligand_name}.csv', index=False)
else:
    print("Failed to retrieve Functional Enrichments data.")


#------------------------------------------------------------------------------------------------------
# Output the results in descending order of data quantity
# The color of the bars is based on the FDR values
# Load the CSV file
csv_file = f'{output_dir}/string_result/{ligand_name}/functional_enrichment_results_of_TOP{top_n}_proteins_for_{ligand_name}.csv'
df = pd.read_csv(csv_file)

# Sort the DataFrame by 'Count_in_Set' in descending order
df = df.sort_values(by='FDR', ascending=True)
df = df.head(30)
df = df.sort_values(['FDR','Count_in_Set'], ascending=[True, False])

# Create a color palette based on the FDR values
norm = plt.Normalize(df['FDR'].min(), df['FDR'].max())
sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
sm.set_array([])

# Draw bar graph
plt.figure(figsize=(30, 20))
barplot = sns.barplot(x='Count_in_Set', y='description', data=df, palette="viridis", hue='FDR', dodge=False, legend=False)

# Add color bar
cbar = plt.colorbar(sm, ax=barplot, orientation='vertical')
cbar.set_label('FDR')

plt.xlabel('Number of Proteins')
plt.ylabel('Function')
plt.title('Number of Proteins for Each Function (Colored by FDR)')

plt.tight_layout()

# Save the image
output_image_path = f'{output_dir}/string_result/{ligand_name}/functional_enrichment_results_of_TOP{top_n}_Proteins_for_{ligand_name}_Colored_by_FDR.png'
plt.savefig(output_image_path)
plt.close()