import os, math, requests, sys
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
    parser.add_argument('--top_n', type=int, default=None )
    parser.add_argument('--thre_score', type=float, default=None)

    return parser

_args = get_parser().parse_args()
output_dir = _args.output_path
vina_csv_path = _args.vina_csv_path
ligand_name = os.path.basename(vina_csv_path).replace('_vina_results.csv', '')

if _args.thre_score is not None and _args.top_n is not None:
    print("error: You cannot specify both --top_n and --thre_score")
    sys.exit()

os.makedirs(f'{output_dir}/string_result', exist_ok=True)
os.makedirs(f'{output_dir}/string_result/{ligand_name}', exist_ok=True)


#------------------------------------------------------------------------------------------------------
# Read vina_result_csv and extract uniprot list and ncbi id
df_vina_result = pd.read_csv(vina_csv_path, index_col=0)
AFs = df_vina_result.index.values.tolist()
UFs = []
uniprots = []
for AF in AFs:
    UFs.append(AF.replace('AF-','').replace('-model_v4',''))
for UF in UFs:
    uniprots.append(UF.split("-")[0])
uniprots = list(set(uniprots))

df_vina_result_new = pd.read_csv(vina_csv_path)
ncbi_id = df_vina_result_new.columns.values[0].split("_")[1]


#------------------------------------------------------------------------------------------------------
# Using the UniProt IDs, retrieve the corresponding STRING IDs through the STRING database API
# Save the results to a text file
string_api_url = "https://version-11-5.string-db.org/api"
output_format = "tsv-no-header"
method = "get_string_ids"

uniprots_split_list = [uniprots[i:i + 1000] for i in range(0, len(uniprots), 1000)]

output_lines = []
uniprot_string = {}

if not os.path.exists(f'{output_dir}/string_result/string_api_mapping.txt'):
    print("retrieving the STRING IDs through the STRING database API.")
    for uniprots_split in uniprots_split_list:
        params = {
            "identifiers": "\r".join(uniprots_split),  # Assuming 'UniprotID'
            "species": ncbi_id,  # Use the NCBI ID corresponding to the selected species
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

else:
    print("The String ID file is already exists.")
    with open(f'{output_dir}/string_result/string_api_mapping.txt') as file:
        lines = file.readlines()
        lines = [line.strip() for line in lines]
        for line in lines:
            l = line.split("\t")
            input_identifier, string_identifier = l[0][7:], l[1][8:]
            uniprot_string[input_identifier] = string_identifier


#------------------------------------------------------------------------------------------------------
# Extract TOP_n or under thre_score proteins
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

if _args.thre_score is None and _args.top_n is None:
    top_n = 100
    uniprot_score = sorted(uniprot_score.items(), key=lambda x:x[1]) # Sort by ascending score
    uniprot_score_extracted = uniprot_score[:top_n]
    print(f"Number of Extracted proteins : {len(uniprot_score_extracted)}")

if _args.top_n is not None:
    top_n = _args.top_n
    uniprot_score = sorted(uniprot_score.items(), key=lambda x:x[1]) # Sort by ascending score
    uniprot_score_extracted = uniprot_score[:top_n]
    print(f"Number of Extracted proteins : {len(uniprot_score_extracted)}")

if _args.thre_score is not None:
    thre_score = _args.thre_score
    uniprot_score_extracted = {uniprot: score for uniprot, score in uniprot_score.items() if score <= thre_score } # Extract under thre_score
    uniprot_score_extracted = sorted(uniprot_score_extracted.items(), key=lambda x:x[1]) # Sort by ascending score
    print(f"Number of Extracted proteins : {len(uniprot_score_extracted)}")

uniprot_string_extracted = {}
for uniprot_score in uniprot_score_extracted:
    uniprot = uniprot_score[0]
    if uniprot in uniprot_string.keys():
        uniprot_string_extracted[uniprot] = uniprot_string[uniprot]

identifiers_extracted = "%0d".join(uniprot_string_extracted.values())


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
    "identifiers": identifiers_extracted,
    "species": ncbi_id,
    "required_score": 400,
    "network_type": "functional"
}

# Fetch network image
network_response = requests.get(network_url, params=params_network)

if network_response.status_code == 200:
    output_image_path = f'{output_dir}/string_result/{ligand_name}/network_image_of_extracted_proteins.png'
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
        "species": ncbi_id, 
        "caller_identity": "www.awesome_app.org"  # Your app name
    }
    response = requests.post(request_url, data=params)
    if response.status_code == 200:
        return response.text
    else:
        print(f"Failed to retrieve data for the protein set, status code: {response.status_code}, response: {response.text}")
        return None
    
# Get and parse Functional annotation data
functional_annotation_data = get_functional_annotation(identifiers_extracted)
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
    grouped_results.to_csv(f'{output_dir}/string_result/{ligand_name}/functional_annotation_of_extracted_proteins.csv', index=False)
else:
    print("Failed to retrieve functional annotation data.")


#------------------------------------------------------------------------------------------------------
# Output the results in descending order of data quantity
# Load the CSV file
csv_file = f'{output_dir}/string_result/{ligand_name}/functional_annotation_of_extracted_proteins.csv'
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
output_image_path = f'{output_dir}/string_result/{ligand_name}/Number_of_proteins_for_Each_Function_in_extracted_Proteins.png'
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
        "species": ncbi_id, 
        "caller_identity": "www.awesome_app.org"  # Your app name
    }
    response = requests.post(request_url, data=params)
    if response.status_code == 200:
        return response.text
    else:
        print(f"Failed to retrieve data for the protein set, status code: {response.status_code}, response: {response.text}")
        return None

# Get and parse Functional Enrichments data
enrichment_data = get_functional_enrichments(identifiers_extracted)
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
    grouped_results = df[['term', 'description', 'number_of_genes', 'number_of_genes_in_background', 'fdr', 'genes_in_network']]
    grouped_results.columns = ['term', 'description', 'Count_in_Set', 'Count_in_Background', 'FDR', 'genes_in_network']

    # Save the results to a CSV file (optional)
    grouped_results.to_csv(f'{output_dir}/string_result/{ligand_name}/functional_enrichment_results_of_extracted_proteins.csv', index=False)
else:
    print("Failed to retrieve Functional Enrichments data.")


#------------------------------------------------------------------------------------------------------
# Output the results in descending order of data quantity
# The color of the bars is based on the FDR values
# Load the CSV file
csv_file = f'{output_dir}/string_result/{ligand_name}/functional_enrichment_results_of_extracted_proteins.csv'
df = pd.read_csv(csv_file)

# Sort the DataFrame by 'Count_in_Set' in descending order
df = df.sort_values(by='Count_in_Set', ascending=False)
df = df.head(30)

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
output_image_path = f'{output_dir}/string_result/{ligand_name}/functional_enrichment_results_of_extracted_Proteins_Colored_by_FDR.png'
plt.savefig(output_image_path)
plt.close()