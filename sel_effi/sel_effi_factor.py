import json
from uncertainties import ufloat

# Function to load JSON data from a file
def load_json(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

# Function to save JSON data to a file
def save_json(data, file_path):
    with open(file_path, 'w') as file:
        json.dump(data, file, indent=4)


beam="pPb"
mode="ana"
selection="oscar_cuts"

# Load the JSON files
file1_data = load_json("/home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/sel_effi/sel_effi_data_"+beam+"_"+mode+"_"+selection+".json")
file2_data = load_json("/home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/sel_effi/sel_effi_mc_"+beam+"_"+mode+"_"+selection+".json")

# Prepare the result dictionary
result_data = {}

# Loop through the data and divide values and errors
for eta_bin, pt_bins in file1_data.items():
    result_data[eta_bin] = {}
    for pt_bin, values in pt_bins.items():
        value1 = file1_data[eta_bin][pt_bin]["sel_effi_data"]["value"]
        error1 = file1_data[eta_bin][pt_bin]["sel_effi_data"]["error"]
        value2 = file2_data[eta_bin][pt_bin]["sel_effi_mc"]["value"]
        error2 = file2_data[eta_bin][pt_bin]["sel_effi_mc"]["error"]

        # Create ufloat objects for division
        u1 = ufloat(value1, error1)
        u2 = ufloat(value2, error2)

        # Perform the division
        result = u1 / u2

        # Store the results
        result_data[eta_bin][pt_bin] = {
            "sel_effi_coeffi": {
                "value": result.nominal_value,
                "error": result.std_dev
            }
        }

# Save the result to a new JSON file
save_json(result_data, "/home3/sara.sellam/RpPb_identified_hadrons_project/efficiencies/sel_effi/sel_effi_coeffi_"+beam+"_"+mode+"_"+selection+".json")
