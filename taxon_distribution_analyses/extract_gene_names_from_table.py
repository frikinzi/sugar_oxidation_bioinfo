import pandas as pd
import os

# Directory containing the tables (CSV files in this example)
directory_path = "significant_hits_garr"

# Output file path for the extracted values
output_file_path = "garl_genes.txt"

# List to store values from the first column of each table
all_values = []

# Iterate through each file in the directory
for filename in os.listdir(directory_path):
    if filename.endswith(".tsv"):  # Adjust the file extension as needed
        file_path = os.path.join(directory_path, filename)

        # Read the table using pandas
        df = pd.read_csv(file_path,sep="\t")
        most_common_value = df['Cluster_ID'].value_counts().idxmax()
        df = df[df['Cluster_ID'] == most_common_value]
        df = df[df['Functions'] == "GarD pathway/GarL"]

        # Extract values from the first column and append to the list
        first_column_values = df.iloc[:, 0].tolist()
        all_values.extend(first_column_values)

# Write the values to the output file
with open(output_file_path, "w") as output_file:
    for value in all_values:
        output_file.write(str(value) + "\n")
