import pandas as pd
import os

# Directory containing the tables (CSV files in this example)
directory_path = "significant_hits_garr"

# Output file path for the extracted values
output_file_path = "ibd_garr_genes.txt"

# List to store values from the first column of each table
all_values = []

# Iterate through each file in the directory
for filename in os.listdir(directory_path):
    if filename.endswith(".tsv"):  # Adjust the file extension as needed
        file_path = os.path.join(directory_path, filename)

        # Read the table using pandas
        df = pd.read_csv(file_path,sep="\t")
        df = df[df['Functions'] == "GarD pathway/GarR"]

        # Extract values from the first column and append to the list
        first_column_values = df.iloc[:, 0].tolist()
        print(first_column_values)
        temp = {}
        if len(first_column_values) != 1:
                for i in first_column_values:
                        blast_file = "gudd_gard_garr/" + filename[:-3] + "blast_m6"
                        blast_df = pd.read_csv(blast_file,sep="\t",header=None)
                        target_rows = blast_df[blast_df.iloc[:, 0] == i]
                        print(target_rows)
                        target_row = target_rows.loc[target_rows.iloc[:, 11].idxmax()]
                        result_value = target_row.iloc[11]
                        temp[i] = result_value
                max_gene = max(temp, key=temp.get)
                all_values.append(max_gene.strip())
        else:
                all_values.extend(first_column_values)

# Write the values to the output file
with open(output_file_path, "w") as output_file:
    for value in all_values:
        output_file.write(str(value) + "\n")
