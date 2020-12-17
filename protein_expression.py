import os, glob
import numpy as np
import pandas as pd

# Function to make dictionary mapping antibody names to their IgG control antibody
# Input
# 	protein_control_df : A Pandas DataFrame with 'Antibody' column containing antibody names and 'Control' column containing correpsonding antibody control
# Output
# 	ret_dict : A dictionary with antibody names as keys and their correspdoning IgG control antibody name as values
def make_control_dict(protein_control_df):
	ret_dict = {}
	for antibody in protein_control_df['Antibody'].to_numpy():
		ctrl_antibody = protein_control_df.loc[protein_control_df['Antibody']==antibody]['Control'].to_numpy()[0]
		ret_dict[antibody] = ctrl_antibody
	return ret_dict

# Function to calculate normalized protein expression according to the formula max(0, log( (expr+1) / (ctrl_expr+1) ) )
# Input
# 	protein_df : A Pandas DataFrame with cell barcodes as columns and antibody UMI counts as rows
# 	protein_control_df : A Pandas DataFrame with 'Antibody' column containing antibody names and 'Control' column containing correpsonding antibody control
# 	barcodes_rna : Cell barcodes recovered from scRNA-seq data as NumPy array
# Output
# 	all_protein_expr : Normalized protein expression for each cell for each antibody
def calculate_protein_expression(protein_df, protein_control_df, barcodes_rna):
	
	control_map_dict = make_control_dict(protein_control_df)
	
	num_antibodies = protein_control_df['Antibody'].to_numpy().size
	expr_labels = protein_df['Antibody'].to_numpy()
	
	all_protein_expr = np.zeros((num_antibodies, barcodes_rna.size))
	for barcode_i, rna_barcode in enumerate(barcodes_rna):

		if rna_barcode in protein_df.columns: # Ensure high coverage of RNA barcodes

			# Antibody expression output by cumulus workflow
			curr_expr = protein_df[rna_barcode].to_numpy()

			assert(curr_expr.size == expr_labels.size)

			norm_expr = np.zeros(num_antibodies)
			for antibody_i, antibody in enumerate(antibody_names):
				targ_ind = np.where(expr_labels == antibody)[0]
				assert(targ_ind.size == 1)
				ctrl_ind = np.where(expr_labels == control_map_dict[antibody])[0]
				assert(ctrl_ind.size == 1)

				log_expr = np.log( (curr_expr[targ_ind]+1) / (curr_expr[ctrl_ind]+1) )
				norm_expr[antibody_i] = np.max([0, log_expr])

			all_protein_expr[:, barcode_i] = norm_expr

	return all_protein_expr
