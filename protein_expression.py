import os, glob
import numpy as np
import pandas as pd

# Function to make dictionary mapping antibody names to their IgG control antibody
# Input
#   protein_control_df : A Pandas DataFrame with 'Antibody' column containing antibody names and 'Control' column containing correpsonding antibody control
# Output
#   ret_dict : A dictionary with antibody names as keys and their correspdoning IgG control antibody name as values
def make_control_dict(protein_control_df):
	ret_dict = {}
	for antibody in protein_control_df['Antibody'].to_numpy():
		ctrl_antibody = protein_control_df.loc[protein_control_df['Antibody']==antibody]['Control'].to_numpy()[0]
		ret_dict[antibody] = ctrl_antibody
	return ret_dict
