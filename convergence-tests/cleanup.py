import os

input_dir = "./convergence_inputs/"
output_dir = "./convergence_outputs/"

print("Removing input files")
os.system("rm -rf " + input_dir)

print("Removing output directories / files")
os.system("rm -rf " + output_dir)
