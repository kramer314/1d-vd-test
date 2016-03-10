import os
import numpy as np

nx_min = 1000
d_nx = 495
n_nx = 200

nt_min = 1000
d_nt = 495
n_nt = 200

output_dir = "./convergence_outputs/"
results_dir = "./convergence_results/"

resid_fname = "vd_resid_analysis.txt"

resid_matrix = np.zeros((n_nx,n_nt))

def ensure_dir(path):
    # Make sure a directory exists.
    os.makedirs(path, exist_ok=True)

ensure_dir(results_dir)

for i_x in range(n_nx):
    nx = nx_min + d_nx * i_x

    for i_t in range(n_nt):
        nt = nt_min + d_nt * i_t

        dir_name = output_dir+str(nx)+"_"+str(nt)+"/"

        f_name = dir_name + resid_fname
        f = open(f_name, "r")

        lines = f.readlines()

        f.close()

        for line in lines:
            if (line.find("Cumulative") >= 0):
                split = line.rstrip("\n").split(":")
                resid_matrix[i_x][i_t] = float(split[1])

resid_matrix_t = np.transpose(resid_matrix)

f = open(results_dir+"rp10-b4-cxm_tx.dat", "w")
for i_x in range(n_nx):
    out_line = []

    for i_t in range(n_nt):
        out_line.append(resid_matrix[i_x][i_t])

    out = "\t".join(str(i) for i in out_line)+"\n"
    f.write(out)
f.close()

g = open(results_dir+"rp10-b4-cm_xt.dat", "w")
for i_t in range(n_nt):
    out_line = []

    for i_x in range(n_nx):
        out_line.append(resid_matrix_t[i_t][i_x])
    out = "\t".join(str(i) for i in out_line)+"\n"
    g.write(out)
g.close()
