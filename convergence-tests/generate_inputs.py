import os

hbar = 1.0
m = 1.0
vd_np = 250

p0 = 10.0
n_sig = 4.0

beta = 4.0
alpha_p = 10.0
x0 = 0.0

vd_semi_classical = True

nx_min = 1000
d_nx = 495
n_nx = 200

nt_min = 1000
d_nt = 495
n_nt = 200

input_dir = "./convergence_inputs/"
output_dir = "./convergence_outputs/"

sig_p = p0 / alpha_p
sig_x = 0.5 * hbar / sig_p

vd_p_min = p0 - n_sig * sig_p
vd_p_max = p0 + n_sig * sig_p

x_min = x0 - n_sig * sig_x
x_max = x0 + n_sig * sig_x + beta * alpha_p / sig_p / m

t_min = 0.0
t_max = m / vd_p_min * (abs(x_min) + abs(x_max))

default_json = {
    "hbar": hbar,
    "m": m,
    "p0": p0,
    "x0": x0,
    "sig_x": sig_x,
    "nxl_external": 1,
    "nxr_external": 1,
    "nxl_vd": 1,
    "nxr_vd": 1,
    "resid_p_eps": 1.0E-6,
    "log_fname": "run.log",
    "psi_xt_fname": "psi_xt.dat",
    "vd_p_fname": "vd_px.dat",
    "vd_resid_fname": "vd_resid.dat",
    "vd_resid_analysis_fname": "vd_resid_analysis.txt",
    "vd_pt_fname": "vd_pxt.dat",
    "write_out_vd_pt": False,
    "write_out_psi_xt": False,
    "print_mod_x": 50,
    "print_mod_t": 10,
    "x_min": x_min,
    "x_max": x_max,
    "t_min": t_min,
    "t_max": t_max,
    "vd_p_min": vd_p_min,
    "vd_p_max": vd_p_max,
    "vd_np": vd_np
}


def ensure_dir(path):
    # Make sure a directory exists.
    os.makedirs(path, exist_ok=True)

ensure_dir(input_dir)

working_json = default_json.copy()

working_json["vd_semi_classical"] = vd_semi_classical

for i_x in range(n_nx):
    nx = nx_min + d_nx * i_x
    working_json["nx"] = nx
    for i_t in range(n_nt):
        nt = nt_min + d_nt * i_t

        working_json["nt"] = nt
        working_json["output_dir"] = output_dir+str(nx)+"_"+str(nt)+"/"

        inp_f = open(input_dir+str(nx)+"_"+str(nt)+".inp", "w")

        for (key, val) in working_json.items():
            inp_f.write(str(key) + " = " + str(val) + "\n")

        inp_f.close()
