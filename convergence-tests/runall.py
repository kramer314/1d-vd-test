import os
import time
import multiprocessing

threads = 4

dev_null = "/dev/null"
input_dir = "./convergence_inputs/"

log_file = dev_null

call = "nice -n 19 ionice -c2 -n7 ../build/main.x "
call_end = " >> " + log_file

syscall_arr = []

input_files = os.listdir(input_dir)

if __name__ == "__main__":
   pool = multiprocessing.Pool(processes=threads)

   for fname in input_files:
      inp_path = input_dir + fname
      syscall = call + inp_path + call_end
      syscall_arr.append(syscall)

   if log_file is not dev_null:
      os.remove(log_file)

   start_time = time.time()

   pool.map(os.system, syscall_arr)

   pool.close()
   pool.join()

   end_time = time.time()

   print("Runtime: ", end_time-start_time)
