import subprocess
import os
import shutil
import filecmp
import csv

def compare_csv_matrices(file1, file2, tolerance=1e-5):
    mismatches = 0
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        reader1 = csv.reader(f1)
        reader2 = csv.reader(f2)
        cmb_fl_ter = zip(reader1, reader2)
        for row_idx, (row1, row2) in enumerate(cmb_fl_ter):
            cmb_row_iter = zip(row1, row2)
            for col_idx, (val1, val2) in enumerate(cmb_row_iter):
                try:
                    num1 = float(val1)
                    num2 = float(val2)
                    if abs(num1 - num2) > tolerance:
                        mismatches += 1
                except ValueError:
                    print(f"Non-numeric data at row {row_idx}, column {col_idx}: {val1} vs {val2}")
    return mismatches


cwd = os.getcwd()
pardir = os.pardir
dfdm_dir = os.path.abspath(os.path.join(cwd,pardir))
exe = os.path.join(dfdm_dir,"build/dfdm")
config_path = os.path.join(dfdm_dir,"tests/ref_config/config_test.toml")
output_path = os.path.join(cwd,"tests/test_out_mpi/")
stdout_capture = os.path.join(cwd,"tests/test_out_mpi/test_2.out")
reference = os.path.join(dfdm_dir, "tests/full_sim_reference/")
domain_out = os.path.join(cwd, "tests/test_out_mpi/domain/")
domain_in = os.path.join(cwd,"tests/ref_mesh_4/")



print(output_path)
if os.path.exists(output_path):
    shutil.rmtree(output_path)
os.mkdir(output_path)
# proc = subprocess.run([exe, config_path, output_path, ">>", stdout_capture], shell=True)
cmd = str("srun -n4 ")+exe+str(" -c ")+config_path+str(" -o ")+output_path+str(" -m ")+domain_in+str(" -d ")+domain_out+str(" ")+">> "+stdout_capture
print(f"executing:{cmd}")
os.system(cmd)

elems = [2,7]
steps = [0,500,1000,1500,2000,2500,3000,3050,3650,3500,4000,4550,5500,5850,6000,7000,8650,9950]
count = 0
for x in elems:
    for z in steps:
        ref_file = os.path.join(reference,"elem_"+str(x)+"_"+str(z)+".out")
        test_file = os.path.join(output_path,"elem_"+str(x)+"_"+str(z)+".out")
        results = compare_csv_matrices(ref_file,test_file)
        if results != 0:
            count = count + 1
            print(f"Last file compared had mismatches:{results}, total {count} mismatched files found in simulation")

if count == 0:
    print("Full Simulation test Passed!")
else:
    print("Full Simulation test Failed!")   

