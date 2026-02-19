import subprocess
import os
import shutil
import filecmp

cwd = os.getcwd()
pardir = os.pardir
dfdm_dir = os.path.abspath(os.path.join(cwd,pardir))
exe = os.path.join(dfdm_dir,"build/dfdm")
config_path = os.path.join(dfdm_dir,"tests/ref_config/config_test.toml")
output_path = os.path.join(cwd,"tests/generated_mats/")
stdout_capture = os.path.join(cwd,"tests/generated_mats/bspline.out")
reference = os.path.join(dfdm_dir, "tests/ref_matrices/")
domain_out = os.path.join(cwd, "tests/generated_mats/domain/")
domain_in = os.path.join(cwd,"tests/ref_mesh_single/")



print(output_path)
if os.path.exists(output_path):
    shutil.rmtree(output_path)
os.mkdir(output_path)
# proc = subprocess.run([exe, config_path, output_path, ">>", stdout_capture], shell=True)
cmd = exe+str(" -c ")+config_path+str(" -o ")+output_path+str(" -m ")+domain_in+str(" -d ")+domain_out+str(" --enable-tests-run true ")+">> "+stdout_capture
print(f"executing:{cmd}")
os.system(cmd)

elems = [2,7]
count = 0
for x in elems:
    ref_file = os.path.join(reference,"B1_X_"+str(x))
    test_file = os.path.join(output_path,"B1_X_"+str(x))
    results = filecmp.cmp(ref_file,test_file)
    if results == False:
        count = count + 1
        print(f"Found {count} mismatches in simulation")

    ref_file = os.path.join(reference,"B2_X_"+str(x))
    test_file = os.path.join(output_path,"B2_X_"+str(x))        
    results = filecmp.cmp(ref_file,test_file)
    if results == False:
        count = count + 1
        print(f"Found {count} mismatches in simulation")

    ref_file = os.path.join(reference,"B1_Z_"+str(x))
    test_file = os.path.join(output_path,"B1_Z_"+str(x))
    results = filecmp.cmp(ref_file,test_file)
    if results == False:
        count = count + 1
        print(f"Found {count} mismatches in simulation")

    ref_file = os.path.join(reference,"B2_Z_"+str(x))
    test_file = os.path.join(output_path,"B2_Z_"+str(x))        
    results = filecmp.cmp(ref_file,test_file)
    if results == False:
        count = count + 1
        print(f"Found {count} mismatches in simulation")


if count == 0:
    print("Full Simulation test Passed!")
else:
    print("Full Simulation test Failed!")   
