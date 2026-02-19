import subprocess
import os
import shutil
import filecmp

cwd = os.getcwd()
pardir = os.pardir
dfdm_dir = os.path.abspath(os.path.join(cwd,pardir))
exe = os.path.join(dfdm_dir,"build/dfdm")
config_path = os.path.join(dfdm_dir,"tests/ref_config/config_test_source_location.toml")
output_path = os.path.join(cwd,"tests/test_source_location_out/")
stdout_capture = os.path.join(cwd,"tests/test_source_location_out/test_source_location.out")
# reference = os.path.join(dfdm_dir, "tests/full_sim_reference/")
domain_out = os.path.join(cwd, "tests/test_source_location_out/domain/")
domain_in = os.path.join(cwd,"tests/ref_mesh_4/")


print(output_path)
if os.path.exists(output_path):
    shutil.rmtree(output_path)
os.mkdir(output_path)
# proc = subprocess.run([exe, config_path, output_path, ">>", stdout_capture], shell=True)
cmd = str("srun -n4 ")+exe+str(" -c ")+config_path+str(" -o ")+output_path+str(" -m ")+domain_in+str(" -d ")+domain_out+str(" ")+">> "+stdout_capture
print(f"executing:{cmd}")
os.system(cmd)

f = open(stdout_capture,'r')
lines=f.readlines()

passed=False

for line in lines: 
    if "[DFDM::LOGGER::] [info] DFDM::Source::get_source_location: xlocal: " in line:
        xlocal = float(line.split("xlocal: ")[1].split(" ")[0])
        zlocal = float(line.split("zlocal: ")[1].split(" ")[0])
        if xlocal == 0.5 and zlocal == 0.5:
            passed=True
        else:   
            print(f"Test failed: xlocal={xlocal}, zlocal={zlocal}")
            passed=False

    if "[DFDM::LOGGER::] [info] DFDM::Source::get_source_location: xglobal: " in line:
        xglobal = float(line.split("xglobal: ")[1].split(" ")[0])
        zglobal = float(line.split("zglobal: ")[1].split(" ")[0])
        if abs(xglobal)<1 and abs(zglobal)<1:
            passed=True
        else:   
            print(f"Test failed: xglobal={xglobal}, zglobal={zglobal}")
            passed=False

    if "[DFDM::LOGGER::] [info] DFDM::Main: Receiver location in rank:" in line:
        xglobal = float(line.split("xglobal:")[1].split(",")[0])
        zglobal = float(line.split("zglobal:")[1].split()[0])
        print(zglobal)
        if abs(xglobal)<1 and abs(zglobal)<1:
            passed=True
        else:   
            print(f"Test failed: xglobal={xglobal}, zglobal={zglobal}")
            passed=False

if passed:
    print("Full Simulation test Passed!")
else:
    print("Full Simulation test Failed!")   