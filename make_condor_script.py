import numpy as np
print('Universe        = vanilla')
print('Notification    = Never')
print('# Run this exe with these args')
print('Executable      = /gpfs01/astro/workarea/lmezini/deblender_tests/new_galsim_deblend.py')
print('Image_Size       =  1000000')
print('GetEnv = True')
print('kill_sig        = SIGINT')
#print('#requirements = (cpu_experiment == "star") || (cpu_experiment == "phenix")')
#print('requirements = (cpu_experiment == "phenix")')
print('+Experiment     = "astro"')

for i in range(1,2501):
    seed = np.random.randint(1, 2**15)
    print('+job_name = "run201-'+str(i).zfill(6)+'"')
    print('Arguments = /gpfs01/astro/workarea/lmezini/scarlet-tests/run201/run201_1-output-'+str(i).zfill(6)+'.fits 400 '+str(seed)+' /gpfs01/astro/workarea/lmezini/deblender_tests/config_files/config_scar_gau_gau_1step_N0.1.yaml')
    print('Queue')
