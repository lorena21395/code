import numpy as np
print('Universe        = vanilla')
print('Notification    = Never')
print('# Run this exe with these args')
print('Executable      = /gpfs01/astro/workarea/lmezini/deblender_tests/code/new_galsim_deblend.py')
print('Image_Size       =  1000000')
print('GetEnv = True')
print('kill_sig        = SIGINT')
#print('#requirements = (cpu_experiment == "star") || (cpu_experiment == "phenix")')
#print('requirements = (cpu_experiment == "phenix")')
print('+Experiment     = "astro"')

for i in range(1,2501):
    seed = np.random.randint(1, 2**15)
    print('+job_name = "run196-'+str(i).zfill(6)+'"')
    print('Arguments = /gpfs01/astro/workarea/lmezini/scarlet-tests/run196/run196_2-output-'+str(i).zfill(6)+'.fits 400 '+str(seed)+' /gpfs01/astro/workarea/lmezini/deblender_tests/config_files/config-scar-sers-2step-N0.1-r07-3band.yaml')
    print('Queue')
