import numpy as np
print('Universe        = vanilla')
print('Notification    = Never')
print('# Run this exe with these args')
print('Executable      = /gpfs01/astro/workarea/lmezini/deblender_tests/new_galsim_deblend.py')
print('Image_Size       =  1000000')
print('GetEnv = True')
print('kill_sig        = SIGINT')
print('#requirements = (cpu_experiment == "star") || (cpu_experiment == "phenix")')
print('requirements = (cpu_experiment == "phenix")')
print('+Experiment     = "astro"')

for i in range(101,201):
    seed = np.random.randint(0, 2**15)
    print('+job_name = "run106-'+str(i).zfill(6)+'"')
    print('Arguments = /gpfs01/astro/workarea/lmezini/scarlet-tests/run106/run106-output-'+str(i).zfill(6)+'.fits 2000 '+str(seed)+' /gpfs01/astro/workarea/lmezini/deblender_tests/config_v6.yaml')
    print('Queue')
