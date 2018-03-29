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

for i in range(401,501):
    seed = np.random.randint(0, 2**15)
    print('+job_name = "run102_5-'+str(i).zfill(6)+'"')
    print('Arguments = /gpfs01/astro/workarea/lmezini/scarlet-tests/run102/run102_5-output-'+str(i).zfill(6)+'.fits 2000 '+str(seed)+' /gpfs01/astro/workarea/lmezini/deblender_tests/config_v2.yaml')
    print('Queue')
