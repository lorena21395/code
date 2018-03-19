#!/usr/bin/env python
import yaml
from minimof import minimof
from scipy.misc import imsave
import fitsio
import galsim
import argparse
import sep
import ngmix
import numpy as np
import scarlet
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib.patches import Ellipse

parser = argparse.ArgumentParser()
parser.add_argument("outfile",help="Output file name and path")
parser.add_argument("ntrials",help="Number of trials to be run",type = int)
parser.add_argument("seed",help="Seed for random number generator",type = int)
parser.add_argument("mode",help="Specify whether scarlet, minimof, or control",type = str)
args = parser.parse_args()
outfile_name = args.outfile
ntrial = args.ntrials
seed = args.seed
mode = args.mode
np.random.seed(seed)
rng = np.random.RandomState(seed=np.random.randint(0,2**30))

Sim_specs = {'Cen': {'Type':'Gaussian','hlr':1.7,'Flux':6000.,'Pos':'Fixed','dx':0.,'dy':0.},
            'Neigh':{'Type':'Exponential','hlr':3.4,'Flux':8000.,'Pos':'Rand_circle','dx':12.,'dy':12.},
            'Image':{'Dims':[50,50],'Bgrms':10.,'Scale':1.},
            'Psf':{'Psf_hlr':1.7,'Scale':1.,'Bgrms_psf':0.0001},
            'Shear':{'Shear1':0.02,'Shear2':0.00}
             }

class Simulation(dict):
    """Simulate galaxy images with noise and psf"""

    def __init__(self,specs,mode):
        self.update(specs)
        self.mode = mode

    def _get_psf_img(self):
        psf_hlr = self['Psf']['Psf_hlr']
        psf = galsim.Gaussian(half_light_radius = psf_hlr)
        psf_gsim = psf.drawImage(scale = self['Psf']['Scale'])
        psf_im = psf_gsim.array
        noise_psf = np.random.normal(scale=self['Psf']['Bgrms_psf'],size=psf_im.shape)
        psf_im += noise_psf

        return psf,psf_im

    def _get_gals(self):
        Cen = galsim.Gaussian(half_light_radius=self['Cen']['hlr'],flux=self['Cen']['Flux'])
        Neigh = galsim.Exponential(half_light_radius=self['Neigh']['hlr'],flux=self['Neigh']['Flux'])

        if self['Cen']['Pos'] == 'Fixed':
            dx1,dy1 = self['Cen']['dx'],self['Cen']['dy']
        elif self['Cen']['Pos'] == 'Rand_Circle':
            theta = 2.*np.pi*np.random.random()
            dx1,dy1 = self['Cen']['dx']*np.cos(theta),self['Cen']['dy']*np.sin(theta)
        if self['Neigh']['Pos'] == 'Fixed':
            dx2,dy2 = self['Neigh']['dx'],self['Neigh']['dy']
        elif self['Neigh']['Pos'] == 'Rand_circle':
            theta = 2.*np.pi*np.random.random()
            dx2,dy2 = self['Neigh']['dx']*np.cos(theta),self['Neigh']['dy']*np.sin(theta)

        Cen = Cen.shift(dx=dx1, dy=dy1)
        Neigh = Neigh.shift(dx=dx2, dy=dy2)
        if mode == 'scarlet' or mode == 'minimof':
            gals = [Cen, Neigh]
            objs = galsim.Add(gals)
        elif mode == 'control':
            gals = [Cen]
            objs = galsim.Add(gals)

        shear1, shear2 = self['Shear']['Shear1'],self['Shear']['Shear2']
        objs = objs.shear(g1=shear1, g2=shear2)

        dims = self['Image']['Dims']
        coord1 = (dy1+(dims[0]-1.)/2.,dx1+(dims[1]-1.)/2.)
        coord2 = (dy2+(dims[0]-1.)/2.,dx2+(dims[1]-1.)/2.)
        coords = [coord1,coord2]
        
        return objs,coords

    def _get_noise(self):
        dims = self['Image']['Dims']
        noise = np.random.normal(scale=self['Image']['Bgrms'],size=(dims[0],dims[1]))
        return noise

    def __call__(self):
        dims = self['Image']['Dims']
        bg_rms = self['Image']['Bgrms']
        psf, psf_im = self._get_psf_img()
        objs,coords = self._get_gals()
        objs = galsim.Convolve(objs, psf)
        gsim = objs.drawImage(nx=dims[1], ny=dims[0], scale=self['Image']['Scale'])
        im = gsim.array
        noise = self._get_noise()
        im += noise

        if mode ==  'scarlet' or mode == 'control':
            im = im.reshape( (1, dims[1], dims[0]) )
        return im,psf_im,coords

class Model(Simulation):
    
    def __init__(self):
        Simulation.__init__(self,Sim_specs,mode)
    
    def _get_model(self):
        im,psf_im,coords = Simulation.__call__(self)
        bg_rms = self['Image']['Bgrms']
        if mode == 'scarlet':
            constraints = {"S": None, "m": {'use_nearest': False}, "+": None}
            sources = [scarlet.ExtendedSource(coord, im, [bg_rms]) for coord in coords]
            blend = scarlet.Blend(sources, im, bg_rms=[bg_rms])
            blend.fit(10000, e_rel=1e-1)
            model = blend.get_model()
            mod1 = blend.get_model(m=0)
            mod2 = blend.get_model(m=1)
            cen_mod = sources[0].get_model()
            neigh_mod = sources[1].get_model()
            steps_used = blend.it
        return im,psf_im,model,mod1,mod2,cen_mod,neigh_mod,steps_used,coords
    
    def _rob_deblend(self,im,model,mod1,mod2):
        C = [np.zeros([50,50]),np.zeros([50,50])]
        I = im
        w = np.array([model[0,:,:],model[0,:,:]])
        T = np.array([mod1[0,:,:],mod2[0,:,:]])
        mod_sum = np.zeros([50,50])
        for r in range(2):
            mod_sum += w[r]*T[r]
        mod_sum += 0.000001
        for r in range(2):
            C[r] = I*np.divide(w[r]*T[r],mod_sum)
            #print(np.linalg.inv(mod_sum))
        return C
    
    def _readd_noise(self,coords,neigh_mod,cen_obj):
        neigh_shape = neigh_mod.shape
        coord2 = coords[1]
        bg_rms = self['Image']['Bgrms']
        region = cen_obj[int(coord2[0]-(neigh_shape[1]/2)+1):int(coord2[0]+(neigh_shape[1]/2)+1),int(coord2[1]-(neigh_shape[2]/2)+1):int(coord2[1]+(neigh_shape[2]/2)+1)]
        extra_noise = np.sqrt(bg_rms**2 - np.var(region))
        cen_obj[int(coord2[0]-(neigh_shape[1]/2)+1):int(coord2[0]+(neigh_shape[1]/2)+1),int(coord2[1]-(neigh_shape[2]/2)+1):int(coord2[1]+(neigh_shape[2]/2)+1)] += rng.normal(scale=extra_noise)
        
        return cen_obj

def observation(image,sigma,row,col,psf_sigma,psf_im):
    # sigma is the standard deviation of the noise
    weight = image*0 + 1.0/sigma**2
    # row,col are the center of the object in the model
    # image, so we need to figure that out somehow
    jacob = ngmix.UnitJacobian(row=row, col=col)
    psf_weight = psf_im*0 + 1.0/psf_sigma**2
    # psf should be centered
    cen = [(np.array(psf_im.shape[0]) - 1.0)/2.0,(np.array(psf_im.shape[1]) - 1.0)/2.0]
    psf_jacob = ngmix.UnitJacobian(row=cen[0], col=cen[1])
    psf_obs = ngmix.Observation(psf_im, weight=psf_weight, jacobian=psf_jacob)
    obs = ngmix.Observation(image,weight=weight,jacobian=jacob,psf=psf_obs)

    return obs

def get_prior():
    cen_sigma = 1.0
    cen_prior = ngmix.priors.CenPrior(
        0.0, 0.0,
        cen_sigma,
        cen_sigma,
    )
    g_prior = ngmix.priors.GPriorBA(0.2)

    T_pars =  [-10.0, 0.03, 1.0e+06, 1.0e+05]
    T_prior = ngmix.priors.TwoSidedErf(*T_pars)

    flux_pars = [-1.0e+04, 1.0, 1.0e+09, 0.25e+08]
    flux_prior = ngmix.priors.TwoSidedErf(*flux_pars)

    prior = ngmix.joint_prior.PriorSimpleSep(
        cen_prior,
        g_prior,
        T_prior,
        flux_prior,
    )
    return prior

def norm_test():
    Mod = Model()
    im,psf_im,model,mod1,mod2,cen_mod,neigh_mod,steps_used,coords = Mod._get_model()
    C = Mod._rob_deblend(im,model,mod1,mod2)
    cen_obj = im[0,:,:] - C[1][0,:,:]#mod2[0,:,:]
    cen_obj = Mod._readd_noise(coords,neigh_mod,cen_obj)
    plt.imshow(cen_obj,interpolation='nearest', cmap='gray',vmin = np.min(cen_obj),vmax= np.max(cen_obj))
    plt.colorbar()
    plt.savefig("test_4.png")
    dobs = observation(cen_obj,Mod['Image']['Bgrms'],coords[0][1],
                       coords[0][0],Mod['Psf']['Bgrms_psf'],psf_im)
    return dobs

def do_metacal(psf_model,gal_model,max_pars,psf_Tguess,prior,
                         ntry,metacal_pars,output,dobs):
    boot = ngmix.bootstrap.MaxMetacalBootstrapper(dobs)
    boot.fit_psfs(
            psf_model,
            psf_Tguess,
            ntry=ntry,
            fit_pars=max_pars['lm_pars'],
            skip_already_done=False,
         )
    boot.fit_metacal(psf_model,gal_model,max_pars,psf_Tguess,prior=prior,
                     ntry=ntry,metacal_pars=metacal_pars,)
    res = boot.get_metacal_result()
    #print("flags:",res['mcal_flags'])                                      
    output['flags'][j] = res['mcal_flags']
    output['pars'][j] = res['noshear']['pars']
    output['pars_1p'][j] = res['1p']['pars']
    output['pars_1m'][j] = res['1m']['pars']
    output['pars_2p'][j] = res['2p']['pars']
    output['pars_2m'][j] = res['2m']['pars']
    #output['noise_std'][j] = reg_std

    return output

dt = [
    ('flags','i4'),
    ('pars','f8',6),
    ('pars_1p','f8',6),
    ('pars_1m','f8',6),
    ('pars_2p','f8',6),
    ('pars_2m','f8',6),
    ('noise_std','f8'),
]

psf_model = 'gauss'
gal_model = 'gauss'
psf_Tguess = 4.0
psf_fit_pars = {'maxfev': 2000}
ntry=2
max_pars = {
    'method': 'lm',
    'lm_pars': {
        'maxfev': 2000,
        'xtol': 5.0e-5,
        'ftol': 5.0e-5,
    }
}

metacal_pars = {
    'symmetrize_psf': True,
    'types': ['noshear','1p','1m','2p','2m'],
}
prior = get_prior()
output = np.zeros(ntrial, dtype=dt)

for j in range(ntrial):
    print(j)
    try:
        dobs = norm_test()
        output = do_metacal(psf_model,gal_model,max_pars,
                         psf_Tguess,prior,ntry,
                         metacal_pars,output,dobs)
    except (np.linalg.linalg.LinAlgError,ValueError):
        print("error")
        output['flags'][j] = 2

fitsio.write(outfile_name, output, clobber=True)
