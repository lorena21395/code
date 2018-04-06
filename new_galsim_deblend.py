#!/usr/bin/env python
import yaml
#from minimof import minimof
from scipy import signal
import fitsio
import galsim
import argparse
import sep
import ngmix
import numpy as np
import scarlet
import matplotlib.pyplot as plt
plt.switch_backend('agg')

parser = argparse.ArgumentParser()
parser.add_argument("outfile",help="Output file name and path")
parser.add_argument("ntrials",help="Number of trials to be run",type = int)
parser.add_argument("seed",help="Seed for random number generator",type = int)
parser.add_argument("config",help="Configuration for the simulation")

args = parser.parse_args()
outfile_name = args.outfile
ntrial = args.ntrials
seed = args.seed
config_file = args.config

np.random.seed(seed)
rng = np.random.RandomState(seed=np.random.randint(0,2**30))

with open(config_file) as fobj:
    Sim_specs =yaml.load(fobj)


class Simulation(dict):
    """Simulate galaxy images with noise and psf"""

    def __init__(self,specs):
        self.update(specs)
    
    def _get_mode(self):
        mode = self['Mode']
        return mode

    def _get_psf_img(self):
        psf_hlr = self['Psf']['Psf_hlr']
        psf = galsim.Gaussian(half_light_radius = psf_hlr)
        psf_gsim = psf.drawImage(scale = self['Psf']['Scale'])
        psf_im = psf_gsim.array
        noise_psf = np.random.normal(scale=self['Psf']['Bgrms_psf'],size=psf_im.shape)
        psf_im += noise_psf

        return psf,psf_im

    def _get_gals(self):
        mode = self['Mode']
        Cen = galsim.Gaussian(half_light_radius=self['Cen']['hlr'],flux=self['Cen']['Flux'])
        Neigh = galsim.Exponential(half_light_radius=self['Neigh']['hlr'],flux=self['Neigh']['Flux'])

        #cen position
        if self['Cen']['Pos'] == 'Fixed':
            dx1,dy1 = self['Cen']['dx'],self['Cen']['dy']
        
        elif self['Cen']['Pos'] == 'Rand_Circle':
            r = self['Cen']['r']
            theta = 2.*np.pi*np.random.random()
            dx1,dy1 = r*np.cos(theta),r*np.sin(theta)
        
        #neigh position
        if self['Neigh']['Pos'] == 'Fixed':
            dx2,dy2 = self['Neigh']['dx'],self['Neigh']['dy']
        
        elif self['Neigh']['Pos'] == 'Rand_circle':
            r = self['Neigh']['r']
            theta = 2.*np.pi*np.random.random()
            dx2,dy2 = r*np.cos(theta),r*np.sin(theta)

        dx1 +=np.random.uniform(low=-0.5, high=0.5)
        dy1 +=np.random.uniform(low=-0.5, high=0.5)
        dx2 +=np.random.uniform(low=-0.5, high=0.5)
        dy2 +=np.random.uniform(low=-0.5, high=0.5)
        
        Cen = Cen.shift(dx=dx1, dy=dy1)
        Neigh = Neigh.shift(dx=dx2, dy=dy2)
        
        if mode == 'scarlet' or mode == 'minimof':
            gals = [Cen]#, Neigh]
            objs = galsim.Add(gals)
        elif mode == 'control':
            gals = [Cen]
            objs = galsim.Add(gals)

        shear1, shear2 = self['Shear']['Shear1'],self['Shear']['Shear2']
        objs = objs.shear(g1=shear1, g2=shear2)
        
        return objs,dx1,dy1,dx2,dy2

    def _get_noise(self,dims,bg_rms):
        noise = np.random.normal(scale=bg_rms,size=(dims[0],dims[1]))
        
        return noise

    def __call__(self):
        mode = self['Mode']
        bg_rms = self['Image']['Bgrms']
        psf, psf_im = self._get_psf_img()
        objs,dx1,dy1,dx2,dy2 = self._get_gals()
        objs = galsim.Convolve(objs, psf)
        gsim = objs.drawImage(scale=self['Image']['Scale'])
        im = gsim.array
        dims = np.shape(im)
        
        #rewrite coords for galsim
        cen =  (np.array(im.shape) - 1.0)/2.0
        coord1 = (dy1+cen[1],dx1+cen[0])
        coord2 = (dy2+cen[1],dx2+cen[0])
        coords = [coord1]#,coord2]

        noise = self._get_noise(dims,bg_rms)
        im += noise
        noise = self._get_noise(dims,bg_rms)

        if mode ==  'scarlet' or mode == 'control':
            im = im.reshape( (1, dims[0], dims[1]) )
        return im,psf_im,coords,dims,dx1,dy1,noise

class Model(Simulation):
    
    def __init__(self):
        Simulation.__init__(self,Sim_specs)
    
    def _get_model(self):
        im,psf_im,coords,dims,dx1,dy1,noise = Simulation.__call__(self)
        bg_rms = self['Image']['Bgrms']
        mode = self['Mode']
        if mode == 'scarlet':
            constraints = {"S": None, "m": {'use_nearest': False}, "+": None}
            sources = [scarlet.ExtendedSource(coord, im, [bg_rms]) for coord in coords]
            config = scarlet.Config(edge_flux_thresh=0.05)
            blend = scarlet.Blend(sources, im, bg_rms=[bg_rms],config=config)
            blend.fit(10000, e_rel=1e-3)
            model = blend.get_model()
            mod1 = blend.get_model(m=0)
            #mod2 = blend.get_model(m=1)
            cen_mod = sources[0].get_model()
            #output['mod_size_flag'][j] = 0
            #if np.shape(cen_mod) != (1,25,25):
            #    output['mod_size_flag'][j] = 3
#neigh_mod = sources[1].get_model()
            #steps_used = blend.it

       # return im,psf_im,model,mod1,mod2,cen_mod,neigh_mod,coords,dx1,dy1
        return im,psf_im,model,mod1,cen_mod,coords,dx1,dy1,noise
    
    def _rob_deblend(self,im,model,mod1,dims):
        C = np.zeros((dims[0],dims[1],1))#,2))
        W = np.zeros((dims[0],dims[1],1))#2))
        I = im
        w = np.array([model[0,:,:]])#,model[0,:,:]])
        T = np.array([mod1[0,:,:]])#,mod2[0,:,:]])
        mod_sum = np.zeros(dims)
        for r in range(1):
            mod_sum += w[r]*T[r] 
        zeros = np.where(mod_sum == 0.)
        mod_sum[zeros] += 0.000001
        for r in range(1):
            W[:,:,r] = np.divide(w[r]*T[r],mod_sum)
            C[:,:,r] = I*np.divide(w[r]*T[r],mod_sum)
        return C,W

    def _readd_noise(self,cen_obj,W):
        bg_rms = self['Image']['Bgrms']
        extra_noise = np.sqrt((bg_rms**2)*(1 - W**2))
        cen_obj += extra_noise*rng.normal(size=cen_obj.shape)
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
    mode = Mod._get_mode()
    if mode == 'scarlet':
        im,psf_im,model,mod1,cen_mod,coords,dx1,dy1,noise = Mod._get_model()
        #im,psf_im,model,mod1,mod2,cen_mod,neigh_mod,coords,dx1,dy1 = Mod._get_model()
        cen_shape = cen_mod.shape
        coord1 = coords[0]
        dims = [np.shape(im)[1],np.shape(im)[2]]
        C,W = Mod._rob_deblend(im,model,mod1,dims)
        Cnoise,Wnoise = Mod._rob_deblend(noise,model,mod1,dims)        
        half1 = cen_shape[1]/2.
        half2 = cen_shape[2]/2.

        beg1 = int(coord1[0]-half1+1)
        end1 = int(coord1[0]+half1+1)

        beg2 = int(coord1[1]-half2+1)
        end2 = int(coord1[1]+half2+1)


        if cen_shape[1] != cen_shape[2]:
            cen_obj = np.zeros((max(cen_shape),max(cen_shape)))
            weights = np.zeros((max(cen_shape),max(cen_shape)))
            noise = np.zeros((max(cen_shape),max(cen_shape)))
            
            cen_obj[0:cen_shape[1],0:cen_shape[2]] = C[beg1:end1,beg2:end2,0]
            weights[0:cen_shape[1],0:cen_shape[2]] = W[beg1:end1,beg2:end2,0]
            noise[0:cen_shape[1],0:cen_shape[2]] = Cnoise[beg1:end1,beg2:end2,0]

        else:
            cen_obj = C[beg1:end1,beg2:end2,0]
            weights = W[beg1:end1,beg2:end2,0]
            noise = Cnoise[beg1:end1,beg2:end2,0]

        cen_obj_w_noise = Mod._readd_noise(cen_obj,weights)
        noise_w_noise = Mod._readd_noise(noise,weights)

        corr = signal.correlate2d(im[0,beg1:end1,beg2:end2],model[0,beg1:end1,beg2:end2])
        #output['dims'][j] = np.array(dims)
    
        plt.imshow(corr,interpolation='nearest', cmap='gray',vmin = np.min(corr),vmax= np.max(corr/2000))
        plt.title("True and Model Corr")
        plt.colorbar();
        plt.tight_layout()
        plt.savefig("test.png")
        
        #coords of object in region
        new_coords = (dx1+(cen_obj.shape[0]-1.0)/2.0,dy1+(cen_obj.shape[1]-1.0)/2.0)
        #dobs = observation(cen_obj_w_noise,Mod['Image']['Bgrms'],new_coords[1],
        #               new_coords[0],Mod['Psf']['Bgrms_psf'],psf_im)
        
        #dobs.noise = noise_w_noise

    elif mode == 'control':
        im,psf_im,coords,dims,dx1,dy1 = Mod.__call__()
        output['dims'][j] = np.array(dims)
        
        dobs = observation(im[0],Mod['Image']['Bgrms'],coords[0][1],
                       coords[0][0],Mod['Psf']['Bgrms_psf'],psf_im)
    
    
    return corr#dobs#,cen_obj,cen_obj_w_noise

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

    return output

dt = [
    ('flags','i4'),
    ('pars','f8',6),
    ('pars_1p','f8',6),
    ('pars_1m','f8',6),
    ('pars_2p','f8',6),
    ('pars_2m','f8',6),
    ('dims','f8',2),
    ('mod_size_flag','i4'),
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
    'use_noise_image': True,
    'types': ['noshear','1p','1m','2p','2m'],
}
prior = get_prior()
#output = np.zeros(ntrial, dtype=dt)

dt = [('corr','f8',(29,29)),('corr_big','f8',(49,49))]      
output = np.zeros(1, dtype=dt)


corr_vals = []
corr_vals_big = []
for j in range(ntrial):
    print(j)
    try:
        corr = norm_test()
        if np.shape(corr) == (29,29):
            corr_vals.append(corr)
        elif np.shape(corr) == (49,49):
            corr_vals_big.append(corr)
        #pix_vals_w_noise.append(cen_obj_w_noise)
        #pix_vals_std = np.std(np.array(pix_vals,axis=0))
        #pix_vals_w_noise_std = np.std(np.array(pix_vals_w_noise,axis=0))
        #output['std'] = pix_vals_std
        #output['noisy_std'] = pix_vals_w_noise_std
        #out = do_metacal(psf_model,gal_model,max_pars,
        #                 psf_Tguess,prior,ntry,
        #                 metacal_pars,output,dobs)
    except (np.linalg.linalg.LinAlgError,ValueError):
        print("error")
        #output['flags'][j] = 2
if len(corr_vals) >= 1:
    avg_corr = np.mean(np.array(corr_vals),axis=0)
    output['corr'] = avg_corr
if len(corr_vals_big) >= 1:
    avg_corr_big = np.mean(np.array(corr_vals_big),axis=0) 
    output['corr_big'] = avg_corr_big
fitsio.write(outfile_name, output, clobber=True)
