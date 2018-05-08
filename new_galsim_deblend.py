#!/usr/bin/env python
#####
#individually shear each galaxy first
#use sersic galaxies
#change minimof gal model to "cm"
#PSF matching

import yaml
import minimof
import fitsio
import galsim
import argparse
import sep
import ngmix
import numpy as np

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
        psf_gsim = psf.drawImage(nx = 25, ny=25, scale = self['Psf']['Scale'])
        psf_im = psf_gsim.array
        noise_psf = np.random.normal(scale=self['Psf']['Bgrms_psf'],size=psf_im.shape)
        psf_im += noise_psf

        return psf,psf_im

    def _get_gals(self):
        mode = self['Mode']
        #Cen = galsim.Gaussian(half_light_radius=self['Cen']['hlr'],flux=self['Cen']['Flux'])
        #Neigh = galsim.Exponential(half_light_radius=self['Neigh']['hlr'],flux=self['Neigh']['Flux'])

        n=rng.uniform(low=0.5, high=4)
        Cen = galsim.Sersic(n,half_light_radius=self['Cen']['hlr'],flux=self['Cen']['Flux'])
        n=rng.uniform(low=0.5, high=4)
        Neigh = galsim.Sersic(n,half_light_radius=self['Neigh']['hlr'],flux=self['Neigh']['Flux'])
        
        Cen = Cen.shear(g1=np.random.normal(scale=0.02),g2=np.random.normal(scale=0.02))
        Neigh = Neigh.shear(g1=np.random.normal(scale=0.02),g2=np.random.normal(scale=0.02))

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
        
        if mode == 'scarlet' or mode == 'mof':
            gals = [Cen, Neigh]
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
        if 'dims' in self['Image']:
            ny,nx=self['Image']['dims']
        else:
            ny, nx = None, None

        gsim = objs.drawImage(
            nx=nx,
            ny=ny,
            scale=self['Image']['Scale'],
        )
        im = gsim.array
        dims = np.shape(im)
        
        #rewrite coords for galsim
        cen =  (np.array(im.shape) - 1.0)/2.0
        coord1 = (dy1+cen[1],dx1+cen[0])
        coord2 = (dy2+cen[1],dx2+cen[0])
        coords = [coord1,coord2]

        noise = self._get_noise(dims,bg_rms)
        im += noise
        noise = self._get_noise(dims,bg_rms)

        if mode ==  'scarlet' or mode == 'control':
            im = im.reshape( (1, dims[0], dims[1]) )

        return im,psf_im,coords,dims,dx1,dy1,noise

class Model(Simulation):
    
    def __init__(self):
        Simulation.__init__(self,Sim_specs)
    
    def _get_scar_model(self):
        import scarlet

        im,psf_im,coords,dims,dx1,dy1,noise = Simulation.__call__(self)
        bg_rms = self['Image']['Bgrms']
        mode = self['Mode']
        constraints = {"S": None, "m": {'use_nearest': False}, "+": None}
        #constraints['l0'] = bg_rms
        psf_dims = np.shape(psf_im)
        psf_im3d = psf_im.reshape( (1, psf_dims[0], psf_dims[1]) )
        target_psf = scarlet.psf_match.fit_target_psf(psf_im3d, scarlet.psf_match.gaussian)
        diff_kernels, psf_blend = scarlet.psf_match.build_diff_kernels(psf_im3d, target_psf)
        sources = [scarlet.ExtendedSource(coord, im, [bg_rms],psf=diff_kernels) for coord in coords]
        #sources = [scarlet.ExtendedSource(coord, im, [bg_rms]) for coord in coords]
        #scarlet.ExtendedSource.shift_center=0.0
        #config = scarlet.Config(edge_flux_thresh=0.05)
        blend = scarlet.Blend(sources, im, bg_rms=[bg_rms])#,config=config)
        blend.fit(10000, e_rel=1e-3)
        model = blend.get_model()
        mod1 = blend.get_model(m=0)
        mod2 = blend.get_model(m=1)
        cen_mod = sources[0].get_model()
        neigh_mod = sources[1].get_model()
        #steps_used = blend.it
        
        return im,psf_im,model,mod1,mod2,cen_mod,neigh_mod,coords,dx1,dy1,noise
        #return im,psf_im,model,mod1,cen_mod,coords,dx1,dy1,noise 


    def _fit_psf_admom(self, obs):
        Tguess=4.0
        am=ngmix.admom.run_admom(obs, Tguess)
        return am.get_gmix()

    def _get_mof_obs(self):
        im,psf_im,coords,dims,dx1,dy1,noise = self()

        bg_rms = self['Image']['Bgrms']
        bg_rms_psf = self['Psf']['Bgrms_psf']

        psf_ccen=(np.array(psf_im.shape)-1.0)/2.0
        psf_jacob = ngmix.UnitJacobian(
            row=psf_ccen[0],
            col=psf_ccen[1],
        )
        psf_weight=psf_im*0 + 1.0/bg_rms_psf**2
        psf_obs = ngmix.Observation(
            psf_im,
            weight=psf_weight,
            jacobian=psf_jacob,
        )

        psf_gmix=self._fit_psf_admom(psf_obs)
        psf_obs.set_gmix(psf_gmix)

        weight=im*0 + 1.0/bg_rms**2
        jacobian=ngmix.UnitJacobian(
            row=0,
            col=0,
        )

        obs = ngmix.Observation(
            im,
            weight=weight,
            jacobian=jacobian,
            psf=psf_obs,
        )
        obs.noise=noise

        return obs, coords

    def _get_mof_guess(self, coord_list):
        npars_per=7
        num=len(coord_list)
        assert num==2,"two objects for now"

        npars_tot = num*npars_per
        guess = np.zeros(npars_tot)

        for i,coords in enumerate(coord_list):
            row, col = coords

            beg=i*npars_per

            if i==0:
                F = self['Cen']['Flux']
            else:
                F = self['Neigh']['Flux']


            # always close guess for center
            guess[beg+0] = row + rng.uniform(low=-0.05, high=0.05)
            guess[beg+1] = col + rng.uniform(low=-0.05, high=0.05)

            # always arbitrary guess for shape
            guess[beg+2] = rng.uniform(low=-0.05, high=0.05)
            guess[beg+3] = rng.uniform(low=-0.05, high=0.05)

            # we could get a better guess from sep
            T = 2.0
            guess[beg+4] = T*(1.0 + rng.uniform(low=-0.05, high=0.05))

            guess[beg+5] = rng.uniform(low=0.4,high=0.6)
            guess[beg+6] = Fguess = F*(1.0 + rng.uniform(low=-0.05, high=0.05))

        return guess


    def _get_mof_prior(self, coord_list):
        """
        prior for N objects.  The priors are the same for
        structural parameters, the only difference being the
        centers
        """

        nobj=len(coord_list)

        cen_priors=[]

        cen_sigma=1.0 #pixels for now
        for coords in coord_list:
            row,col=coords
            p=ngmix.priors.CenPrior(
                row,
                col,
                cen_sigma, cen_sigma,
                rng=rng,
            )
            cen_priors.append(p)

        g_prior=ngmix.priors.GPriorBA(
            0.2,
            rng=rng,
        )
        T_prior = ngmix.priors.TwoSidedErf(
            -1.0, 0.1, 1.0e6, 1.0e5,
            rng=rng,
        )

        fracdev_prior = ngmix.priors.Normal(0.0, 0.1, rng=rng)

        F_prior = ngmix.priors.TwoSidedErf(
            -100.0, 1.0, 1.0e9, 1.0e8,
            rng=rng,
        )

        return minimof.priors.PriorBDFSepMulti(
            cen_priors,
            g_prior,
            T_prior,
            fracdev_prior,
            F_prior,
        )

    def get_mof_model(self):
        import images
        obs, coords = self._get_mof_obs()

        prior=self._get_mof_prior(coords)

        nobj=len(coords)
        mm = minimof.MOF(obs, "bdf", nobj, prior=prior)

        for i in range(2):
            guess=self._get_mof_guess(coords)
            ngmix.print_pars(guess, front="    guess:")
            mm.go(guess)
            res=mm.get_result()
            if res['flags']==0:
                break

        #mim=mm.make_image() 
        #images.compare_images(
        #    obs.image, mim,
        #)
        #images.multiview(obs.psf.image)
        #images.multiview(obs.image)
        #stop


        if res['flags'] == 0:
            ngmix.print_pars(res['pars'], front="      fit:")
            center_obs = mm.make_corrected_obs(0)
            center_obs.noise = obs.noise
        else:
            center_obs=None
        
        return center_obs


    def _rob_deblend(self,im,model,mod1,mod2,dims):
        C = np.zeros((dims[0],dims[1],2))
        W = np.zeros((dims[0],dims[1],2))
        I = im
        w = np.array([model[0,:,:],model[0,:,:]])
        T = np.array([mod1[0,:,:],mod2[0,:,:]])
        mod_sum = np.zeros(dims)
        for r in range(2):
            mod_sum += w[r]*T[r] 
        zeros = np.where(mod_sum == 0.)
        mod_sum[zeros] += 0.000001
        for r in range(2):
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
        #im,psf_im,model,mod1,cen_mod,coords,dx1,dy1,noise = Mod._get_scar_model()
        im,psf_im,model,mod1,mod2,cen_mod,neigh_mod,coords,dx1,dy1,noise = Mod._get_scar_model()
        cen_shape = cen_mod.shape
        im_shape = np.shape(im)

        #check if model dimensions are ever larger than image dims
        #if larger, trim the dimensions
        #This is sometimes an issue for low edge flux test
        if cen_shape[1] > im_shape[1]:
            cen_shape = (1,im_shape[1],cen_shape[2])
        if cen_shape[2] > im_shape[2]:
            cen_shape = (1,cen_shape[1],im_shape[2])

        coord1 = coords[0]
        dims = [np.shape(im)[1],np.shape(im)[2]]
        C,W = Mod._rob_deblend(im,model,mod1,mod2,dims)
        Cnoise,Wnoise = Mod._rob_deblend(noise,model,mod1,mod2,dims)

        half1 = cen_shape[1]/2.
        half2 = cen_shape[2]/2.

        beg1 = int(coord1[0]-half1+1)
        end1 = int(coord1[0]+half1+1)

        beg2 = int(coord1[1]-half2+1)
        end2 = int(coord1[1]+half2+1)

        #metacal needs symmetric image
        if cen_shape[1] != cen_shape[2]:
            cen_obj = np.zeros((max(cen_shape),max(cen_shape)))
            weights = np.zeros((max(cen_shape),max(cen_shape)))
            mod_noise = np.zeros((max(cen_shape),max(cen_shape)))
            
            cen_obj[0:cen_shape[1],0:cen_shape[2]] = C[beg1:end1,beg2:end2,0]
            weights[0:cen_shape[1],0:cen_shape[2]] = W[beg1:end1,beg2:end2,0]
            mod_noise[0:cen_shape[1],0:cen_shape[2]] = Cnoise[beg1:end1,beg2:end2,0]
            shape = C[beg1:end1,beg2:end2,0].shape
            new_coords = (dx1+(shape[1]-1.0)/2.0,dy1+(shape[0]-1.0)/2.0)

        else:
            cen_obj = C[beg1:end1,beg2:end2,0]
            weights = W[beg1:end1,beg2:end2,0]
            mod_noise = Cnoise[beg1:end1,beg2:end2,0]
            new_coords = (dx1+(cen_obj.shape[1]-1.0)/2.0,dy1+(cen_obj.shape[0]-1.0)/2.0)

        
        cen_obj_w_noise = Mod._readd_noise(cen_obj,weights)
        shape = np.shape(cen_obj_w_noise)
        #cen_obj_w_noise += noise[0:shape[0],0:shape[1]]
        noise_w_noise = Mod._readd_noise(mod_noise,weights)
        tot_noise = noise_w_noise#noise[0:shape[0],0:shape[1]]

        dobs = observation(cen_obj_w_noise,Mod['Image']['Bgrms'],
                           new_coords[1],new_coords[0],
                           Mod['Psf']['Bgrms_psf'],psf_im)
        dobs.noise = tot_noise
    
    elif mode=='mof':
        dobs = Mod.get_mof_model()

    elif mode == 'control':
        im,psf_im,coords,dims,dx1,dy1,noise = Mod.__call__()
        output['dims'][j] = np.array(dims)
        
        dobs = observation(im[0],Mod['Image']['Bgrms'],coords[0][0],
                           coords[0][1],Mod['Psf']['Bgrms_psf'],psf_im)
    
        #dobs.noise = noise
    return dobs#,cen_obj,cen_obj_w_noise

def do_metacal(psf_model,gal_model,max_pars,psf_Tguess,prior,
               ntry,metacal_pars,dobs):

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

    return res


dt = [
    ('flags','i4'),
    ('pars','f8',6),
    ('pars_1p','f8',6),
    ('pars_1m','f8',6),
    ('pars_2p','f8',6),
    ('pars_2m','f8',6),
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
output = np.zeros(ntrial, dtype=dt)

for j in range(ntrial):
    print(j)
    try:
        dobs = norm_test()
        if dobs is None:
            output['flags'][j] = 2
        else:
            res = do_metacal(
                psf_model,gal_model,max_pars,
                psf_Tguess,prior,ntry,
                metacal_pars,dobs,
            )

            output['flags'][j] = res['mcal_flags']
            if res['mcal_flags'] == 0:
                output['pars'][j] = res['noshear']['pars']
                output['pars_1p'][j] = res['1p']['pars']
                output['pars_1m'][j] = res['1m']['pars']
                output['pars_2p'][j] = res['2p']['pars']
                output['pars_2m'][j] = res['2m']['pars']

    except (np.linalg.linalg.LinAlgError,ValueError,ngmix.gexceptions.BootGalFailure):
        print("error")
        output['flags'][j] = 2

fitsio.write(outfile_name, output, clobber=True)
