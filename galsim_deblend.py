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

def make_image(gal1_flux,gal2_flux,gal1_hlr,gal2_hlr,psf_hlr,dims,scale,bg_rms,bg_rms_psf,seed):
    """a quick example with two objects convolved by a point spread function """
    psf = galsim.Gaussian(half_light_radius = psf_hlr)
    gal1 = galsim.Gaussian(half_light_radius = gal1_hlr, flux=gal1_flux)
    gal2 = galsim.Exponential(half_light_radius = gal2_hlr, flux=gal2_flux)
    dx1,dy1 = 0.0,0.0
    theta = 2*np.pi*np.random.random()
    dx2,dy2 = 12.0*np.cos(theta),12.0*np.sin(theta)
    coord1 = (dy1+(dims[0]-1.)/2.,dx1+(dims[1]-1.)/2.)
    coord2 = (dy2+(dims[0]-1.)/2.,dx2+(dims[1]-1.)/2.)
    coords = [coord1,coord2]
    gal1 = gal1.shift(dx=dx1, dy=dy1)
    gal2 = gal2.shift(dx=dx2, dy=dy2)
    if mode == 'scarlet' or mode == 'minimof':
        gals = [gal1, gal2]
        objs = galsim.Add(gals)
    elif mode == 'control':
        gals = [gal1]
        objs = galsim.Add(gals)
    #Add shear
    shear1, shear2 = 0.02, 0.00
    objs  = objs.shear(g1=shear1, g2=shear2)    

    objs = galsim.Convolve(objs, psf)
    gsim = objs.drawImage(nx=dims[1], ny=dims[0], scale=scale)
    psf_gsim = psf.drawImage(scale=scale)
    psf_im = psf_gsim.array

    noise_psf = np.random.normal(scale=bg_rms_psf,size=psf_im.shape)
    psf_im += noise_psf
    
    # galsim numpy array is in the .array attrubute
    im = gsim.array
    
    #np.savez('image_file.npz',image=im)
    #imsave('linalg_image.png',im)
    #add noise
    noise = np.random.normal(scale=bg_rms,size=(dims[0],dims[1])) 
    im += noise
    
    #m = im.mean()
    #s = im.std()
    #err = noise
    #print("image err:", err)
    #print("image mean:",m)
    #print("image standard dev:",s)

    # add extra dimension for scarlet
    if mode == 'scarlet' or mode == 'control':
        im = im.reshape( (1, dims[1], dims[0]) )
    return im,coords,psf_im


def make_model(img,bg_rms,B,coords):
    #constraints on morphology:
    # "S": symmetry
    # "m": monotonicity (with neighbor pixel weighting)
    # "+": non-negativity
    constraints = {"S": None, "m": {'use_nearest': False}, "+": None}
    #constraints['l0'] = bg_rms
    # initial size of the box around each object
    # will be adjusted as needed
    shape = (B, 15, 15)
    sources = [scarlet.Source(coord, shape, constraints=constraints) for coord in coords]
    blend = scarlet.Blend(sources, img, bg_rms=[bg_rms])
    # if you have per-pixel weights:
    #weights = np.empty((B,Ny,Nx))
    #blend = scarlet.Blend(sources, img, bg_rms=bg_rms, weights=weights)

    # if you have per-band PSF kernel images:
    # Note: These need to be difference kernels to a common minimum
    #pdiff = [PSF[b] for b in range(B)]
    #psf = scarlet.transformations.GammaOp(shape, pdiff)
    #blend = scarlet.Blend(sources, img, bg_rms=bg_rms, psf=psf)
    # run the fitter for 200 steps (or until convergence)

    blend.fit(200)#, e_rel=1e-1)
    # render the multi-band model: has same shape as img
    model = blend.get_model()
    mod2 = blend.get_model(m=1)
    cen_mod = sources[0].get_model()
    neigh_mod = sources[1].get_model()
    return model,mod2,cen_mod,neigh_mod

def observation(image,sigma,row,col,psf_sigma,psf_im):
    # sigma is the standard deviation of the noise
    weight = image*0 + 1.0/sigma**2

    # row,col are the center of the object in the model
    # image, so we need to figure that out somehow
    jacob = ngmix.UnitJacobian(row=row, col=col)

    ############Done above already, modify this
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
dt = [
    ('flags','i4'),
    ('pars','f8',6),
    ('pars_1p','f8',6),
    ('pars_1m','f8',6),
    ('pars_2p','f8',6),
    ('pars_2m','f8',6),
    ('noise_std','f8'),
]


output = np.zeros(ntrial, dtype=dt)

scale=1.0
psf_hlr = 1.7
gal1_hlr = 1.7
gal2_hlr = 3.4
gal1_flux = 6000.0
gal2_flux = 8000.0
dims = [50,50]
bg_rms = 10
bg_rms_psf = 0.0001
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
#allconf=yaml.load(open('/astro/u/esheldon/git/nsim/config/run-nbr01-mcal-02.yaml'))
#config = allconf['mof']
#config['psf_pars'] = {'model':'gauss','ntry':2}
subtr_reg_stds = []

for j in range(ntrial):
    print(j)
    try:
        img,coords,psf_im = make_image(gal1_flux,gal2_flux,gal1_hlr,gal2_hlr,psf_hlr,dims,scale,bg_rms,bg_rms_psf,seed)
        coord1,coord2 = coords[0],coords[1]
        if mode == 'minimof':
            print('mini')
            allobs = []
            for coord in coords:
                row,col = coord
                #row,col = int(row),int(col)
                obs = observation(img,bg_rms,row,col,bg_rms_psf,psf_im)
                this_jacob = obs.jacobian.copy()
                this_jacob.set_cen(row=row, col=col)
                this_obs = ngmix.Observation(
                    obs.image.copy(),
                    weight=obs.weight.copy(),
                    jacobian=this_jacob,
                    psf=obs.psf,
                    )
                allobs.append(this_obs)
            mm = minimof.MiniMOF(config, allobs, rng = rng)
            mm.go()
            res=mm.get_result()
            if not res['converged']:
                output['flags'][j] = 2
            else:
                dobs = mm.get_corrected_obs(0)
        
        elif mode == 'scarlet':
            B,Ny,Nx = img.shape
            model,mod2,cen_mod,neigh_mod = make_model(img,bg_rms,B,coords)
            
            #isolate central object
            cen_obj = img[0,:,:]-mod2[0,:,:]
            
            #subtract full model from original image
            orig_minus_model = img[0,:,:]-model[0,:,:]
            #find shape of neighbor object
            neigh_shape = neigh_mod.shape
            #identify region in remainder image associated with neighbor
            region = orig_minus_model[int(coord2[0]-(neigh_shape[1]/2)-1):int(coord2[0]+(neigh_shape[1]/2)-1),int(coord2[1]-(neigh_shape[2]/2)-1):int(coord2[1]+(neigh_shape[2]/2)-1)]
            
            """
            plt.imshow(region[:,:])
            plt.savefig("fig_1.png")

            #subtract model from original image
            orig_minus_model = img[0,:,:]-model[0,:,:]
            
            #find shape of central object model
            cen_obj_shape = cen_mod.shape

            #identify region where central object originally was
            region = orig_minus_model[int(coord1[0]-(cen_obj_shape[1]/2)-1):int(coord1[0]+(cen_obj_shape[1]/2)-1),int(coord1[1]-(cen_obj_shape[2]/2)-1):int(coord1[1]+(cen_obj_shape[2]/2)-1)]
            """

            #calculate std of region
            reg_std = np.std(region)
            subtr_reg_stds.append(reg_std)

            #calculate extra noise
            extra_noise = np.sqrt(np.abs(bg_rms**2 - np.var(region)))
            
            """
            #check if central object is symmetric, if not, make it
            if cen_obj_shape[1] != cen_obj_shape[2]:
                print('asym')
                sym_cen_mod = np.zeros((max(cen_obj_shape),max(cen_obj_shape)))
                sym_cen_mod[0:cen_obj_shape[1],0:cen_obj_shape[2]]=cen_mod[0,:,:]
            elif cen_obj_shape[1] == cen_obj_shape[2]:
                sym_cen_mod = cen_mod[0]
            
            """

            plt.imshow(cen_obj,interpolation='nearest', cmap='gray', vmin=np.min(cen_obj), vmax = np.max(cen_obj))                                        
            plt.colorbar();
            plt.savefig('fig_before.png')
            #print("extra noise: ",extra_noise)
            #print("mean before: ",np.mean(cen_obj))
            cen_obj[int(coord2[0]-(neigh_shape[1]/2)-1):int(coord2[0]+(neigh_shape[1]/2)-1),int(coord2[1]-(neigh_shape[2]/2)-1):int(coord2[1]+(neigh_shape[2]/2)-1)] += rng.normal(scale=extra_noise)
            plt.clf()
            plt.imshow(cen_obj,interpolation='nearest', cmap='gray', vmin=np.min(cen_obj), vmax = np.max(cen_obj))
            plt.colorbar();
            plt.savefig("fig_after.png")
            plt.clf()
#print("mean after: ",np.mean(cen_obj))
            #new_coords = (sym_cen_mod.shape[0]/2-1,sym_cen_mod.shape[1]/2-1)
            dobs = observation(cen_obj,bg_rms,coord1[0],coord1[1],bg_rms_psf,psf_im)
        

        elif mode == 'control':
            cen_obj = img[0,:,:]
            dobs = observation(cen_obj,bg_rms,coord1[0],coord1[1],bg_rms_psf,psf_im)

        #Create container
        #obs = observation(cen_obj,bg_rms,int(coord1[0]),int(coord1[1]),bg_rms_psf,psf_im)
        #fit
        boot = ngmix.bootstrap.MaxMetacalBootstrapper(dobs)
        boot.fit_psfs(
                psf_model,
                psf_Tguess,
                ntry=ntry,
                fit_pars=max_pars['lm_pars'],
                skip_already_done=False,
            )
        boot.fit_metacal(psf_model,gal_model,max_pars,psf_Tguess,prior=prior,ntry=ntry,metacal_pars=metacal_pars,)
        res = boot.get_metacal_result()
        #print("flags:",res['mcal_flags'])
        output['flags'][j] = res['mcal_flags']
        output['pars'][j] = res['noshear']['pars']
        output['pars_1p'][j] = res['1p']['pars']
        output['pars_1m'][j] = res['1m']['pars']
        output['pars_2p'][j] = res['2p']['pars']
        output['pars_2m'][j] = res['2m']['pars']
        output['noise_std'][j] = reg_std
        """
        #print("obects found:",len(model))
        #mod_m = np.mean(model)
        #mod_s = np.std(model)
        #print("model mean:",mod_m)
        #print("model stand dev:",mod_s)
        #print("coord1,coord2:",coord1,coord2)
        #print("Neighbor mean:",np.mean(mod2))
    
        plt.figure(figsize=(11,8))
        plt.plot(2,3,5)
            
        plt.subplot(231)
        plt.imshow(img[0,:,:],interpolation='nearest', cmap='gray', vmin=np.min(img[0,:,:]), vmax = np.max(img[0,:,:]))
        plt.colorbar();
        plt.title("Original Image")
        
        plt.subplot(232)
        plt.imshow(model[0,:,:],interpolation='nearest', cmap='gray', vmin=np.min(model[0,:,:]),vmax=np.max(model[0,:,:]))
        plt.title("Model")
        plt.colorbar();
        
        diff = img[0,:,:] - model[0,:,:]
        plt.subplot(233)
        plt.imshow(diff,interpolation='nearest', cmap='gray',vmin = np.min(diff),vmax = np.max(diff))
        plt.colorbar();
        plt.title("Original - Model")

        plt.subplot(234)
        plt.imshow(mod2[0,:,:],interpolation='nearest', cmap='gray',vmin = np.min(diff),vmax= np.max(diff))
        plt.colorbar()
        plt.title("Model of Neighbor")

        plt.subplot(235)
        plt.imshow(cen_obj,interpolation='nearest', cmap='gray',vmin = np.min(diff),vmax = np.max(diff))
        plt.colorbar()
        plt.title("Original - Model of Neighbor")

        plt.tight_layout()
        plt.savefig("/gpfs01/astro/workarea/lmezini/scarlet-tests/test_"+str(j)+".png")
        #plt.show()
        #plt.clf()
        """
    except (np.linalg.linalg.LinAlgError,ValueError):
        output['flags'][j] = 2
        #print("flags: 1")


noise_mean = np.mean(subtr_reg_stds)
std = np.std(subtr_reg_stds)
avg_std_err = std/np.sqrt(len(subtr_reg_stds))

print(noise_mean,avg_std_err)
fitsio.write(outfile_name, output, clobber=True)
