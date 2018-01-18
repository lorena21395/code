import fitsio
import galsim
import argparse
import sep
import ngmix
import numpy as np
import scarlet
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


def make_image(gal1_flux,gal2_flux,gal1_hlr,gal2_hlr,psf_hlr,dims,scale,bg_rms,bg_rms_psf):
    """a quick example with two objects convolved by a point spread function """
    psf = galsim.Gaussian(half_light_radius = psf_hlr)
    gal1 = galsim.Gaussian(half_light_radius = gal1_hlr, flux=gal1_flux)
    gal2 = galsim.Exponential(half_light_radius = gal2_hlr, flux=gal2_flux)
    dx1,dy1 = 0.0,0.0
    theta = 2*np.pi*np.random.random()
    dx2,dy2 = 12.0*np.cos(theta),12.0*np.sin(theta)
    coord1 = (dx1+(dims[0]-1.)/2.,dy1+(dims[1]-1.)/2.)
    coord2 = (dx2+(dims[0]-1.)/2.,dy2+(dims[1]-1.)/2.)
    coords = [coord1,coord2]
    gal1 = gal1.shift(dx=dx1, dy=dy1)
    gal2 = gal2.shift(dx=dx2, dy=dy2)
    gals = [gal1, gal2]
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
    im = im.reshape( (1, dims[1], dims[0]) )
    return im,coords,psf_im


def make_model(img,bg_rms,B,coords):
    #constraints on morphology:
    # "S": symmetry
    # "m": monotonicity (with neighbor pixel weighting)
    # "+": non-negativity
    #if len(catalog) <2:
    #    print("cat wrong")
    #    catalog = [catalog[0],catalog[0]]
        #print(catalog)
    constraints = {"S": None, "m": {'use_nearest': False}, "+": None}

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
    return model,mod2

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
]

ntrial = 1
output = np.zeros(ntrial, dtype=dt)

scale=1.0
psf_hlr = 1.7
gal1_hlr = 1.7
gal2_hlr = 3.4
gal1_flux = 6000.0
gal2_flux = 8000.0
dims = [50,50]
bg_rms = 100
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

for j in range(ntrial):
    print(j)
    try:
        img,coords,psf_im = make_image(gal1_flux,gal2_flux,gal1_hlr,gal2_hlr,psf_hlr,dims,scale,bg_rms,bg_rms_psf)
        coord1,coord2 = coords[0],coords[1]
        B,Ny,Nx = img.shape
        model,mod2 = make_model(img,bg_rms,B,coords)
        cen_obj = img[0,:,:]-mod2[0,:,:]

        #create container
        obs = observation(cen_obj,bg_rms,int(coord1[0]),int(coord1[1]),bg_rms_psf,psf_im)
        #fit
        boot = ngmix.bootstrap.MaxMetacalBootstrapper(obs)
        boot.fit_psfs(
                psf_model,
                psf_Tguess,
                ntry=ntry,
                fit_pars=max_pars['lm_pars'],
                skip_already_done=False,
            )
        boot.fit_metacal(psf_model,gal_model,max_pars,psf_Tguess,prior=prior,ntry=ntry,
                         metacal_pars=metacal_pars,)
        res = boot.get_metacal_result()
        print("flags:",res['mcal_flags'])
        output['flags'][j] = res['mcal_flags']
        output['pars'][j] = res['noshear']['pars']
        output['pars_1p'][j] = res['1p']['pars']
        output['pars_1m'][j] = res['1m']['pars']
        output['pars_2p'][j] = res['2p']['pars']
        output['pars_2m'][j] = res['2m']['pars']
        
        #print("obects found:",len(model))
        mod_m = np.mean(model)
        mod_s = np.std(model)
        #print("model mean:",mod_m)
        #print("model stand dev:",mod_s)
        #print("coord1,coord2:",coord1,coord2)
        #print("Neighbor mean:",np.mean(mod2))
        
        plt.figure(figsize=(11,8))
        plt.plot(2,3,5)
            
        plt.subplot(231)
        plt.imshow(img[0,:,:],interpolation='nearest', cmap='gray', vmin=np.min(img[0,:,:]), vmax = np.max(img[0,:,:]))
        plt.colorbar();
        plt.title("Original Image (lambda ="+str(i)+")")

        plt.subplot(232)
        plt.imshow(model[0,:,:],interpolation='nearest', cmap='gray', vmin=np.min(model[0,:,:]),vmax=np.max(model[0,:,:]))
        plt.title("Model (lambda = " + str(i) + ")")
        plt.colorbar();

        diff = img[0,:,:] - model[0,:,:]
        plt.subplot(233)
        plt.imshow(diff,interpolation='nearest', cmap='gray',vmin = np.min(diff),vmax = np.max(diff))
        plt.colorbar();
        plt.title("Original - Model (lambda = "+str(i)+")")

        plt.subplot(234)
        plt.imshow(mod2[0,:,:],interpolation='nearest', cmap='gray',vmin = np.min(diff),vmax= np.max(diff))
        plt.colorbar()
        plt.title("Model of Neighbor")

        plt.subplot(235)
        plt.imshow(cen_obj,interpolation='nearest', cmap='gray',vmin = np.min(diff),vmax = np.max(diff))
        plt.colorbar()
        plt.title("Original - Model of Neighbor")

        plt.tight_layout()
        #plt.savefig("/Users/lorena/git/test/blender/new/lam"+str(i)+"_"+str(j)+".png")
        plt.show()
        #plt.clf()
    except (ValueError, np.linalg.linalg.LinAlgError):
        output['flags'][j] = 1
        print("flags: 1")

        
filename = 'data2.fits'
fitsio.write(filename, output, clobber=True)
