import argparse
import sep
import ngmix
import numpy as np
import scarlet
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

parser = argparse.ArgumentParser()
parser.add_argument("lam",help="lambda to be used for nosie",type=int)
args = parser.parse_args()

def make_image(lam,flux1,flux2,dims):
    """a quick example with two objects convolved by a point spread function """
    # image center is (dims-1)/2
    # [rowcenter, colcenter, g1, g2, T, flux]
    # for a gaussian, T = 2*sigma**2
    psf_pars = [0.0, 0.0, 0.0, 0.0, 4.0, 1.0]

    pars1 = [3.0, 5.0, 0.1, 0.2, 8.0, flux1]
    pars2 = [-4.0, -3.0, -0.05, 0.1, 4.0, flux2]
    
    coord1 = (pars1[0]+(dims[0]-1.)/2.,pars1[1]+(dims[1]-1.)/2.)
    coord2 = (pars2[0]+(dims[0]-1.)/2.,pars2[1]+(dims[1]-1.)/2.)

    psf = ngmix.GMixModel(psf_pars, "gauss")
    gmix1_0 = ngmix.GMixModel(pars1, "gauss")
    gmix2_0 = ngmix.GMixModel(pars2, "exp")

    gmix1 = gmix1_0.convolve(psf)
    gmix2 = gmix2_0.convolve(psf)

    im1 = gmix1.make_image(dims)
    im2 = gmix2.make_image(dims)

    im = im1 + im2
    #add noise
    s = np.random.poisson(lam, (dims[0], dims[1]))
    im += s

    m = im.mean()
    s = im.std()
    err = np.sqrt(lam)

    print("image mean:",m)
    print("image standard dev:",s)

    # add extra dimension for scarlet
    im = im.reshape( (1, dims[0], dims[1]) )
    return im,m,s,err,coord1,coord2

def makeCatalog(img,error):
    detect = img.mean(axis=0)
    #bkg = sep.Background(img.mean(axis=0))
    #bkg = sep.Background(img[0,:,:])
    catalog = sep.extract(detect, 1.5,err=error)# err=bkg.globalrms)
    bg_rms = np.array([error for i in img])#np.array([sep.Background(i).globalrms for i in img])
    return catalog, bg_rms

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
    
    blend = scarlet.Blend(sources, img, bg_rms=bg_rms)
    
    # if you have per-pixel weights:
    #weights = np.empty((B,Ny,Nx))
    #blend = scarlet.Blend(sources, img, bg_rms=bg_rms, weights=weights)

    # if you have per-band PSF kernel images:
    # Note: These need to be difference kernels to a common minimum
    #pdiff = [PSF[b] for b in range(B)]
    #psf = scarlet.transformations.GammaOp(shape, pdiff)
    #blend = scarlet.Blend(sources, img, bg_rms=bg_rms)#, psf=psf)
    # run the fitter for 200 steps (or until convergence)

    blend.fit(200)#, e_rel=1e-1)

    # render the multi-band model: has same shape as img
    model = blend.get_model()
    
    return model

lam = [args.lam]
flux1 = 6000
flux2 = 6000
dims=[50,50]
num_realiz = 10
j = 0
while j < num_realiz:
    for i in lam:
        print("num, lambda:",j,i)
        img,m,s,err,coord1,coord2 = make_image(i,flux1,flux2,dims)
        coords = [coord1,coord2]
        B,Ny,Nx = img.shape
        catalog,bg_rms = makeCatalog(img,err)
        model = make_model(img,bg_rms,B,coords)

        print("coord1,coord2:",coord1,coord2)
    
        plt.figure(figsize=(15,5))
        plt.plot(1,3,3)
    
        plt.subplot(131)
        plt.imshow(img[0,:,:],interpolation='nearest', cmap='gray', vmin=np.min(img[0,:,:]), vmax = np.max(img[0,:,:]))
        plt.colorbar();
        plt.title("Original Image (lambda ="+str(i)+")")

        plt.subplot(132)
        plt.imshow(model[0,:,:],interpolation='nearest', cmap='gray', vmin=np.min(model[0,:,:]),vmax=np.max(model[0,:,:]))
        plt.title("Model (lambda = " + str(i) + ")")
        plt.colorbar();

        diff = img[0,:,:] - model[0,:,:]
        plt.subplot(133)
        plt.imshow(diff,interpolation='nearest', cmap='gray',vmin = np.min(diff),vmax = np.max(diff))
        plt.colorbar();
        plt.title("Original - Model (lambda = "+str(i)+")")
    
        plt.tight_layout()
        #plt.savefig("lam"+str(i)+"_"+str(j)+".png")
        #plt.show()
        plt.clf()
        
        j = j +1

"""
    objects = coords
    fig, ax = plt.subplots()
    m, s = np.mean(model[0,:,:]), np.std(model[0,:,:])
    implot = ax.imshow(model[0,:,:],interpolation='nearest', cmap='gray')#, vmin=m-\s, vmax=m+s, origin='lower')
    cbar = fig.colorbar(implot, ticks=[np.min(model[0,:,:]),0,np.max(model[0,:,:])])

    for i in range(len(objects)):
        plt.scatter(coords[i], coords[i])
        e = Ellipse(xy=(coords[i], coords[i]),
                width=6*coords[i],
                height=6*coords[i],
                angle=catalog['theta'][i] * 180. / np.pi)
        e.set_facecolor('none')
        e.set_edgecolor('red')
        ax.add_artist(e)
    plt.title("Model with Positions (lam = "+str(i)+")")
    plt.savefig("mod_w_pos_lam"+str(i)+".png")
    plt.show()
    plt.clf()
"""
