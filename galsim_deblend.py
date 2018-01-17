import galsim
import argparse
import sep
import ngmix
import numpy as np
import scarlet
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

parser = argparse.ArgumentParser()
parser.add_argument("lam",help="lambda to be used for nosie",type=float)
args = parser.parse_args()

def make_image(lam,gal1_flux,gal2_flux,gal1_hlr,gal2_hlr,psf_hlr,dims,scale):
    """a quick example with two objects convolved by a point spread function """
    psf = galsim.Gaussian(half_light_radius = psf_hlr)
    gal1 = galsim.Gaussian(half_light_radius = gal1_hlr, flux=gal1_flux)
    gal2 = galsim.Exponential(half_light_radius = gal2_hlr, flux=gal2_flux)
    dx1,dy1 = 0.0,0.0
    theta = 2*np.pi*np.random.random()
    dx2,dy2 = 4.0*np.cos(theta),4.0*np.sin(theta)
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
    #psf_gsim = psf.drawImage(scale=scale)

    #add noise
    noise = galsim.GaussianNoise(sigma=np.sqrt(lam)[0])
    #print("noise:",type(noise))
    #noise = np.random.normal(scale=scale)
    #print("noise",type(noise))
    gsim.addNoise(noise)

    # galsim numpy array is in the .array attrubute
    im = gsim.array
    m = im.mean()
    s = im.std()
    err = np.sqrt(lam)
    print("image err:", err)
    print("image mean:",m)
    print("image standard dev:",s)

    # add extra dimension for scarlet
    im = im.reshape( (1, dims[1], dims[0]) )
    return im,m,s,err[0],coords

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
    sources2 = [sources[1],sources[1]]
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
    mod2 = sources[1].get_model()
    return model,mod2

lam = [args.lam]
scale=1.0
psf_hlr = 1.7
gal1_hlr = 1.7
gal2_hlr = 3.4
gal1_flux = 6000.0
gal2_flux = 8000.0
dims = [50,50]
num_realiz = 10
j = 0

while j < num_realiz:
    for i in lam:
        print("num, lambda:",j,i)
        img,m,s,err,coords = make_image(lam,gal1_flux,gal2_flux,gal1_hlr,gal2_hlr,psf_hlr,dims,scale)
        coord1,coord2 = coords[0],coords[1]
        B,Ny,Nx = img.shape
        catalog,bg_rms = makeCatalog(img,err)
        print("objects found:",len(catalog))
        model,mod2 = make_model(img,bg_rms,B,coords)

        mod_m = np.mean(model)
        mod_s = np.std(model)
        print("model mean:",mod_m)
        print("model stand dev:",mod_s)
        print("coord1,coord2:",coord1,coord2)
    
        print("Neighbor mean:",np.mean(mod2))
        
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
        wall = np.zeros((1,50,50))
        x = coord2[0]
        y = coord2[1]
        wall[0,int(y-mod2.shape[1]/2):int(y+mod2.shape[1]/2),int(x-mod2.shape[2]/2):int(x+mod2.shape[2]/2)] = mod2
        plt.imshow(wall[0,:,:],interpolation='nearest', cmap='gray',vmin = np.min(diff),vmax= np.max(diff))
        plt.colorbar()
        plt.title("Model of Neighbor")

        plt.subplot(235)
        im_minus_mod2 = img[0,:,:]-wall[0,:,:]
        plt.imshow(im_minus_mod2,interpolation='nearest', cmap='gray',vmin = np.min(diff),vmax = np.max(diff))
        plt.colorbar()
        plt.title("Original - Model of Neighbor")

        plt.tight_layout()
        #plt.savefig("/Users/lorena/git/test/blender/new/lam"+str(i)+"_"+str(j)+".png")
        plt.show()
        #plt.clf()
        
        j = j +1
