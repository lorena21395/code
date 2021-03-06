#!/usr/bin/env python
#####
import time
import yaml
import minimof
import fitsio
import galsim
import argparse
import sep
import ngmix
import esutil as eu
from ngmix.observation import Observation, ObsList, MultiBandObsList
import numpy as np
#import matplotlib.pyplot as plt
#plt.switch_backend('agg')
#import logging
#logger = logging.getLogger('scarlet')
#logger.setLevel(logging.DEBUG)

from numpy.linalg import LinAlgError
from ngmix.gexceptions import BootGalFailure

#from astropy.io import fits

parser = argparse.ArgumentParser()
parser.add_argument("outfile",help="Output file name and path")
parser.add_argument("ntrials",help="Number of trials to be run",type = int)
parser.add_argument("seed",help="Seed for random number generator",type = int)
parser.add_argument("config",help="Configuration for the simulation")

parser.add_argument("--profile",action='store_true',help="run the profiler")
parser.add_argument("--show",action='store_true',help="show images")

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


class Simulation(dict):
    """Simulate galaxy images with noise and psf"""

    def __init__(self,specs, rng):
        self.update(specs)
        self.specs = specs
        self.rng=rng
        self.galsim_rng = galsim.BaseDeviate(self.rng.randint(0,2**30))

    def get_mode(self):
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

    def _add_ellipticity(self, obj, spec):
        if 'g' in spec:
            gsigma=spec['g']['sigma']
            g1,g2 = self.rng.normal(scale=gsigma, size=2)
            obj = obj.shear(g1=g1, g2=g2)

        return obj

    def _split_flux_with_knots(self, spec, flux):
        krng=spec['knots']['fluxfrac']['range']
        knot_frac = self.rng.uniform(
            low=krng[0],
            high=krng[1],
        )
        knot_flux = flux*knot_frac
        remaining_flux = flux*(1.0 - knot_frac)

        return remaining_flux, knot_flux

    def _get_knots(self, spec, hlr, knot_flux):
        return galsim.RandomWalk(
            npoints=spec['knots']['num'],
            half_light_radius=hlr,
            flux=knot_flux,
            rng=self.galsim_rng,
        )

    def _get_bd_model(self, spec):
        rng=self.rng

        flux = spec['Flux']

        fracdev=rng.uniform(low=0.0, high=1.0)

        bulge_flux = fracdev*flux
        disk_flux = (1.0 - fracdev)*flux

        disk_hlr = spec['hlr']
        if 'bulge_size_factor' in spec:
            frange=spec['bulge_size_factor']['range']
            bulge_hlr = disk_hlr*rng.uniform(
                low=frange[0],
                high=frange[1],
            )
        else:
            bulge_hlr=disk_hlr

        if 'knots' in spec:
            disk_flux, knot_flux = self._split_flux_with_knots(spec, disk_flux)

        disk = galsim.Exponential(
            half_light_radius=disk_hlr,
            flux=disk_flux,
        )
        bulge = galsim.DeVaucouleurs(
            half_light_radius=bulge_hlr,
            flux=bulge_flux,
        )

        if 'knots' in spec:
            knots = self._get_knots(spec, disk_hlr, knot_flux)
            disk = galsim.Add([disk, knots])


        disk = self._add_ellipticity(disk, spec)
        bulge = self._add_ellipticity(bulge, spec)


        return galsim.Add([disk, bulge])



    def _get_sersic_model(self, spec):
        rng=self.rng

        nrange=spec['n']['range']
        n=rng.uniform(low=nrange[0], high=nrange[1])

        obj = galsim.Sersic(
            n,
            half_light_radius=spec['hlr'],
            flux=spec['Flux'],
        )

        return obj

    def _get_exp_model(self, spec):
        obj = galsim.Exponential(
            half_light_radius=spec['hlr'],
            flux=spec['Flux'],
        )

        return obj


    def _get_dev_model(self, spec):
        obj = galsim.DeVaucouleurs(
            half_light_radius=spec['hlr'],
            flux=spec['Flux'],
        )

        return obj


    def _get_gauss_model(self, spec):
        obj = galsim.Gaussian(
            half_light_radius=spec['hlr'],
            flux=spec['Flux'],
        )

        return obj

    def _get_pure_knots_model(self, spec):
        return self._get_knots(
            spec,
            spec['hlr'],
            spec['Flux'],
        )

    def _get_single_component_model(self, spec):
        mod=spec['Type']
        if mod=='Sersic':
            obj=self._get_sersic_model(spec)
        elif mod=='Gaussian':
            obj=self._get_gauss_model(spec)
        elif mod=='Exponential':
            obj=self._get_exp_model(spec)
        elif mod=='DeVaucouleurs':
            obj=self._get_dev_model(spec)
        else:
            raise RuntimeError("bad model: '%s'" % mod)


        return obj

    def _get_galaxy_model(self, spec):

        if spec['Type']=='bulge+disk':
            obj=self._get_bd_model(spec)
        elif spec['Type']=='Knots':
            obj=self._get_pure_knots_model(spec)
            obj = self._add_ellipticity(obj, spec)
        else:
            obj=self._get_single_component_model(spec)

            if 'knots' in spec:
                flux=obj.flux
                hlr=obj.half_light_radius

                remaining_flux, knot_flux = \
                        self._split_flux_with_knots(spec, flux)

                knots = self._get_knots(spec, hlr, knot_flux)
                obj = galsim.Add([obj, knots])

            obj = self._add_ellipticity(obj, spec)

        return obj

    def _get_cen_model(self):
        spec=self['Cen']
        return self._get_galaxy_model(spec)

    def _get_nbr_model(self):
        spec=self['Neigh']
        return self._get_galaxy_model(spec)

    def _get_norm_sed(self):
        cen_sed = np.array(self['Cen']['SED'])/np.sum(np.array(self['Cen']['SED']))
        neigh_sed = np.array(self['Neigh']['SED'])/np.sum(np.array(self['Neigh']['SED']))
        return cen_sed, neigh_sed

    def _get_shifts(self):
        rng=self.rng

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
    
        return dx1,dy1,dx2,dy2

    def _get_gals(self):
        dx1,dy1,dx2,dy2 = self._get_shifts()

        Cen = self._get_cen_model()
        Neigh = self._get_nbr_model()
        Cen = Cen.shift(dx=dx1, dy=dy1)
        Neigh = Neigh.shift(dx=dx2, dy=dy2)

        return Cen, Neigh, dx1,dy1,dx2,dy2

    def _get_band_obj(self,Cen,Neigh,cen_sed,neigh_sed):
        mode = self['Mode']
        if mode == 'scarlet' or mode == 'mof':
            Cen = Cen*cen_sed
            Neigh = Neigh*neigh_sed
            gals = [Cen,Neigh]
            objs = galsim.Add(gals)
        elif mode == 'control':
            Cen = Cen*cen_sed
            gals = [Cen]
            objs = galsim.Add(gals)
        return objs
    
    def _get_shear_obj(self,objs):
        shear_spec = self['Shear']
        shear1, shear2 = shear_spec['Shear1'],shear_spec['Shear2']
        objs = objs.shear(g1=shear1, g2=shear2)

        if 'Extra_Shear' in shear_spec:
            if 'shear+' in shear_spec['Extra_Shear']:
                objs = objs.shear(g1=shear_spec['Extra_Shear']['shear+'],g2=0.0)
            elif 'shear-' in shear_spec['Extra_Shear']:
                objs = objs.shear(g1=shear_spec['Extra_Shear']['shear-'],g2=0.0)
            
        return objs

    def _get_noise(self,dims,bg_rms):
        noise = np.random.normal(scale=bg_rms,size=(dims[0],dims[1]))
        
        return noise


    def __call__(self):
        mode = self['Mode']
        bg_rms = self['Image']['Bgrms']

        if 'Nbands' in self.specs:
            nband = self['Nbands']
            cen_sed,neigh_sed = self._get_norm_sed()
        else:
            nband = 1
            cen_sed = [1.]
            neigh_sed = [1.]
        
        bg_rms = bg_rms/np.sqrt(nband)

        if 'dims' in self['Image']:
            ny,nx=self['Image']['dims']
        else:
            ny, nx = None,None

        psf, psf_im = self._get_psf_img()
        Cen,Neigh,dx1,dy1,dx2,dy2 = self._get_gals()
        
        ims = []
        noises = []
        for i in range(nband):
            objs = self._get_band_obj(Cen,Neigh,cen_sed[i],neigh_sed[i])
            objs = self._get_shear_obj(objs)
            objs = galsim.Convolve(objs, psf)
            gsim = objs.drawImage(
                nx=nx,
                ny=ny,
                scale=self['Image']['Scale'],
            )
            im = gsim.array
            dims = np.shape(im)
            noise = self._get_noise(dims,bg_rms)
            #im += noise
            noise = self._get_noise(dims,bg_rms)
            noises.append(noise)
            ims.append(im)
        ims = np.array(ims)

        #if using sep
        """
        sep_im = np.sum(ims,axis=0)
        objects = sep.extract(sep_im, 1.5, err=10.)
        s_coords = [(objects['y'][0],objects['x'][0]),(objects['y'][1],objects['x'][1])]
        """
        """
        #create fits file
        hdu1 = fits.PrimaryHDU()
        hdu1.data = ims
        hdu1.writeto('multi_band_im.fits',overwrite=True)
        """

        #rewrite coords for scarlet
        
        cen =  (np.array(dims) - 1.0)/2.0
        coord1 = (dy1+cen[1],dx1+cen[0])
        coord2 = (dy2+cen[1],dx2+cen[0])
        
        #figure out which sep coord belong to each object
        """
        if np.abs(coord1[0]-s_coords[0][0]) < np.abs(coord2[0]-s_coords[0][0]) and np.abs(coord1[1]-s_coords[0][1]) < np.abs(coord2[1]-s_coords[0][1]):
            s_coord1 = s_coords[0]
            s_coord2 = s_coords[1]
        else:
            s_coord1 = s_coords[1]
            s_coord2 = s_coords[0]
        """
        if mode == 'control':
            coords = [coord1]
            #coords = [s_coord1]
        else:
            coords = [coord1,coord2]
            #coords = [s_coord1,s_coord2]

        return ims,psf_im,coords,dims,dx1,dy1,noises

class NGal_Simulation(dict):
    def __init__(self,specs,seed,rng):
        self.specs = specs
        self.seed = seed
        self.rng = rng

    def get_mode(self):
        mode = 'Ngal_mof'
        
        return mode

    def _get_sep_objs(self,sim):
        objects, segmap = sep.extract(
            sim.image,
            1.5,
            #.5,                                                               
            err=sim['noise_sigma'],
            #minarea=2,                                                         
            deblend_cont=0.0001,
            segmentation_map=True,
        )
        
        return objects, segmap
    
    def _get_coord_list(self,objects):
        coord_list = []
        for i in range(len(objects)):
            coord_list.append([objects['y'][i],objects['x'][i]])

        return coord_list

    def _get_fit_model(self):
        fit_model = self.specs['fit_model']
        
        return fit_model

    def _get_noise_im(self,obs):
        wt = obs.weight
        noise = np.random.normal(size=wt.shape)*np.sqrt(1./wt)
        return noise

    def __call__(self):
        sim = minimof.moftest.Sim(self.specs,self.seed)
        sim.make_obs()
        objects, segmap = self._get_sep_objs(sim)
        noise = self._get_noise_im(sim.obs)
        sim.obs.noise = noise
        fit_model = self._get_fit_model()
        dims = self.specs['dims']
        nfail = 0
        nobj = 0
        nobj = len(objects)
        if nobj == 0:
            print("failed to find any objects")
            nfail += 1
            coord_list = []
        else:
            coord_list = self._get_coord_list(objects)

        return sim,nobj,nfail,objects,coord_list,fit_model,noise,dims

class Model(object):
    
    def __init__(self, sim, show=False):
        self.sim=sim
        self.show=show

    def get_scar_model(self):
        import scarlet
        import scarlet.constraint as sc
        im,psf_im,coords,dims,dx1,dy1,noise = self.sim()
        bg_rms = self.sim['Image']['Bgrms']
        bg_rms = bg_rms/np.sqrt(len(im))
        mode = self.sim['Mode']
        
        #create multiband background im
        bg = np.zeros((len(im)))
        bg += bg_rms
        
        #choose constraints for PSF and sources (optional)
        #psf_constraints = ()#sc.SimpleConstraint())
        #source_constraints = ()#sc.SimpleConstraint(),sc.DirectSymmetryConstraint()) #sc.DirectMonotonicityConstraint(use_nearest=False)
        
        config = scarlet.Config(source_sizes = [25])
        
        #create multiband PSF image
        #psf_dims = np.shape(psf_im)
        #psf_im3d = psf_im.reshape( (1, psf_dims[0], psf_dims[1]) )
        #psfs = np.zeros((len(im), psf_dims[0],psf_dims[1]))
        #psfs += psf_im3d

        #PSF Matching
        #target_psf = scarlet.psf_match.fit_target_psf(psfs, 
        #                                scarlet.psf_match.gaussian)
        #diff_kernels, psf_blend = scarlet.psf_match.build_diff_kernels(psfs,
        #                                    target_psf)
        #if use PSF matching, set psf = diff_kernels below
        sources = [scarlet.ExtendedSource(coord, im, bg_rms = bg,
                psf=None,config=config) for coord in coords]
        #config = scarlet.Config(edge_flux_thresh=0.05)
        blend = scarlet.Blend(sources)
        blend.set_data(im, bg_rms = bg,config=config)
        blend.fit(10000, e_rel=1e-3)
        
        #get full sized scarlet models
        model = blend.get_model()
        mod1 = blend.get_model(0)
        mod2 = blend.get_model(1)
        
        #get postage stamp models
        cen_mod = sources[0].get_model()
        neigh_mod = sources[1].get_model()
        
        #get scarlet coordinates
        cen_cen_pos = blend[0][0].center
        neigh_cen_pos = blend[1][0].center

        return im,psf_im,model,mod1,mod2,cen_mod,neigh_mod,coords,dx1,dy1,noise,cen_cen_pos,neigh_cen_pos
        #return im,psf_im,model,mod1,cen_mod,coords,dx1,dy1,noise 


    def _fit_psf_admom(self, obs):
        Tguess=4.0
        am=ngmix.admom.run_admom(obs, Tguess)
        return am.get_gmix()

    def _get_mof_obs(self):
        im,psf_im,coords,dims,dx1,dy1,noises = self.sim()

        bg_rms = self.sim['Image']['Bgrms']/np.sqrt(len(im))
        bg_rms_psf = self.sim['Psf']['Bgrms_psf']

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

        mb = MultiBandObsList()
        for i in range(len(im)):
            if self.show:
                import images
                tim = im[i]/im[i].max()
                tim = np.log10( tim-tim.min() + 1.0 )

                images.view(tim)
                if 'q'==input('hit a key: (q to quit) '):
                    stop


            weight=im[i]*0 + 1.0/bg_rms**2
            jacobian=ngmix.UnitJacobian(
                row=0,
                col=0,
            )

            obs = ngmix.Observation(
                im[i],
                weight=weight,
                jacobian=jacobian,
                psf=psf_obs,
            )
            obs.noise=noises[i]
            olist = ObsList()
            olist.append(obs)
            mb.append(olist)

        return mb, coords

    def _get_mof_guess(self, coord_list, mb, jacobian = None, objects = None):
        rng=self.sim.rng
        if jacobian == None:
            nbands = len(mb)
        else:
            nbands = 1
        npars_per=6+nbands
        num=len(coord_list)

        #assert num==2,"two objects for now"

        npars_tot = num*npars_per
        guess = np.zeros(npars_tot)

        for i,coords in enumerate(coord_list):
            if jacobian != None:
                row=objects['y'][i]
                col=objects['x'][i]
                v, u = jacobian(row, col)
                scale=jacobian.get_scale()
                pos_range = 0.1*scale
                T=scale**2 * (objects['x2'][i] + objects['y2'][i])
                F=scale**2 * objects['flux'][i]

                beg=i*npars_per
                
                guess[beg+0] = v + rng.uniform(low=-pos_range, high=pos_range)
                guess[beg+1] = u + rng.uniform(low=-pos_range, high=pos_range)
            
            else:
                row,col = coords
                beg=i*npars_per
                T = 10.0
                if i==0:
                    F = self.sim['Cen']['Flux']
                else:
                    F = self.sim['Neigh']['Flux' ]

                # always close guess for center
                guess[beg+0] = row + rng.uniform(low=-0.05, high=0.05)
                guess[beg+1] = col + rng.uniform(low=-0.05, high=0.05)

            # always arbitrary guess for shape
            guess[beg+2] = rng.uniform(low=-0.05, high=0.05)
            guess[beg+3] = rng.uniform(low=-0.05, high=0.05)

            # we could get a better guess from sep
            guess[beg+4] = T*(1.0 + rng.uniform(low=-0.05, high=0.05))

            guess[beg+5] = rng.uniform(low=0.4,high=0.6)
            
            for i in range(1,nbands+1):
                guess[beg+5+i] = F*(1.0 + rng.uniform(low=-0.05, high=0.05))

        return guess

    def _get_mof_prior(self, coord_list, nband, jacobian = None, objects = None):
        """
        prior for N objects.  The priors are the same for
        structural parameters, the only difference being the
        centers
        """
        rng=self.sim.rng

        nobj=len(coord_list)

        cen_priors=[]

        cen_sigma=1.0 #pixels for now
        for i,coords in enumerate(coord_list):
            if jacobian != None:
                row=objects['y'][i]#-1
                col=objects['x'][i]#-1
                v, u = jacobian(row, col)
                cen_sigma = jacobian.get_scale()
                p=ngmix.priors.CenPrior(
                    v,
                    u,
                    cen_sigma, cen_sigma,
                    rng=rng,
                )
                cen_priors.append(p)
            else:
                row, col = coords
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
            [F_prior]*nband,
        )

    def get_mof_ngal_model(self):

        sim,nobj,nfail,objects,coord_list,fit_model,noise,dims = self.sim()

        ####Created simple simulation for testing
        scale = 0.263
        psf = galsim.Gaussian(fwhm = 0.9)
        psf_gsim = psf.drawImage(nx = 25, ny=25, scale = scale)
        psf_im = psf_gsim.array

        psf_err=0.0001
        noise_psf = np.random.normal(scale=psf_err,size=psf_im.shape)
        psf_im += noise_psf
        psf_weight = psf_im*0 + 1.0/psf_err**2
        
        obj = galsim.Gaussian(
            half_light_radius=.5,
            flux=6000.,
        )
        dx1 = np.random.uniform(low=-0.5, high=0.5)
        dy1 = np.random.uniform(low=-0.5, high=0.5)
        obj = obj.shift(dx=dx1*scale, dy=dy1*scale)
        obj = galsim.Convolve(obj, psf)
        gsim = obj.drawImage(
                nx=dims[0],
                ny=dims[1],
                scale=scale,
            )
        im = gsim.array

        err=0.001
        noise=np.random.normal(scale=err,size=im.shape)
        im += noise
        weight = im*0 + 1.0/err**2

        psf_cen =  (np.array(psf_im.shape) - 1.0)/2.0
        psf_jac=ngmix.DiagonalJacobian(scale=scale,row=psf_cen[0],col=psf_cen[1])
        psf_obs = ngmix.Observation(
            psf_im,
            weight=psf_weight,
            jacobian=psf_jac,
        )

        cen =  (np.array(dims) - 1.0)/2.0
        coord1 = (dy1+cen[1],dx1+cen[0])
        jac=ngmix.DiagonalJacobian(
			scale=scale,row=cen[0]+dy1,col=cen[1]+dx1)
        obs = ngmix.Observation(
            im,
            weight=weight,
            jacobian=jac,
            psf=psf_obs,
        )
        noise=np.random.normal(scale=err,size=im.shape)
        obs.noise=noise
        #coord2 = (dy2+cen[1],dx2+cen[0])
        coord_list = [coord1]
            
#        obs = observation(im,0.001,coord1[0],coord1[1],0.0001,psf_im)
 #       obs.noise = noise
  
 #noise = np.random.normal(scale=0.001,size=(dims[0],dims[1]))
        
        return [obs], coord_list

        jac=sim.obs.jacobian
        jac.set_cen(row=objects['y'][0],col=objects['x'][0])
        sim.obs.jacobian=jac
        return [sim.obs], coord_list
        
        ###Orginal start of get ngal model
        prior = self._get_mof_prior(coord_list,1,sim.obs.jacobian,objects)
        fitter = minimof.MOF(
            sim.obs,
            fit_model,
            len(coord_list),
            prior=prior,
        )
        nobs = []
        if len(coord_list) != 0:
            for i in range(2):
                guess = self._get_mof_guess(coord_list, sim.obs,
                                        sim.obs.jacobian, objects)
                fitter.go(guess)
                res = fitter.get_result()
                if res['flags']==0:
                    break
        
            if res['flags'] == 0:
                ngmix.print_pars(res['pars'], front="      fit:")
                for i in range(len(coord_list)):
                    obs = fitter.make_corrected_obs(i, band=None,
                                               obsnum=None, recenter=True)
                    obs[0][0].noise = sim.obs.noise
                    nobs.append(obs)
                
        return nobs, coord_list

    def get_mof_model(self):
        mb, coords = self._get_mof_obs()
        prior=self._get_mof_prior(coords, len(mb))

        nobj=len(coords)
        mm = minimof.MOF(mb, "bdf", nobj, prior=prior)
        for i in range(2):
            guess=self._get_mof_guess(coords,mb)
            ngmix.print_pars(guess, front="    guess:")
            mm.go(guess)
            res=mm.get_result()
            if res['flags']==0:
                break


        if res['flags'] == 0:
            ngmix.print_pars(res['pars'], front="      fit:")
            center_obs = mm.make_corrected_obs(0, band=None, 
                                               obsnum=None, recenter=True)
            for i in range(len(center_obs)):
                center_obs[i][0].noise = mb[i][0].noise

            if self.show:
                import images
                cim=center_obs.image
                s=cim.shape
                tim=np.zeros((s[0], 2*s[1]))

                tim[:, 0:s[1]] = obs.image
                tim[:, s[1]:] = cim
                tim *= 1.0/tim.max()
                tim = np.log10( tim-tim.min() + 1.0 )

                images.view(
                    #tim-tim.min(),
                    tim,
                    title='MOF',
                )
                if 'q'==input('hit a key: (q to quit) '):
                    stop

        else:
            center_obs=None
        
        return center_obs,coords


    def rob_deblend(self,im,model,mod1,mod2,dims,i):
        #SDSS style second step for deblending
        C = np.zeros((dims[0],dims[1],2))
        W = np.zeros((dims[0],dims[1],2))
        I = im[i]
        w = np.array([model[i,:,:],model[i,:,:]]).clip(1.0e-15)
        T = np.array([mod1[i,:,:],mod2[i,:,:]]).clip(1.0e-15)
        mod_sum = np.zeros(dims)
        for r in range(2):
            mod_sum += w[r]*T[r] 
        for r in range(2):
            W[:,:,r] = np.divide(w[r]*T[r],mod_sum)
            C[:,:,r] = I*np.divide(w[r]*T[r],mod_sum)
        return C,W

    def readd_noise(self,cen_obj,W):
        rng=self.sim.rng
        bg_rms = self.sim['Image']['Bgrms']
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

def get_prior(nband, row, col):
    """
    metacal prior
    """
    cen_sigma = 1.0
    cen_prior = ngmix.priors.CenPrior(
        row, col,
        cen_sigma,
        cen_sigma,
    )
    g_prior = ngmix.priors.GPriorBA(0.2)

    T_pars =  [-10.0, 0.03, 1.0e+06, 1.0e+05]
    T_prior = ngmix.priors.TwoSidedErf(*T_pars)

    flux_pars = [-1.0e+04, 1.0, 1.0e+09, 0.25e+08]
    flux_prior = ngmix.priors.TwoSidedErf(*flux_pars)

    flux_prior=[flux_prior]*nband

    prior = ngmix.joint_prior.PriorSimpleSep(
        cen_prior,
        g_prior,
        T_prior,
        flux_prior,
    )
    return prior

def create_mb(mb_mod,bg_rms,coord1,coord2,psf_bg_rms,psf_im,noise):
    #creats multiband object
    mb = MultiBandObsList()
    for i in range(len(mb_mod[:])):
        o = observation(mb_mod[i],bg_rms,coord1,coord2,
                               psf_bg_rms,psf_im)
        o.noise = noise[i]
        olist = ObsList()
        olist.append(o)
        mb.append(olist)

    return mb

def norm_test(args, sim):

    mode = sim.get_mode()

    if mode == 'control':
        im,psf_im,coords,dims,dx1,dy1,noise = sim()
        output['dims'][j] = np.array(dims)
        bg_rms = Mod['Image']['Bgrms']/np.sqrt(len(im))
        dobs = observation(im[0],coords[0][0],
                           coords[0][1],Mod['Psf']['Bgrms_psf'],psf_im)

    else:
 
        mod = Model(sim, show=args.show)
        if mode == 'scarlet':
            im,psf_im,model,mod1,mod2,cen_mod,neigh_mod,coords,dx1,dy1,noise,cen_cen_pos,neigh_cen_pos = mod.get_scar_model()
            bg_rms = sim['Image']['Bgrms']/np.sqrt(len(im))
            dims = [np.shape(im)[1],np.shape(im)[2]]
            #create coord list including scarlet and true coords
            #can later calculate how scarlet coordinates deviate from true
            coord_list = []
            coord_list.append(np.array(coords)[0])
            coord_list.append(np.array(coords)[1])
            coord_list.append(cen_cen_pos)
            coord_list.append(neigh_cen_pos)

            if sim['Steps'] == 'one':
                coord1 = coords[0]
                #isolate central object                                       
                cen_obj = im[:,:,:] - mod2[:,:,:]
                dobs = create_mb(cen_obj,bg_rms,
                                coord1[0],coord1[1],
                                sim['Psf']['Bgrms_psf'],psf_im,noise)

            if sim['Steps'] == 'two':
                cen_shape = cen_mod.shape
                im_shape = np.shape(im)
                #check if model dimensions are ever larger than image dims
                #if larger, trim the dimensions
                #This is sometimes an issue for low edge flux test
                if cen_shape[1] > im_shape[1]:
                    cen_shape = (len(im),im_shape[1],cen_shape[2])
                if cen_shape[2] > im_shape[2]:
                    cen_shape = (len(im),cen_shape[1],im_shape[2])
                coord1 = coords[0]
                half1 = cen_shape[1]/2.
                half2 = cen_shape[2]/2.

                beg1 = int(coord1[0]-half1+1)
                end1 = int(coord1[0]+half1+1)

                beg2 = int(coord1[1]-half2+1)
                end2 = int(coord1[1]+half2+1)

                dims = [np.shape(im)[1],np.shape(im)[2]]

                #metacal needs symmetric image
                if cen_shape[1] != cen_shape[2]:
                    cen_obj = np.zeros((len(im),max(cen_shape),max(cen_shape)))
                    weights = np.zeros((len(im),max(cen_shape),max(cen_shape)))
                    mod_noise = np.zeros((len(im),max(cen_shape),max(cen_shape)))
                    for i in range(len(im)):
                        C,W = mod.rob_deblend(im,model,mod1,mod2,dims,i)
                        Cnoise,Wnoise = mod.rob_deblend(noise,model,mod1,mod2,dims,i)
                        cen_obj[i,0:cen_shape[1],0:cen_shape[2]] = C[beg1:end1,beg2:end2,0]
                        weights[i,0:cen_shape[1],0:cen_shape[2]] = W[beg1:end1,beg2:end2,0]
                        mod_noise[i,0:cen_shape[1],0:cen_shape[2]] = Cnoise[beg1:end1,beg2:end2,0]
                    shape = C[beg1:end1,beg2:end2,0].shape
                    new_coords = (dx1+(shape[1]-1.0)/2.0,dy1+(shape[0]-1.0)/2.0)

                else:
                    cen_obj = np.zeros(cen_shape)
                    weights = np.zeros(cen_shape)
                    mod_noise = np.zeros(cen_shape)
                    for i in range(len(im)):
                        C,W = mod.rob_deblend(im,model,mod1,mod2,dims,i)
                        Cnoise,Wnoise = mod.rob_deblend(noise,model,mod1,mod2,dims,i)
                        cen_obj[i] = C[beg1:end1,beg2:end2,0]
                        weights[i] = W[beg1:end1,beg2:end2,0]
                        mod_noise[i] = Cnoise[beg1:end1,beg2:end2,0]
                    new_coords = (dx1+(cen_obj.shape[1]-1.0)/2.0,dy1+(cen_obj.shape[2]-1.0)/2.0)
                #cen_obj_w_noise = mod.readd_noise(cen_obj,weights)
                #shape = np.shape(cen_obj_w_noise)
                #cen_obj_w_noise += noise[0:shape[0],0:shape[1]]
                #noise_w_noise = mod.readd_noise(mod_noise,weights)
                #tot_noise = noise_w_noise#noise[0:shape[0],0:shape[1]]

                """
                dobs = observation(cen_obj,sim['Image']['Bgrms'],
                               new_coords[1],new_coords[0],
                               sim['Psf']['Bgrms_psf'],psf_im)
                dobs.noise = mod_noise
                """
                dobs = create_mb(cen_obj,bg_rms,
                                 new_coords[1],new_coords[0],
                                sim['Psf']['Bgrms_psf'],psf_im,mod_noise)
        elif mode == 'mof':
            dobs,coord_list = mod.get_mof_model()
        
        elif mode == 'Ngal_mof':
            dobs,coord_list = mod.get_mof_ngal_model()
    
    return dobs, coord_list

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

def get_output(n, nband):
    npar=5+nband
    dt = [
        ('flags','i4'),
        ('pars','f8',npar),
        ('pars_1p','f8',npar),
        ('pars_1m','f8',npar),
        ('pars_2p','f8',npar),
        ('pars_2m','f8',npar),
        ('nobj','i4'),
        ('nobj_found','i4'),
    ]


    return np.zeros(n, dtype=dt)

def main(args):

    with open(args.config) as fobj:
        Sim_specs = yaml.load(fobj)

    
    nband = Sim_specs.get('Nbands',1)
    nobj = Sim_specs.get('nobj',2)
    prior = get_prior(nband,0.0,0.0)

    if nobj == 2:
        output = get_output(args.ntrials, nband)
        np.random.seed(args.seed)
        rng = np.random.RandomState(seed=np.random.randint(0,2**30))
        sim=Simulation(Sim_specs, rng)
    else:
        output = []
        rng = np.random.RandomState(seed=np.random.randint(0,2**30))
        sim = NGal_Simulation(Sim_specs,args.seed,rng)
    tm_deblend=0.0
    tm_metacal=0.0

    for j in range(args.ntrials):
        print(j)
        try:
            tm0=time.time()
            dobs,coord_list = norm_test(args,sim)
            tm_deblend += time.time()-tm0
            if dobs is None:
                if nobj == 2:
                    output['flags'][j] = 2
                else:
                    ngal_output = get_output(1, nband)
                    ngal_output['flags'] = 2
                    output.append(ngal_output)
            else:
                if nobj == 2:
                    tm0=time.time()
                    res = do_metacal(
                        psf_model,gal_model,max_pars,
                        psf_Tguess,prior,ntry,
                        metacal_pars,dobs,
                    )
                    tm_metacal += time.time()-tm0

                    output['flags'][j] = res['mcal_flags']
                    if res['mcal_flags'] == 0:
                        output['pars'][j] = res['noshear']['pars']
                        output['pars_1p'][j] = res['1p']['pars']
                        output['pars_1m'][j] = res['1m']['pars']
                        output['pars_2p'][j] = res['2p']['pars']
                        output['pars_2m'][j] = res['2m']['pars']

                else:
                    ngal_res = []
                    for i in range(len(coord_list)):
                        tm0=time.time()
                        res = do_metacal(
                            psf_model,gal_model,max_pars,
                            psf_Tguess,prior,ntry,
                            metacal_pars,dobs[i],
                        )
                        tm_metacal += time.time()-tm0
                        ngal_output = get_output(1, nband)
                        ngal_output['flags'][0] = res['mcal_flags']
                        if res['mcal_flags'] == 0:
                            ngal_output['pars'][0] = res['noshear']['pars']
                            ngal_output['pars_1p'][0] = res['1p']['pars']
                            ngal_output['pars_1m'][0] = res['1m']['pars']
                            ngal_output['pars_2p'][0] = res['2p']['pars']
                            ngal_output['pars_2m'][0] = res['2m']['pars']
                            ngal_output['nobj'][0] = 1000
                            ngal_output['nobj_found'][0] = 1000
                            ngal_res.append(ngal_output)
                    ngal_res = eu.numpy_util.combine_arrlist(ngal_res)
                    if len(coord_list) != 0:
                        ngal_res['nobj'][-1] = nobj
                        ngal_res['nobj_found'][-1] = len(coord_list)
                    output.append(ngal_res)
        #except (LinAlgError,ValueError,BootGalFailure) as err:
        #except ValueError as err:
        except BootGalFailure as err:
            print(str(err))
            output['flags'][j] = 2
    tm_deblend_per = tm_deblend/args.ntrials
    tm_metacal_per = tm_metacal/args.ntrials

    print("time sim+deblend: %g (%g per object)" % (tm_deblend,tm_deblend_per))
    print("time metacal: %g (%g per object)" % (tm_metacal,tm_metacal_per))
    print("writing:",args.outfile)
    if type(output) == list:
        output = eu.numpy_util.combine_arrlist(output)
    fitsio.write(args.outfile, output, clobber=True)

def do_profile(args):
    # don't want to see the JIT
    import cProfile
    import pstats

    import ngmix

    cProfile.runctx('main(args)',
                    globals(),locals(),
                    'profile_stats')
    
    p = pstats.Stats('profile_stats')
    p.sort_stats('time').print_stats()


if __name__=='__main__':
    args = parser.parse_args()
    if args.profile:
        do_profile(args)
    else:
        main(args)
