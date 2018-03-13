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

    def __init__(self,specs,mode,rng):
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
        Neigh = galsim.Gaussian(half_light_radius=self['Neigh']['hlr'],flux=self['Neigh']['Flux'])

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

        Cen.shift(dx=dx1, dy=dy1)
        Neigh.shift(dx=dx2, dy=dy2)

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
        Simulation.__init__(self,Sim_specs,scarlet,rng)
        #im,psf_im,coords = Simulation.__call__(self)
    
    def _get_model(self):
        im,psf_im,coords = Simulation.__call__(self)
        print(im,coords)
        bg_rms = self['Image']['Bgrms']
        if mode == 'scarlet':
            constraints = {"S": None, "m": {'use_nearest': False}, "+": None}
            sources = [scarlet.ExtendedSource(coord, im, [bg_rms]) for coord in coords]
            blend = scarlet.Blend(sources, im, bg_rms=[bg_rms])
            blend.fit(10000, e_rel=1e-1)
            model = blend.get_model()
            mod2 = blend.get_model(m=1)
            cen_mod = sources[0].get_model()
            neigh_mod = sources[1].get_model()
            steps_used = blend.it
        
            return model,mod2,cen_mod,neigh_mod,steps_used


Mod = Model()
model = Mod._get_model()
print(model)

