"""Independent quadrature and scrambled-Sobol evaluator of the stated FSO model.

This module does not import the manuscript's original simulator, its random
generator, noncentral-chi-square coupling, or tail-statistics functions.
It validates a numerical implementation of an assumed model, not an atmosphere.
"""
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import PchipInterpolator
from scipy.special import ndtri
from scipy.stats import qmc

APERTURE_RADIUS=.025
angles=np.arange(6)*np.pi/3
POSITIONS=np.vstack(([0.,0.],.075*np.column_stack((np.cos(angles),np.sin(angles)))))
PAIRS=[(),(1,2),(2,3),(3,4),(4,5),(5,6),(6,1)]
phi=np.pi/6
rotation=np.array([[np.cos(phi),-np.sin(phi)],[np.sin(phi),np.cos(phi)]])
POINTING_COV=rotation@np.diag([.035**2,.015**2])@rotation.T
separations=np.linalg.norm(POSITIONS[:,None,:]-POSITIONS[None,:,:],axis=2)
LATENT_COV=np.exp(-(separations/.055)**(5/3))+1e-12*np.eye(7)

def symmetric_root(matrix):
    val,vec=np.linalg.eigh(matrix)
    if val.min()<=0:raise ValueError('Covariance must be positive definite')
    return (vec*np.sqrt(val))@vec.T

POINTING_ROOT=symmetric_root(POINTING_COV)
LATENT_ROOT=symmetric_root(LATENT_COV)

def disk_integral(distances,width,n_radial=32,n_angle=128):
    """Direct 2D intensity integral over a circular receiver aperture.

    Substitution u=rho^2/a^2 removes the radial Jacobian. Gauss-Legendre
    integration in u and a periodic angular rule integrate the 2D Gaussian
    directly, with no Marcum Q, Bessel, or noncentral-chi-square evaluation.
    """
    distances=np.asarray(distances,float)
    shape=distances.shape;flat=distances.ravel();result=np.empty_like(flat)
    nodes,weights=leggauss(n_radial);u=(nodes+1)/2;weights=weights/2
    radii=APERTURE_RADIUS*np.sqrt(u)
    cosine=np.cos(2*np.pi*np.arange(n_angle)/n_angle)
    cross=radii[:,None]*cosine[None,:]
    radius2=radii[:,None]**2
    for start in range(0,len(flat),96):
        dd=flat[start:start+96,None,None]
        exponent=-2*(dd**2+radius2[None,:,:]-2*dd*cross[None,:,:])/width**2
        angular=np.exp(exponent).mean(axis=2)
        result[start:start+96]=2*APERTURE_RADIUS**2/width**2*(angular@weights)
    return result.reshape(shape)

class QuadratureCoupling:
    def __init__(self,theta_urad,step=.0002,n_radial=32,n_angle=128):
        self.width=np.sqrt(.010**2+(float(theta_urad)*1e-6*1000.)**2)
        self.grid=np.arange(0,.5500001,step)
        self.values=disk_integral(self.grid,self.width,n_radial,n_angle)
        self.interpolator=PchipInterpolator(self.grid,self.values,extrapolate=False)
    def __call__(self,distance):
        distance=np.asarray(distance)
        if np.any(distance<0) or np.any(distance>self.grid[-1]):
            raise ValueError('A draw exceeds the independently validated distance grid')
        return self.interpolator(distance)

def sobol_channel(n,seed,pair):
    m=int(np.log2(n))
    if 2**m!=n:raise ValueError('Sobol balance requires a power of two')
    uniforms=qmc.Sobol(d=30,scramble=True,seed=seed).random_base2(m)
    gaussian=ndtri(np.clip(uniforms[:,:9],1e-12,1-1e-12))
    pointing=np.array([.030,-.015])+gaussian[:,:2]@POINTING_ROOT.T
    logfield=gaussian[:,2:9]@LATENT_ROOT.T
    scintillation=np.exp(.38*logfield-.5*.38**2)
    retained=np.ones((n,7))
    for aperture in pair:retained[:,aperture]=.12*uniforms[:,9+aperture]
    retained*=np.where(uniforms[:,16:23]<.01,.4*uniforms[:,23:30],1.)
    return pointing,retained*scintillation

def collected(vector,channel,coupling):
    pointing,transmission_fading=channel
    center=pointing+np.asarray(vector[:2])
    difference=center[:,None,:]-POSITIONS[None,:,:]
    distance=np.sqrt(np.sum(difference*difference,axis=2))
    return np.sum(transmission_fading*coupling(distance),axis=1)

def lower_tail(values,alpha=.05):
    # Full sorting provides a separate implementation from the partition path.
    ordered=np.sort(values)
    count=int(np.ceil(alpha*len(ordered)))
    return float(np.sum(ordered[:count])/count)
