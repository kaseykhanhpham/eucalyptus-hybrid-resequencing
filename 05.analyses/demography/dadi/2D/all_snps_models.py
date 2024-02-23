# models to test for E. globulus x E. cordata hybridization scenarios
# referenced Models_2D.py from dadi_pipeline and the dadi manual

import numpy
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum

def schange(params, ns, pts):
    # Divergence followed by exponential population size change over time
    # parameters to optimize: nu1i, nu2i, nu1f, nu2f, T
    nu1i, nu2i, nu1f, nu2f, T = params
    
    # Make grid
    xx = Numerics.default_grid(pts)
    
    # Divergence
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    # Define functions for exponential size change from nuXi to nuXf
    nu1_fun = lambda t: nu1i * (nu1f/nu1i) ** (t/T)
    nu2_fun = lambda t: nu2i * (nu2f/nu2i) ** (t/T)
    # Integrate over time T with exponential size change
    phi = Integration.two_pops(phi, xx, T, nu1 = nu1_fun, nu2 = nu2_fun, m12=0, m21=0)
    
    # Generate final expected SFS under model
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def bottle_schange(params, ns, pts):
    # Divergence followed by exponential population size change over time
    # parameters to optimize: nu1i, nu2i, nu1f, nu2f, T
    nu1i, nu2i, nu1m, nu2m, nu1f, nu2f, T1, T2 = params
    
    # Make grid
    xx = Numerics.default_grid(pts)
    
    # Divergence
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    # Integrate over initial period of stable size
    phi = Integration.two_pops(phi, xx, T1, nu1i, nu2i, m12=0, m21=0)
    
    # Define functions for exponential size change from nuXi to nuXf
    nu1_fun = lambda t: nu1m * (nu1f/nu1m) ** (t/T2)
    nu2_fun = lambda t: nu2m * (nu2f/nu2m) ** (t/T2)
    # Integrate over time T2 with exponential size change after a bottleneck
    phi = Integration.two_pops(phi, xx, T2, nu1 = nu1_fun, nu2 = nu2_fun, m12=0, m21=0)
    
    # Generate final expected SFS under model
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
 
def bottle_schange_thr_epoch(params, ns, pts):
    # Divergence followed by instantaneous size change, a period of neutral evolution,
    # and then a period of gradual exponential size change
    nu1i, nu2i, nu1m, nu2m, nu1f, nu2f, T1, T2, T3= params
    
    # Make grid
    xx = Numerics.default_grid(pts)
    
    # Divergence
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    # Integrate over initial period of stable size
    phi = Integration.two_pops(phi, xx, T1, nu1i, nu2i, m12=0, m21=0)

    # Integrate over period after instant size change
    phi = Integration.two_pops(phi, xx, T2, nu1m, nu2m, m12=0, m21=0)

    # Define functions for exponential size change from nuXi to nuXf
    nu1_fun = lambda t: nu1m * (nu1f/nu1m) ** (t/T3)
    nu2_fun = lambda t: nu2m * (nu2f/nu2m) ** (t/T3)
    # Integrate over time T3 with exponential size change after a bottleneck
    phi = Integration.two_pops(phi, xx, T3, nu1 = nu1_fun, nu2 = nu2_fun, m12=0, m21=0)
    
    # Generate final expected SFS under model
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def sec_contact_schange(params, ns, pts):
    # Divergence followed by isolation, then asymmetrical secondary contact
    # and exponential population size change over time
    nu1i, nu2i, nu1f, nu2f, m12, m21, T1, T2 = params
    
    # Make grid
    xx = Numerics.default_grid(pts)
    
    # Divergence
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    # Integrate over time spent in isolation
    phi = Integration.two_pops(phi, xx, T1, nu1i, nu2i, m12=0, m21=0)
    
    # Define functions for exponential size change from nuXi to nuXf
    nu1_fun = lambda t: nu1i * (nu1f/nu1i) ** (t/T2)
    nu2_fun = lambda t: nu2i * (nu2f/nu2i) ** (t/T2)
    # Integrate over time where species are hybridizing and population size is changing
    phi = Integration.two_pops(phi, xx, T2, nu1 = nu1_fun, nu2 = nu2_fun, m12 = m12, m21 = m21)

    # Generate final expected SFS under model
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def sec_contact_bottle_schange(params, ns, pts):
    # Divergence followed by isolation, then instantaneous size change
    # and a period of asymmetrical secondary contact and gradual size change
    nu1i, nu2i, nu1m, nu2m, nu1f, nu2f, m12, m21, T1, T2 = params

    # Make grid
    xx = Numerics.default_grid(pts)
    
    # Divergence
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    # Integrate over time spent in isolation
    phi = Integration.two_pops(phi, xx, T1, nu1i, nu2i, m12=0, m21=0)

    # Define functions for exponential size change from nuXm to nuXf
    nu1_fun = lambda t: nu1m * (nu1f/nu1m) ** (t/T2)
    nu2_fun = lambda t: nu2m * (nu2f/nu2m) ** (t/T2)
    # Integrate over time where species are hybridizing and population size is changing
    phi = Integration.two_pops(phi, xx, T2, nu1 = nu1_fun, nu2 = nu2_fun, m12 = m12, m21 = m21)

    # Generate final expected SFS under model
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def sec_contact_bottle_schange_thr_epoch(params, ns, pts):
    # Divergence followed by instantaneous size change, a period of neutral evolution,
    # and then a period of gradual exponential size change
    nu1i, nu2i, nu1m, nu2m, nu1f, nu2f, m12, m21, T1, T2, T3= params
    
    # Make grid
    xx = Numerics.default_grid(pts)
    
    # Divergence
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    # Integrate over initial period of stable size
    phi = Integration.two_pops(phi, xx, T1, nu1i, nu2i, m12=0, m21=0)

    # Integrate over period after instant size change
    phi = Integration.two_pops(phi, xx, T2, nu1m, nu2m, m12=m12, m21=m21)

    # Define functions for exponential size change from nuXi to nuXf
    nu1_fun = lambda t: nu1m * (nu1f/nu1m) ** (t/T3)
    nu2_fun = lambda t: nu2m * (nu2f/nu2m) ** (t/T3)
    # Integrate over time T3 with exponential size change after a bottleneck
    # and asymmetric migration
    phi = Integration.two_pops(phi, xx, T3, nu1 = nu1_fun, nu2 = nu2_fun, m12=0, m21=0)
    
    # Generate final expected SFS under model
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs