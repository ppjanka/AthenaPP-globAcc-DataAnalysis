'''
[diagnostics_vars.py]
Handles data requirements specifications for visualization globAcc project simulation results. Note: spherical_polar coordinates only!
Author: Patryk Pjanka, ppjanka@princeton.edu
'''
import numpy as np

exec(open("./diagnostics_header.py").read())

def identity (x):
    return x
def SQR (x):
    return x*x

# Account for frame rotation around the barycenter and move the velocities to a proper inertial frame.
#  - NOTE: This will give fully inertial velocities relative to the system's barycenter(!). If you instead wish to obtain frame-rotation-subtracted velocities relative to the coordinate center (M1), simply perform vel3 += omega*r.
#  - NOTE: make sure that this matches your implementation!
def inertial_frame (GM, omega, r, theta, phi, vel1, vel2, vel3):

    # transform to cylindrical coordinates
    r_cyl = r * np.sin(theta)
    # move to the barycenter system of coordinates
    r_bary = np.sqrt(SQR(r_cyl) + SQR(1.-GM) - 2.*r_cyl*(1.-GM)*np.cos(phi))
    phi_bary = np.arcsin(r_cyl*np.sin(phi)/r_bary)
    phi_bary = np.where(r_cyl*np.cos(phi) < (1.-GM), np.pi - phi_bary, phi_bary)
    # useful angle for coordinate transformations
    sin_alpha = np.sin(phi - phi_bary + 0.5*np.pi)
    cos_alpha = np.cos(phi - phi_bary + 0.5*np.pi)
    
    # frame rotation around the barycenter
    Mphi_bary = omega * r_bary

    # move back to the LAB coordinates
    dMr = Mphi_bary * cos_alpha
    dMphi = Mphi_bary * sin_alpha
    
    #update momenta
    res_vel1 = vel1 + dMr * np.sin(theta)
    res_vel2 = vel2 + dMr * np.cos(theta)
    res_vel3 = vel3 + dMphi
    
    return res_vel1, res_vel2, res_vel3

#-------------------------------------------------------------------------

class Vars ():

    def __init__ (self):
        self.quantities = []

    def quantity (self, data):
        pass
    def post_profile (self, quantity):
        # operations to be done on quantity after space-averaging and other profile-generation operations
        # e.g., division for PlasmaBeta
        return quantity

# for rho, press, vel1, etc., already in the athdf file
class Default(Vars):

    def __init__(self, var):
        super().__init__()
        self.absolute = (var[:3] == 'abs')
        if self.absolute:
            var = var[3:]
            self.process = np.abs
        else:
            self.process = identity
        self.quantities = [var,]
        self.var = var
        self.label = var

    def quantity (self,data):
        return self.process(data[self.var])

# azimuthal velocity with frame rotation and Keplerian motion subtracted
class Vel3Fluct(Vars):

    def __init__ (self, GM=GM, omega=frame_omega):
        super().__init__()
        self.quantities = ['vel3',]
        self.label = 'vel3-fluct'
        self.GM = GM
        self.omega = omega

    def quantity (self, data):
        _, _, r = np.meshgrid(data['x3v'], data['x2v'], data['x1v'], indexing='ij')
        return data['vel3'] + self.omega*r - np.sqrt(self.GM / r)

class Bfield(Vars):

    def __init__(self):
        super().__init__()
        self.quantities = ['Bcc1', 'Bcc2', 'Bcc3']
        self.label = 'Bmag'

    def quantity (self, data):
        return np.sqrt(data['Bcc1']**2 + data['Bcc2']**2 + data['Bcc3']**2)

class B2_over_B1(Vars):

    def __init__(self):
        super().__init__()
        self.quantities=['Bcc1', 'Bcc2']
        self.label='B2_over_B1'

    def quantity (self, data):
        return np.moveaxis(np.array([data['Bcc2'], data['Bcc1']]),0,-1)

    def post_profile (self, quantity):
        quantity = np.moveaxis(quantity,-1,0)
        return quantity[0] / quantity[1]

class B3_over_B1(Vars):

    def __init__(self):
        super().__init__()
        self.quantities=['Bcc1', 'Bcc3']
        self.label='B3_over_B1'

    def quantity (self, data):
        return np.moveaxis(np.array([data['Bcc3'], data['Bcc1']]),0,-1)

    def post_profile (self, quantity):
        quantity = np.moveaxis(quantity,-1,0)
        return quantity[0] / quantity[1]

class B2_over_B3(Vars):

    def __init__(self):
        super().__init__()
        self.quantities=['Bcc2', 'Bcc3']
        self.label='B2_over_B3'

    def quantity (self, data):
        return np.moveaxis(np.array([data['Bcc2'], data['Bcc3']]),0,-1)

    def post_profile (self, quantity):
        quantity = np.moveaxis(quantity,-1,0)
        return quantity[0] / quantity[1]

class PressMag(Vars):

    def __init__(self):
        super().__init__()
        self.quantities = ['Bcc1', 'Bcc2', 'Bcc3']
        self.label = 'Pmag'

    def quantity (self, data):
        return 0.25 * (data['Bcc1']**2 + data['Bcc2']**2 + data['Bcc3']**2)

class PressTot(Vars):

    def __init__(self):
        super().__init__()
        self.quantities = ['Bcc1', 'Bcc2', 'Bcc3', 'press']
        self.label = 'Ptot'

    def quantity (self, data):
        return 0.25 * (data['Bcc1']**2 + data['Bcc2']**2 + data['Bcc3']**2) \
         + data['press']

class StressReynolds(Vars):

    def __init__(self):
        super().__init__()
        self.quantities = ['rho', 'vel1', 'vel3']
        self.label = 'T_Reynolds'

    def quantity (self, data, GM=GM):
        rho, vr, vphi = map(np.array, [data['rho'], data['vel1'], data['vel3']])
        _, _, r = np.meshgrid(data['x3v'], data['x2v'], data['x1v'], indexing='ij')
        return rho * vr * (vphi + frame_omega*r - np.sqrt(GM / r))

class TurbKineticEnergy(Vars):

    def __init__(self):
        super().__init__()
        self.quantities = ['rho', 'vel1', 'vel2', 'vel3']
        self.label = 'Turb_Ekin'

    def quantity (self, data, GM=GM):
        rho, vr, vtheta, vphi = map(np.array, [data['rho'], data['vel1'], data['vel2'], data['vel3']])
        _, _, r = np.meshgrid(data['x3v'], data['x2v'], data['x1v'], indexing='ij')
        return 0.5 * rho * ( vr**2 + vtheta**2  + (vphi + frame_omega*r - np.sqrt(GM / r))**2 )

class EMF1(Vars):

    def __init__(self):
        super().__init__()
        self.quantities = ['vel2', 'vel3', 'Bcc2', 'Bcc3']
        self.label = 'EMF1'

    def quantity (self, data, GM=GM):
        _, _, r = np.meshgrid(data['x3v'], data['x2v'], data['x1v'], indexing='ij')
        return data['vel2'] * data['Bcc3'] - (data['vel3'] + frame_omega*r - np.sqrt(GM / r)) * data['Bcc2']

class EMF2(Vars):

    def __init__(self):
        super().__init__()
        self.quantities = ['vel1', 'vel3', 'Bcc1', 'Bcc3']
        self.label = 'EMF2'

    def quantity (self, data, GM=GM):
        _, _, r = np.meshgrid(data['x3v'], data['x2v'], data['x1v'], indexing='ij')
        return (data['vel3'] + frame_omega*r - np.sqrt(GM / r)) * data['Bcc1'] - data['vel1'] * data['Bcc3']

class EMF3(Vars):

    def __init__(self):
        super().__init__()
        self.quantities = ['vel1', 'vel2', 'Bcc1', 'Bcc2']
        self.label = 'EMF3'

    def quantity (self, data):
        return data['vel1'] * data['Bcc2'] - data['vel2'] * data['Bcc1']

class StressMaxwell(Vars):

    def __init__(self):
        super().__init__()
        self.quantities = ['Bcc1', 'Bcc3']
        self.label = 'T_Maxwell'

    def quantity (self, data):
        br, bphi = map(np.array, [data['Bcc1'], data['Bcc3']])
        return br * bphi

class AlphaReynolds(Vars):

    def __init__(self):
        super().__init__()
        self.quantities = ['rho', 'vel1', 'vel3', 'press']
        self.label = 'alpha_Reynolds'

    def quantity (self, data, GM=GM):
        rho, vr, vphi, press = map(np.array, [data['rho'], data['vel1'], data['vel3'], data['press']])
        _, _, r = np.meshgrid(data['x3v'], data['x2v'], data['x1v'], indexing='ij')
        r.shape = rho.shape
        return np.moveaxis(np.array([(2./3.) * rho * vr * (vphi + frame_omega*r - np.sqrt(GM / r)), press]),0,-1)

    def post_profile (self, quantity):
        quantity = np.moveaxis(quantity,-1,0)
        return quantity[0] / quantity[1]

class AlphaMaxwell(Vars):

    def __init__(self):
        super().__init__()
        self.quantities = ['Bcc1', 'Bcc3', 'press']
        self.label = 'alpha_Maxwell'

    def quantity (self, data):
        br, bphi, press = map(np.array, [data['Bcc1'], data['Bcc3'], data['press']])
        return np.moveaxis(np.array([- (2./3.) * br * bphi, press]),0,-1)

    def post_profile (self, quantity):
        quantity = np.moveaxis(quantity,-1,0)
        return quantity[0] / quantity[1]
    
class Mdot(Vars):

    def __init__(self, rho_crit=0., vel_crit=-1.0e6, ignore_floor=False, vertical=False):
        super().__init__()
        if vertical:
            self.quantities = ['rho', 'vel1', 'vel2'] # Mdot through cylindrical surface
        else:
            self.quantities = ['rho', 'vel1'] # Mdot through radial surface
        self.label = 'Mdot'
        self.rho_crit = rho_crit
        self.vel_crit = vel_crit
        self.ignore_floor = ignore_floor
        self.vertical = vertical
        if vertical:
            self.label += 'Vert'

    # WARNING: this is 3D Mdot
    if False:
        def quantity_approx (self, data):
            rho, vr, r = map(np.array, [data['rho'], data['vel1'], data['x1v']])
            dx2 = data['x2f'][1:] - data['x2f'][:-1]
            dx3 = data['x3f'][1:] - data['x3f'][:-1]
            dx3, dx2, _ = np.meshgrid(dx3, dx2, data['x1v'], indexing='ij')
            x3_mesh, x2_mesh, r = np.meshgrid(data['x3v'], data['x2v'], data['x1v'], indexing='ij')
            dS = r*dx2 * r*np.sin(x2_mesh)*dx3
            # reshape if 2D and needed
            dS.shape = rho.shape
            return - rho * vr * dS
            # sum to get total mdot over given area

    def quantity (self, data):
        rho, vr = map(np.array, [data['rho'], data['vel1']])
        if self.vertical:
            vtheta = np.array(data['vel2'])
        dx2 = - np.cos(data['x2f'][1:]) + np.cos(data['x2f'][:-1]) # delta -cos theta
        dx3 = data['x3f'][1:] - data['x3f'][:-1] # delta phi
        dx3, dx2, _ = np.meshgrid(dx3, dx2, data['x1v'], indexing='ij')
        _, theta, r = np.meshgrid(data['x3v'], data['x2v'], data['x1v'], indexing='ij')
        if self.vertical:
            dS = r*dx2*np.sin(theta) * r*dx3
        else:
            dS = r*r * dx2 * dx3
        # reshape if 2D and needed
        dS.shape = rho.shape
        if self.vertical:
            result = rho * (-vr*np.sin(theta) - vtheta*np.cos(theta)) * dS
        else:
            result = - rho * vr * dS
        if not self.ignore_floor:
            return result
        else:
            return np.where(np.logical_or(rho > self.rho_crit, vr > self.vel_crit), result, np.zeros(result.shape))
        # sum to get total mdot over given area

class Csound(Vars):

    def __init__(self):
        super().__init__()
        self.quantities = ['rho', 'press']
        self.label = 'csound'

    def csound (self, data, adiab_idx=adiab_idx):
        return np.sqrt(adiab_idx * np.array(data['press']) / np.array(data['rho']))

    def quantity (self, data, adiab_idx=adiab_idx):
        return self.csound(data, adiab_idx)

class ScaleHeight(Csound):

    def __init__(self):
        super().__init__()
        self.label = 'ScaleHeight'

    def quantity (self, data, GM=GM, adiab_idx=adiab_idx):
        _, _, r = np.meshgrid(data['x3v'], data['x2v'], data['x1v'], indexing='ij')
        r.shape = data['rho'].shape # needed for slices
        return r**1.5 * super().csound(data, adiab_idx) / np.sqrt(GM)

class PlasmaBeta(Vars):

    def __init__(self):
        super().__init__()
        self.quantities = ['press', 'Bcc1', 'Bcc2', 'Bcc3']
        self.label = 'PlasmaBeta'

    def quantity (self, data):
        return np.moveaxis(np.array([data['press'] * 4., data['Bcc1']**2 + data['Bcc2']**2 + data['Bcc3']**2 ]),0,-1)

    def post_profile (self, quantity):
        quantity = np.moveaxis(quantity,-1,0)
        return quantity[0] / quantity[1]

# TODO: not sure how to implement it at the moment
def alpha_eff (data, GM=GM, adiab_idx=adiab_idx, cutoff=np.inf):
    csound_avg = np.ma.array(csound(data, adiab_idx), mask=disk_mask(data, cutoff)).mean(axis=1)
    _, r = np.meshgrid(data['x3v'], data['x1v'], indexing='ij')
    return mdot_2D(data, cutoff=cutoff) / (3.*np.pi * sigma(data, cutoff=cutoff) * csound_avg**2 * r**1.5 / np.sqrt(GM))

class MRIstability(Csound):

    def __init__(self):
        super().__init__()
        self.quantities += ['Bcc1','Bcc2','Bcc3']
        self.label = '$H_{th}/\\lambda_{MRI,min},2\\pi/\\omega_{MRI,max}$'

    def quantity (self, data):
        _, _, r = np.meshgrid(data['x3v'], data['x2v'], data['x1v'], indexing='ij')
        r.shape = data['rho'].shape # needed for slices
        result = super().csound(data, adiab_idx)*np.sqrt( (15./(4.*np.pi)) * data['rho'] / (data['Bcc1']**2 + data['Bcc2']**2 + data['Bcc3']**2))
        return result

class KinHelicity (Vars):

    def __init__(self, GM=None, omega=None):
        super().__init__()
        self.quantities = ['vel1', 'vel2', 'vel3']
        self.label = 'KinHelicity'
        # central object mass (assumed G(M1+M2) = 1.0, and frame rotation speed)
        self.GM = GM
        self.omega = omega

    def quantity (self, data):

        from scipy.ndimage.filters import convolve

        r, theta, phi = data['x1v'], data['x2v'], data['x3v']

        if len(r) == 1 or len(theta) == 1 or len(phi) == 1:
            print("Calculation of helicity requires 3D data structure. Aborting.")
            return [None, None, None, None]

        rm, tm, pm = np.meshgrid(r, theta, phi, indexing='ij')
        dr = r[1:] - r[:-1]
        dtheta = theta[1:] - theta[:-1]
        dphi = phi[1:] - phi[:-1]
        drm, dtm, dpm = np.meshgrid(dr, dtheta, dphi, indexing='ij')
        ravg = 0.5 * (r[1:] + r[:-1])
        tavg = 0.5 * (theta[1:] + theta[:-1])
        pavg = 0.5 * (phi[1:] + phi[:-1])
        ravgm, tavgm, pavgm = np.meshgrid(ravg, tavg, pavg, indexing='ij')
        dV = ravgm**2 * np.sin(tavgm) * drm*dtm*dpm
        del drm, dtm, dpm

        vel1, vel2, vel3 = data['vel1'], data['vel2'], data['vel3']
        # ensure indices order [r, theta, phi]
        vel1 = np.transpose(vel1, (2,1,0))
        vel2 = np.transpose(vel2, (2,1,0))
        vel3 = np.transpose(vel3, (2,1,0))
        # account for frame rotation
        #vel1, vel2, vel3 = inertial_frame(self.GM, self.omega, rm, tm, pm, vel1, vel2, vel3)
        vel3 += self.omega * rm
        # average to match dimensions (n-1,n-1,n-1)
        vel1avg = convolve(vel1, np.ones([2,2,2]) * (0.5**3), mode='constant', cval=np.nan)[:-1,:-1,:-1]
        vel2avg = convolve(vel2, np.ones([2,2,2]) * (0.5**3), mode='constant', cval=np.nan)[:-1,:-1,:-1]
        vel3avg = convolve(vel3, np.ones([2,2,2]) * (0.5**3), mode='constant', cval=np.nan)[:-1,:-1,:-1]

        # -------------------------------------------------

        # nabla x vel

        # - component r
        vphi_sint = vel3 * np.sin(tm)
        dtheta_vphi_sint = (vphi_sint[:,1:,:] - vphi_sint[:,:-1,:]) / (np.tile(dtheta, [len(r), len(phi),1]).transpose(0,2,1))
        # average to match dimensions (n-1,n-1,n-1)
        dtheta_vphi_sint = 0.5 * (dtheta_vphi_sint[:,:,1:] + dtheta_vphi_sint[:,:,:-1])
        dtheta_vphi_sint = 0.5 * (dtheta_vphi_sint[1:,:,:] + dtheta_vphi_sint[:-1,:,:])
        del vphi_sint

        dphi_vtheta = (vel2[:,:,1:] - vel2[:,:,:-1]) / np.tile(dphi, [len(r), len(theta),1])
        # average to match dimensions (n-1,n-1,n-1)
        dphi_vtheta = 0.5 * (dphi_vtheta[:,1:,:] + dphi_vtheta[:,:-1,:])
        dphi_vtheta = 0.5 * (dphi_vtheta[1:,:,:] + dphi_vtheta[:-1,:,:])

        nabxvel_r = (1./(ravgm * np.sin(tavgm))) * (dtheta_vphi_sint - dphi_vtheta)
        del dtheta_vphi_sint, dphi_vtheta

        # - component theta
        dphi_vr = (vel1[:,:,1:] - vel1[:,:,:-1]) / np.tile(dphi, [len(r), len(theta),1])
        # average to match dimensions (n-1,n-1,n-1)
        dphi_vr = 0.5 * (dphi_vr[:,1:,:] + dphi_vr[:,:-1,:])
        dphi_vr = 0.5 * (dphi_vr[1:,:,:] + dphi_vr[:-1,:,:])

        rvphi = rm * vel3
        dr_rvphi = (rvphi[1:,:,:] - rvphi[:-1,:,:]) / np.tile(dr, [len(theta), len(phi),1]).transpose(2,0,1)
        # average to match dimensions (n-1,n-1,n-1)
        dr_rvphi = 0.5 * (dr_rvphi[:,1:,:] + dr_rvphi[:,:-1,:])
        dr_rvphi = 0.5 * (dr_rvphi[:,:,1:] + dr_rvphi[:,:,:-1])
        del rvphi

        nabxvel_t = (1./ravgm) * ((1./np.sin(tavgm)) * dphi_vr - dr_rvphi)
        del dphi_vr, dr_rvphi

        # - component phi
        rvtheta = rm * vel2
        dr_rvtheta = (rvtheta[1:,:,:] - rvtheta[:-1,:,:]) / np.tile(dr, [len(theta), len(phi),1]).transpose(2,0,1)
        # average to match dimensions (n-1,n-1,n-1)
        dr_rvtheta = 0.5 * (dr_rvtheta[:,1:,:] + dr_rvtheta[:,:-1,:])
        dr_rvtheta = 0.5 * (dr_rvtheta[:,:,1:] + dr_rvtheta[:,:,:-1])
        del rvtheta

        dtheta_vr = (vel1[:,1:,:] - vel1[:,:-1,:]) / (np.tile(dtheta, [len(r), len(phi),1]).transpose(0,2,1))
        # average to match dimensions (n-1,n-1,n-1)
        dtheta_vr = 0.5 * (dtheta_vr[:,:,1:] + dtheta_vr[:,:,:-1])
        dtheta_vr = 0.5 * (dtheta_vr[1:,:,:] + dtheta_vr[:-1,:,:])

        nabxvel_p = (1./ravgm) * (dr_rvtheta - dtheta_vr)
        del dr_rvtheta, dtheta_vr

        #vorticity = np.array([nabxvel_r, nabxvel_t, nabxvel_p])

        # --------------------------------
        # calculate kinetic helicity density

        helicity_dV = (vel1avg * nabxvel_r + vel2avg * nabxvel_t + vel3avg * nabxvel_p) # / dV

        # --------------------------------

        # transpose to (phi, theta, r, quantity)
        result = np.array([nabxvel_r, nabxvel_t, nabxvel_p, helicity_dV])
        result = result.transpose(3,2,1,0)
        
        # save both vorticity and helicity density
        return result
