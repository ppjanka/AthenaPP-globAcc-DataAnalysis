'''
[diagnostics_ops.py]
Operations on data for globAcc project simulation results. Note: spherical_polar coordinates only!
Author: Patryk Pjanka, ppjanka@princeton.edu

'''
import numpy as np
from scipy.optimize import curve_fit
import pickle as pkl
import copy as cp

import matplotlib.pyplot as plt

import athena_read as ath
import diagnostics_vars as divars

exec(open("./diagnostics_header.py").read())

def del_volume (data):
    dx1 = data['x1f'][1:] - data['x1f'][:-1]
    dx2 = data['x2f'][1:] - data['x2f'][:-1]
    dx3 = data['x3f'][1:] - data['x3f'][:-1]
    x3, x2, x1 = np.meshgrid(data['x3v'], data['x2v'], data['x1v'], indexing='ij')
    dx3, dx2, dx1 = np.meshgrid(dx3, dx2, dx1, indexing='ij')
    return x1**2 * np.sin(x2) * dx1*dx2*dx3

def ignore_nan (test_array, out_arrays):
    mask = np.isnan(test_array)
    mask = np.where(test_array == np.inf, True, mask)
    mask = np.where(test_array == -np.inf, True, mask)
    return [np.ma.array(out_array, mask=mask) for out_array in out_arrays]

def cyl2sph (r_cyl, z):
    r_sph = np.sqrt(r_cyl**2+z**2)
    return [np.arccos(z/r_sph), r_sph] # theta, r_sph

# fit a step function to a radial profile
def step_function (r, rstep, val_lowr, val_highr):
    return (np.arctan((r-rstep)*200.)/np.pi+0.5) * (val_highr-val_lowr) + val_lowr
def fit_step (r, val, rstep0=0.2):
    try:
        # pre-fit levels with default rstep
        popt, pcov = curve_fit(lambda rr,val_lowr,val_highr : step_function(rr,rstep0,val_lowr,val_highr), r, val)
        # refit rstep
        popt, pcov = curve_fit(step_function, r, val, p0=[rstep0,popt[0],popt[1]])
        return lambda rr : step_function(rr, *popt)
    except Exception as e:
        return lambda rr : (rr/rr) * np.mean(val)

def get_time (filename):
    return ath.athdf(filename, quantities=[])['Time']

class Ops ():

    def __init__ (self, vars_object, \
            x1min=0.,     x2min=0.,    x3min=0., \
            x1max=np.inf, x2max=np.pi, x3max=2.*np.pi, \
            ignore_nan=True):
        self.x1min = x1min
        self.x1max = x1max
        self.x2min = x2min
        self.x2max = x2max
        self.x3min = x3min
        self.x3max = x3max
        self.vars_object = vars_object
        self.data = {}
        self.ignore_nan = ignore_nan

    # read from file and perform operations to do after read (e.g., reshape for slices)
    def read (self, filename):
        self.data = ath.athdf(filename, quantities=self.vars_object.quantities, \
        x1_min=self.x1min, x1_max=self.x1max, \
        x2_min=self.x2min, x2_max=self.x2max, \
        x3_min=self.x3min, x3_max=self.x3max)

    # return dictionary with relevant dimensions and the quantity array
    def output (self):
        return self.data

    # plot the quantity
    def plot (self, fig, subplot_spec=[1,1,1]):
        pass

    def save (self, filename):
        with open(filename, 'wb') as f:
            pkl.dump(self, f)
    def load (self, filename):
        with open(filename, 'rb') as f:
            self.__dict__ = pkl.load(f).__dict__.copy()

# --------------------------------------------------------------------

# FULL 3D (for averaging)

class Full3D (Ops):

    def __init__(self, vars_object, \
            x1min=0.,     x2min=0.,    x3min=0., \
            x1max=np.inf, x2max=np.pi, x3max=2.*np.pi):
        super().__init__(vars_object, x1min, x2min, x3min, x1max, x2max, x3max)
        self.label='Full3D'
        self.r = []; self.theta = []; self.phi = []; self.val = []
        self.time = None

    def read (self, filename):
        # turn filename to list if not list (for averaging)
        if isinstance(filename, str):
            filenames = [filename,]
        else: filenames = filename
        # read and average
        counter = 0
        for filename in filenames:
            super().read(filename)
            if len(self.phi) == 0 or len(self.r) == 0 or len(self.theta) == 0:
                self.phi = np.array(1. * self.data['x3v'])
                self.theta = np.array(1. * self.data['x2v'])
                self.r = np.array(1. * self.data['x1v'])
            if self.time == None:
                self.time = self.data['Time']
            if len(self.phi) != len(self.data['x3v']) or len(self.r) != len(self.data['x1v']) or len(self.theta) != len(self.data['x2v']):
                print(" > coordinates don't match. Ignoring frame in averaging.")
                continue
            if len(self.val) == 0:
                self.val = np.array(self.vars_object.post_profile(1. * self.vars_object.quantity(self.data)))
            else:
                self.val += np.array(self.vars_object.post_profile(self.vars_object.quantity(self.data)))
            del self.data # free memory
            counter += 1
        self.val = (np.array(self.val) / counter).transpose()
        self.data = {'phi':self.phi, 'r':self.r, 'theta':self.theta, 'val':self.val}


# --------------------------------------------------------------------

# SLICES

class EquatorialSlice (Ops):

    def __init__ (self, vars_object, \
            x1min=0.,     x2min=0.,    x3min=0., \
            x1max=np.inf, x2max=np.pi, x3max=2.*np.pi, \
            intersect=0.5*np.pi):
        super().__init__(vars_object, x1min, x2min, x3min, x1max, x2max, x3max)
        self.intersect = 0.5*np.pi
        self.x2min = self.intersect
        self.x2max = self.intersect
        self.label = 'EquatorialSlice'
        self.phi = []; self.r = []; self.val = []
        self.time = None

    def read (self, filename):
        # turn filename to list if not list (for averaging)
        if isinstance(filename, str):
            filenames = [filename,]
        else: filenames = filename
        # read and average
        counter = 0
        for filename in filenames:
            super().read(filename)
            for var in self.vars_object.quantities:
                shape = self.data[var].shape
                self.data[var].shape = [shape[0], shape[2]]
            if len(self.phi) == 0 or len(self.r) == 0:
                self.phi = 1. * self.data['x3v']
                self.r = 1. * self.data['x1v']
            if self.time == None:
                self.time = self.data['Time']
            if len(self.phi) != len(self.data['x3v']) or len(self.r) != len(self.data['x1v']) or not ((self.data['x3v'] == self.phi).all() and (self.data['x1v'] == self.r).all()):
                print(" > coordinates don't match. Ignoring frame in averaging.")
                continue
            if len(self.val) == 0:
                self.val = self.vars_object.post_profile(1. * self.vars_object.quantity(self.data))
            else:
                self.val += self.vars_object.post_profile(self.vars_object.quantity(self.data))
            del self.data # free memory
            counter += 1
        self.val = (self.val / counter).transpose()
        self.data = {'phi':self.phi, 'r':self.r, 'val':self.val}

    def plot (self, fig, ax=None, subplot_spec=111, log_scale=False, levels=10, kwargs={}, vmin=None, vmax=None, title=None):
        if ax == None:
            ax = fig.add_subplot(subplot_spec, projection='polar')
        z = self.val
        if log_scale:
            z = np.log10(self.val)
            if vmin != None: vmin = np.log10(vmin)
            if vmax != None: vmax = np.log10(vmax)
        if vmin != None and vmax != None:
            levels = np.linspace(vmin, vmax, levels)
        plot = ax.contourf(self.phi, self.r, z, levels=levels, vmin=vmin, vmax=vmax, extend='both', **kwargs)
        if title == 'auto':
            label = self.vars_object.label + ': ' + self.label
            if log_scale:
                label = 'log10 ' + label
            ax.set_title(label)
        elif title != None:
            ax.set_title(title)
        plt.setp( ax.get_xticklabels(), visible=False)
        return ax, plot

'''
    def equatorial_slice (self, data, variable):
        if len(data['x2v']) < 3:
            return data[variable][:,0,:]
        for key in ['x2v', variable]:
            data[key] = np.array(data[key])
        idx_below = max(np.where(data['x2v'] < intersect)[0])
        x1, x2 = data['x2v'][idx_below], data['x2v'][idx_below+1]
        y1, y2 = data[variable][:,idx_below,:], data[variable][:,idx_below+1,:]
        result = (y2-y1) * (intersect - x1) / (x2-x1) + y1
        return result
'''

class PoloidalSlice (Ops):

    def __init__ (self, vars_object, \
            x1min=0.,     x2min=0.,    x3min=0., \
            x1max=np.inf, x2max=np.pi, x3max=2.*np.pi, \
            intersect=1.0e-3):
        super().__init__(vars_object, x1min, x2min, x3min, x1max, x2max, x3max)
        self.intersect = intersect
        self.x3min = intersect
        self.x3max = intersect
        self.label = 'PoloidalSlice'
        self.theta = []; self.r = []; self.val = []
        self.time = None

    def poloidal_slice (data, variable):
        idx_below = min(np.where(data['x3v'] > intersect)[0]) - 1
        x1, x2 = data['x3v'][idx_below], data['x3v'][idx_below+1]
        if idx_below < 0: x1 -= 2.*np.pi
        y1, y2 = data[variable][idx_below,:,:], data[variable][idx_below+1,:,:]
        return (y2-y1) * (intersect - x1) / (x2-x1) + y1

    def read (self, filename):
        # turn filename to list if not list (for averaging)
        if isinstance(filename, str):
            filenames = [filename,]
        else: filenames = filename
        # read and average
        counter = 0
        for filename in filenames:
            super().read(filename)
            for var in self.vars_object.quantities:
                shape = self.data[var].shape
                self.data[var].shape = [shape[1], shape[2]]
            if len(self.theta) == 0 or len(self.r) == 0:
                self.theta = 1. * self.data['x2v']
                self.r = 1. * self.data['x1v']
            if self.time == None:
                self.time = self.data['Time']
            if len(self.theta) != len(self.data['x2v']) or len(self.r) != len(self.data['x1v']) or not ((self.data['x2v'] == self.theta).all() and (self.data['x1v'] == self.r).all()):
                print(" > coordinates don't match. Ignoring frame in averaging.")
                continue
            if len(self.val) == 0:
                self.val = self.vars_object.post_profile(1. * self.vars_object.quantity(self.data))
            else:
                self.val += self.vars_object.post_profile(self.vars_object.quantity(self.data))
            del self.data # free memory
            counter += 1
        self.val = (self.val / counter).transpose()
        self.data = {'theta':self.theta, 'r':self.r, 'val':self.val}

    def plot (self, fig, ax=None, subplot_spec=111, flip=False, crop=False, log_scale=False, levels=10, kwargs={}, vmin=None, vmax=None, title=None, plot_type='poloidal', grid=True, rticks=[], tick_ang=0.125*np.pi, tick_box=False, tick_linecolor='grey', rticklabels=True):
        if plot_type == 'poloidal':
            if ax == None:
                ax = fig.add_subplot(subplot_spec, projection='polar')
            z = self.val
            if log_scale:
                z = np.log10(z)
                if vmin != None: vmin = np.log10(vmin)
                if vmax != None: vmax = np.log10(vmax)
            if flip:
                ang = 0.5*np.pi+self.theta
            else:
                ang = 0.5*np.pi-self.theta
            if vmin != None and vmax != None:
                levels = np.linspace(vmin, vmax, levels)
            plot = ax.contourf(ang, self.r, z, levels=levels, vmin=vmin, vmax=vmax, extend='both', **kwargs)
            #print('PlSl: levels=%i, vmin=%f, vmax=%f' % (levels, vmin, vmax))
            if crop:
                ax.set_thetamin(min(ang))
                ax.set_thetamax(max(ang))
            del z, ang
            plt.setp( ax.get_xticklabels(), visible=False)
            if title == 'auto':
                label = self.vars_object.label + ': ' + self.label
                if log_scale:
                    label = 'log10 ' + label
                ax.set_title(label)
            elif title != None:
                ax.set_title(title)
        elif plot_type == 'box':
            if ax == None:
                ax = fig.add_subplot(subplot_spec)
            z = self.val
            if log_scale:
                z = np.log10(z)
                if vmin != None: vmin = np.log10(vmin)
                if vmax != None: vmax = np.log10(vmax)
            tick_ang = tick_ang
            if flip:
                ang = 0.5*np.pi+self.theta
                tick_ang += np.pi
            else:
                ang = 0.5*np.pi-self.theta
            if vmin != None and vmax != None:
                levels = np.linspace(vmin, vmax, levels)
            # transform to x,y coordinates
            ang_mesh, r_mesh = np.meshgrid(ang, self.r)
            ang_mesh = ang_mesh.flatten(); r_mesh = r_mesh.flatten()
            x_mesh = r_mesh * np.cos(ang_mesh)
            y_mesh = r_mesh * np.sin(ang_mesh)
            # plot
            plot = ax.tricontourf(x_mesh, y_mesh, z.flatten(), levels=levels, vmin=vmin, vmax=vmax, extend='both', **kwargs)
            # add a circle to mark the inner edge of the grid
            if grid:
                ax.add_artist(plt.Circle((0.,0.), min(self.r), color='white'))
                for rtick in rticks:
                    ax.add_artist(plt.Circle((0.,0.), rtick, color=tick_linecolor, fill=False, linewidth=0.5))
                    if rticklabels:
                        if tick_box:
                            ax.text(rtick * np.cos(tick_ang), rtick * np.sin(tick_ang), rtick, color='k', \
                                bbox=dict(boxstyle='square', ec='k', fc='white'))
                        else:
                            ax.text(rtick * np.cos(tick_ang), rtick * np.sin(tick_ang), rtick, color='k')
                x_mesh = np.array([-max(self.r), max(self.r)])
                for angle in np.linspace(0., np.pi, 5)[:-1]:
                    ax.plot(x_mesh, np.tan(angle)*x_mesh, color=tick_linecolor, linewidth=0.5)
                ax.set_xticks([])
                ax.set_yticks([])
            ax.set_aspect(1.0)
            del z, ang
            if title == 'auto':
                label = self.vars_object.label + ': ' + self.label
                if log_scale:
                    label = 'log10 ' + label
                ax.set_title(label)
            elif title != None:
                ax.set_title(title)
        return ax, plot

# --------------------------------------------------------------------

# INTEGRALS

#def disk_mask (data, cutoff=np.inf):
#    x3, x2, x1 = np.meshgrid(data['x3v'], data['x2v'], data['x1v'], indexing='ij')
#    return (np.abs(x2 - 0.5*np.pi) > cutoff)

def disk_mask (ops_object, cutoff=np.inf):
    dx2_up = 0.5*np.pi - ops_object.x2min
    if dx2_up > cutoff:
        ops_object.x2min = 0.5*np.pi - cutoff
    dx2_down = ops_object.x2max - 0.5*np.pi
    if dx2_down > cutoff:
        ops_object.x2max = 0.5*np.pi + cutoff

class RadialProfile (Ops):

    def __init__ (self, vars_object, \
            x1min=0.,     x2min=0.,    x3min=0., \
            x1max=np.inf, x2max=np.pi, x3max=2.*np.pi, \
            type='avg_vol'):
        super().__init__(vars_object, x1min, x2min, x3min, x1max, x2max, x3max)
        self.type = type
        self.label = 'RadialProfile'
        self.r = []; self.val = []
        self.time = None

    
    def radial_profile (self):
        if self.type == 'avg': # radial profiles from data by averaging
            # mask values if requested
            if self.ignore_nan:
                self.val = ignore_nan(self.val, [self.val,])[0]
            self.val = np.array(np.mean(self.val, axis=0))
            if len(self.val.shape) > 1:
                self.val = np.mean(self.val, axis=0)
        elif self.type == 'avg_vol': # volume-weighted avg
            # calculate cell volumes
            dV = del_volume(self.data)
            if len(self.val.shape) > len(dV.shape): #multidim quantities
                dV = np.repeat(np.array([dV,]), self.val.shape[-1], axis=0)
                dV = np.moveaxis(dV,0,-1)
            # mask values if requested
            if self.ignore_nan:
                self.val, dV = ignore_nan(self.val, [self.val, dV])
            # perform average
            self.val = np.sum(self.val * dV, axis=0)
            V = np.sum(dV, axis=0)
            if len(self.val.shape) > 1:
                self.val = np.sum(self.val, axis=0)
                V = np.sum(V, axis=0)
            self.val /= V
        elif self.type == 'sum': # radial profiles from data by averaging
            # mask values if requested
            if self.ignore_nan:
                self.val = ignore_nan(self.val, [self.val,])[0]
            self.val = np.array(np.sum(self.val, axis=0))
            if len(self.val.shape) > 1:
                self.val = np.sum(self.val, axis=0)
        elif self.type == 'intdS': # radial profiles from data by surface integration (for fluxes)
            r = np.array(self.data['x1v'])
            dx2 = self.data['x2f'][1:] - self.data['x2f'][:-1]
            dx3 = self.data['x3f'][1:] - self.data['x3f'][:-1]
            x3, x2, x1 = np.meshgrid(self.data['x3v'], self.data['x2v'], self.data['x1v'], indexing='ij')
            dx3, dx2, x1 = np.meshgrid(dx3, dx2, r, indexing='ij')
            self.val = np.array(np.sum(self.val*np.sin(x2)*dx3, axis=0)) # sum over phi
            self.val = np.sum(self.val*dx2, axis=0) # sum over theta
            self.val *= r**2
        elif self.type == 'surface_density':
            r = np.array(self.data['x1v'])
            L = r * (self.data['x2f'][-1] - self.data['x2f'][1])
            # calculate cell volumes
            dV = del_volume(self.data)
            if len(self.val.shape) > len(dV.shape): #multidim quantities
                dV = np.repeat(np.array([dV,]), self.val.shape[-1], axis=0)
                dV = np.moveaxis(dV,0,-1)
            # mask values if requested
            if self.ignore_nan:
                self.val, dV = ignore_nan(self.val, [self.val, dV])
            # perform average
            self.val = np.sum(self.val * dV, axis=0)
            V = np.sum(dV, axis=0)
            if len(self.val.shape) > 1:
                self.val = np.sum(self.val, axis=0)
                V = np.sum(V, axis=0)
            self.val *= L/V

    def read (self, filename):
        # turn filename to list if not list (for averaging)
        if isinstance(filename, str):
            filenames = [filename,]
        else: filenames = filename
        # read and average
        counter = 0; profile_sum = []
        for filename in filenames:
            super().read(filename)
            if len(self.r) == 0:
                self.r = 1. * self.data['x1v']
            if self.time == None:
                self.time = self.data['Time']
            if len(self.r) != len(self.data['x1v']) or not (self.r == self.data['x1v']).all():
                print(" > coordinates don't match. Ignoring frame in averaging.")
                continue
            self.val = 1. * self.vars_object.quantity(self.data)
            self.radial_profile()
            if len(profile_sum) == 0:
                profile_sum = self.vars_object.post_profile(1. * self.val)
            else:
                profile_sum += self.vars_object.post_profile(self.val)
            del self.data, self.val # free memory
            counter += 1
        self.val = profile_sum / counter
        self.data = {'r':self.r, 'val':self.val}

    def plot (self, fig, subplot_spec=111, log_scale=False, kwargs={}, vmin=None, vmax=None, title=None, fit=False, yunit=1.0):
        ax = fig.add_subplot(subplot_spec)
        z = self.val / yunit
        if vmin != None: vmin /= yunit
        if vmax != None: vmax /= yunit
        if log_scale:
            z = np.log10(z)
            if vmin != None: vmin = np.log10(vmin)
            if vmax != None: vmax = np.log10(vmax)
        plot = ax.plot(self.r, z, **kwargs)
        ax.grid()
        ax.set_ylim(vmin, vmax)
        if fit:
            z_fit = fit_step(self.r, self.val)(self.r)
            if log_scale: z_fit = np.log10(z_fit)
            plot = ax.plot(self.r, z_fit/yunit, **kwargs, linestyle=':')
        if title == 'auto':
            label = self.vars_object.label + ': ' + self.label
            if log_scale:
                label = 'log10 ' + label
            ax.set_title(label)
        elif title != None:
            ax.set_title(title)
        ax.set_xlabel('r [sim.u.]')
        return ax, plot

# radial profile averaged over zones symmetric wrt equator, i.e.,  (x2min, x2max), (pi-x2max, pi-x2min)
class RadialProfile_Straddle ():

    def __init__ (self, vars_object, \
            x1min=0.,     x2min=0.,    x3min=0., \
            x1max=np.inf, x2max=0.5*np.pi, x3max=2.*np.pi, \
            type='avg_vol'):
        self.type=type
        self.label = 'RadialProfile_Straddle'
        self.time=None
        self.r = []; self.val = []
        self.upper_ops = RadialProfile(vars_object, \
            x1min, x2min, x3min, \
            x1max, x2max, x3max, type)
        self.lower_ops = RadialProfile(vars_object, \
            x1min, np.pi-x2max, x3min, \
            x1max, np.pi-x2min, x3max, self.type)
        self.combined_ops = None

    def read (self, filename):
        self.upper_ops.read(filename)
        self.lower_ops.read(filename)
        self.time=self.upper_ops.time
        self.combined_ops = self.upper_ops
        self.combined_ops.val = 0.5 * (self.upper_ops.val + self.lower_ops.val)
        self.r = self.combined_ops.r
        self.val = self.combined_ops.val
        self.data = {'r':self.r, 'val':self.val}

    def plot (self, fig, subplot_spec=111, log_scale=False, kwargs={}, vmin=None, vmax=None, title=None, fit=False, yunit=1.0):
        # just in case, make sure the vals are the same (probably unnecessary)
        self.combined_ops.val = self.val
        return self.combined_ops.plot(fig, subplot_spec=subplot_spec, log_scale=log_scale, kwargs=kwargs, vmin=vmin, vmax=vmax, title=title, fit=fit, yunit=yunit)

class Profile_Theta (Ops):

    def __init__ (self, vars_object, \
            x1min=0.,     x2min=0.,    x3min=0., \
            x1max=np.inf, x2max=np.pi, x3max=2.*np.pi, \
            type='avg_vol'):
        super().__init__(vars_object, x1min, x2min, x3min, x1max, x2max, x3max)
        self.type = type
        self.label = 'Profile_Theta'
        self.theta = []; self.val = []
        self.time = None

    
    def profile (self):
        if self.type == 'avg': # profiles from data by averaging
            # mask values if requested
            if self.ignore_nan:
                self.val = ignore_nan(self.val, [self.val,])[0]
            self.val = np.array(np.mean(self.val, axis=0))
            if len(self.val.shape) > 1:
                self.val = np.mean(self.val, axis=1)
        elif self.type == 'avg_vol': # volume-weighted avg
            dV = del_volume(self.data)
            if len(self.val.shape) > len(dV.shape): #multidim quantities
                dV = np.repeat(np.array([dV,]), self.val.shape[-1], axis=0)
                dV = np.moveaxis(dV,0,-1)
            # mask values if requested
            if self.ignore_nan:
                self.val, dV = ignore_nan(self.val, [self.val, dV])
            self.val = np.sum(self.val * dV, axis=0) # r sum
            V = np.sum(dV, axis=0)
            if len(self.val.shape) > 1:
                self.val = np.sum(self.val, axis=1) # phi sum
                V = np.sum(V, axis=1)
            self.val /= V
        elif self.type == 'sum': # profiles from data by sum
            # mask values if requested
            if self.ignore_nan:
                self.val = ignore_nan(self.val, [self.val,])[0]
            self.val = np.array(np.sum(self.val, axis=0))
            if len(self.val.shape) > 1:
                self.val = np.sum(self.val, axis=1)
        elif self.type == 'intdS': # profiles from data by surface integration (for fluxes), returns dF/rdTheta
            dx1 = self.data['x1f'][1:] - self.data['x1f'][:-1]
            dx2 = self.data['x2f'][1:] - self.data['x2f'][:-1]
            dx3 = self.data['x3f'][1:] - self.data['x3f'][:-1]
            x3, x2, x1 = np.meshgrid(self.data['x3v'], self.data['x2v'], self.data['x1v'], indexing='ij')
            dx3, dx2, dx1 = np.meshgrid(dx3, dx2, dx1, indexing='ij')
            self.val = np.sum(self.val*x1*np.sin(x2)*dx3, axis=0) # sum over phi
            self.val = np.sum(self.val*dx1, axis=1) # sum over r

    def read (self, filename):
        # turn filename to list if not list (for averaging)
        if isinstance(filename, str):
            filenames = [filename,]
        else: filenames = filename
        # read and average
        counter = 0; profile_sum = []
        for filename in filenames:
            super().read(filename)
            if len(self.theta) == 0:
                self.theta = 1. * self.data['x2v']
            if self.time == None:
                self.time = self.data['Time']
            if len(self.theta) != len(self.data['x2v']) or not (self.theta == self.data['x2v']).all():
                print(" > coordinates don't match. Ignoring frame in averaging.")
                continue
            self.val = 1. * self.vars_object.quantity(self.data)
            self.profile()
            if len(profile_sum) == 0:
                profile_sum = self.vars_object.post_profile(1. * self.val)
            else:
                profile_sum += self.vars_object.post_profile(self.val)
            del self.data, self.val # free memory
            counter += 1
        self.val = profile_sum / counter
        self.data = {'theta':self.theta, 'val':self.val}

    def plot (self, fig, subplot_spec=111, log_scale=False, vertical=False, kwargs={}, vmin=None, vmax=None, title=None, theta_H=0.2):
        ax = fig.add_subplot(subplot_spec)
        x = self.theta
        y = self.val
        if log_scale:
            y = np.log10(self.val)
            if vmin != None: vmin = np.log10(vmin)
            if vmax != None: vmax = np.log10(vmax)
        if vertical:
            x,y = y,x
        plot = ax.plot(x, y, **kwargs)
        _ = ax.grid()
        if vertical:
            _ = ax.set_ylim(0., np.pi)
            _ = ax.axhline(0.5*np.pi, color='k', linewidth=1.5, linestyle='-')
            _ = ax.axhline(0.5*np.pi-5.*theta_H, color='r', linewidth=1.5, linestyle=':')
            _ = ax.axhline(0.5*np.pi+5.*theta_H, color='r', linewidth=1.5, linestyle=':')
            _ = ax.set_xlim(vmin, vmax)
            ax.set_ylabel('Theta')
        if title == 'auto':
            label = self.vars_object.label + ': ' + self.label
            if log_scale:
                label = 'log10 ' + label
            ax.set_title(label)
        elif title != None:
            ax.set_title(title)
        return ax, plot

class Profile_PhiR (Ops):

    def __init__ (self, vars_object, \
            x1min=0.,     x2min=0.,    x3min=0., \
            x1max=np.inf, x2max=np.pi, x3max=2.*np.pi, \
            type='avg_vol'):
        super().__init__(vars_object, x1min, x2min, x3min, x1max, x2max, x3max)
        self.type = type
        self.label = 'Profile_PhiR'
        self.phi = []; self.r = []; self.val = []
        self.time = None
    
    def profile (self):
        if self.type == 'avg': # equatorial profile averaged over theta
            # mask values if requested
            if self.ignore_nan:
                self.val = ignore_nan(self.val, [self.val,])[0]
            self.val = np.array(np.mean(self.val, axis=1)) # avg over theta
        elif self.type == 'avg_vol': # volume-weighted avg
            dV = del_volume(self.data)
            if len(self.val.shape) > len(dV.shape): #multidim quantities
                dV = np.repeat(np.array([dV,]), self.val.shape[-1], axis=0)
                dV = np.moveaxis(dV,0,-1)
            # mask values if requested
            if self.ignore_nan:
                self.val, dV = ignore_nan(self.val, [self.val, dV])
            self.val = np.sum(self.val * dV, axis=1) / np.sum(dV, axis=1)
        elif self.type == 'sum': # equatorial profile summed over theta
            # mask values if requested
            if self.ignore_nan:
                self.val = ignore_nan(self.val, [self.val,])[0]
            self.val = np.array(np.sum(self.val, axis=1)) # sum over theta
        elif self.type == 'intdTheta': # equatorial profile integrated over r*dtheta, e.g., surface density
            dx2 = self.data['x2f'][1:] - self.data['x2f'][:-1]
            x3, dx2, x1 = np.meshgrid(self.data['x3v'], dx2, self.data['x1v'], indexing='ij')
            self.val = np.sum(self.val*x1*dx2, axis=1)
        elif self.type == 'intdV': # equatorial profile integrated over volume, e.g., total mass per pixel
            dx1 = self.data['x1f'][1:] - self.data['x1f'][:-1]
            dx2 = self.data['x2f'][1:] - self.data['x2f'][:-1]
            dx3 = self.data['x3f'][1:] - self.data['x3f'][:-1]
            x3, x2, x1 = np.meshgrid(self.data['x3v'], self.data['x2v'], self.data['x1v'], indexing='ij')
            dx3, dx2, dx1 = np.meshgrid(dx3, dx2, dx1, indexing='ij')
            dV = x1**2 * np.sin(x2) * dx1*dx2*dx3
            self.val = np.sum(self.val*dV, axis=1)

    def read (self, filename):
        # turn filename to list if not list (for averaging)
        if isinstance(filename, str):
            filenames = [filename,]
        else: filenames = filename
        # read and average
        counter = 0; profile_sum = []
        for filename in filenames:
            super().read(filename)
            if len(self.phi) == 0 or len(self.r) == 0:
                self.phi = 1. * self.data['x3v']
                self.r = 1. * self.data['x1v']
            if self.time == None:
                self.time = self.data['Time']
            if len(self.phi) != len(self.data['x3v']) or len(self.r) != len(self.data['x1v']) or not ((self.phi == self.data['x3v']).all() and (self.r == self.data['x1v']).all()):
                print(" > coordinates don't match. Ignoring frame in averaging.")
                continue
            self.val = (1. * self.vars_object.quantity(self.data))
            self.profile()
            if len(profile_sum) == 0:
                profile_sum = self.vars_object.post_profile(1. * self.val)
            else:
                profile_sum += self.vars_object.post_profile(self.val)
            del self.data, self.val # free memory
            counter += 1
        self.val = (profile_sum / counter).transpose()
        self.data = {'phi':self.phi, 'r':self.r, 'val':self.val}

    def plot (self, fig, subplot_spec=111, log_scale=False, levels=10, kwargs={}, vmin=None, vmax=None, title=None):
        ax = fig.add_subplot(subplot_spec, projection='polar')
        z = self.val
        if log_scale:
            z = np.log10(z)
            if vmin != None: vmin = np.log10(vmin)
            if vmax != None: vmax = np.log10(vmax)
        if vmin != None and vmax != None:
            levels = np.linspace(vmin, vmax, levels)
        plot = ax.contourf(self.phi, self.r, z, levels=levels, vmin=vmin, vmax=vmax, extend='both', **kwargs)
        if title == 'auto':
            label = self.vars_object.label + ': ' + self.label
            if log_scale:
                label = 'log10 ' + label
            ax.set_title(label)
        elif title != None:
            ax.set_title(title)
        plt.setp( ax.get_xticklabels(), visible=False)
        return ax, plot

class Profile_ThetaR (Ops):

    def __init__ (self, vars_object, \
            x1min=0.,     x2min=0.,    x3min=0., \
            x1max=np.inf, x2max=np.pi, x3max=2.*np.pi, \
            type='avg'):
        super().__init__(vars_object, x1min, x2min, x3min, x1max, x2max, x3max)
        self.type = type
        self.label = 'Profile_ThetaR'
        self.theta = []; self.r = []; self.val = []
        self.time = None

    def profile (self):
        if self.type == 'avg' or self.type == 'avg_vol': # poloidal profile averaged over phi
            # mask values if requested
            if self.ignore_nan:
                self.val = ignore_nan(self.val, [self.val,])[0]
            self.val = np.array(np.mean(self.val, axis=0)) # avg over phi
        elif self.type == 'sum': # poloidal profile summed over phi
            # mask values if requested
            if self.ignore_nan:
                self.val = ignore_nan(self.val, [self.val,])[0]
            self.val = np.array(np.sum(self.val, axis=0)) # sum over phi
        elif self.type == 'intdPhi': # poloidal profile integrated over r sinTheta dphi, e.g., surface density
            dx3 = self.data['x3f'][1:] - self.data['x3f'][:-1]
            dx3, x2, x1 = np.meshgrid(dx3, self.data['x2v'], self.data['x1v'], indexing='ij')
            self.val = np.array(np.sum(self.val*x1*np.sin(x2)*dx3, axis=0)) # integrate over phi
        elif self.type == 'intdV': # poloidal profile integrated over volume, e.g., total mass per pixel
            dx1 = self.data['x1f'][1:] - self.data['x1f'][:-1]
            dx2 = self.data['x2f'][1:] - self.data['x2f'][:-1]
            dx3 = self.data['x3f'][1:] - self.data['x3f'][:-1]
            x3, x2, x1 = np.meshgrid(self.data['x3v'], self.data['x2v'], self.data['x1v'], indexing='ij')
            dx3, dx2, dx1 = np.meshgrid(dx3, dx2, dx1, indexing='ij')
            dV = x1**2 * np.sin(x2) * dx1*dx2*dx3
            self.val = np.sum(self.val*dV, axis=0)

    def read (self, filename):
        # turn filename to list if not list (for averaging)
        if isinstance(filename, str):
            filenames = [filename,]
        else: filenames = filename
        # read and average
        counter = 0; profile_sum = []
        for filename in filenames:
            super().read(filename)
            if len(self.theta) == 0 or len(self.r) == 0:
                self.theta = 1. * self.data['x2v']
                self.r = 1. * self.data['x1v']
            if self.time == None:
                self.time = self.data['Time']
            if len(self.theta) != len(self.data['x2v']) or len(self.r) != len(self.data['x1v']) or not ((self.theta == self.data['x2v']).all() and (self.r == self.data['x1v']).all()):
                print(" > coordinates don't match. Ignoring frame in averaging.")
                continue
            self.val = (1. * self.vars_object.quantity(self.data))
            self.profile()
            if len(profile_sum) == 0:
                profile_sum = self.vars_object.post_profile(1. * self.val)
            else:
                profile_sum += self.vars_object.post_profile(self.val)
            del self.data, self.val # free memory
            counter += 1
        self.val = (profile_sum / counter).transpose()
        self.data = {'theta':self.theta, 'r':self.r, 'val':self.val}

    def plot (self, fig, ax=None, subplot_spec=111, flip=False, crop=False, log_scale=False, levels=10, kwargs={}, vmin=None, vmax=None, title=None, plot_type='poloidal', grid=True, rticks=[], tick_ang=0.125*np.pi, rticklabels=True):
        if plot_type == 'poloidal':
            if ax == None:
                ax = fig.add_subplot(subplot_spec, projection='polar')
            z = self.val
            if log_scale:
                z = np.log10(z)
                if vmin != None: vmin = np.log10(vmin)
                if vmax != None: vmax = np.log10(vmax)
            if flip:
                ang = 0.5*np.pi+self.theta
            else:
                ang = 0.5*np.pi-self.theta
            if vmin != None and vmax != None:
                levels = np.linspace(vmin, vmax, levels)
            plot = ax.contourf(ang, self.r, z, levels=levels, vmin=vmin, vmax=vmax, extend='both', **kwargs)
            if crop:
                ax.set_thetamin(min(ang)*180./np.pi)
                ax.set_thetamax(max(ang)*180./np.pi)
            ax.set_xticks(np.linspace(-0.5*np.pi, 1.5*np.pi, 9))
            plt.setp( ax.get_xticklabels(), visible=False)
            del z, ang
            if title == 'auto':
                label = self.vars_object.label + ': ' + self.label
                if log_scale:
                    label = 'log10 ' + label
                ax.set_title(label)
            elif title != None:
                ax.set_title(title)
        elif plot_type == 'box':
            if ax == None:
                ax = fig.add_subplot(subplot_spec)
            z = self.val
            if log_scale:
                z = np.log10(z)
                if vmin != None: vmin = np.log10(vmin)
                if vmax != None: vmax = np.log10(vmax)
            tick_ang = tick_ang
            if flip:
                ang = 0.5*np.pi+self.theta
                tick_ang += np.pi
            else:
                ang = 0.5*np.pi-self.theta
            if vmin != None and vmax != None:
                levels = np.linspace(vmin, vmax, levels)
            # transform to x,y coordinates
            ang_mesh, r_mesh = np.meshgrid(ang, self.r)
            ang_mesh = ang_mesh.flatten(); r_mesh = r_mesh.flatten()
            x_mesh = r_mesh * np.cos(ang_mesh)
            y_mesh = r_mesh * np.sin(ang_mesh)
            # plot
            plot = ax.tricontourf(x_mesh, y_mesh, z.flatten(), levels=levels, vmin=vmin, vmax=vmax, extend='both', **kwargs)
            # add a circle to mark the inner edge of the grid
            if grid:
                ax.add_artist(plt.Circle((0.,0.), min(self.r), color='white'))
                for rtick in rticks:
                    ax.add_artist(plt.Circle((0.,0.), rtick, color='grey', fill=False, linewidth=0.5))
                    if rticklabels:
                        ax.text(rtick * np.cos(tick_ang), rtick * np.sin(tick_ang), rtick, color='k')
                ax.add_artist(plt.Circle((0.,0.), max(self.r), color='k', fill=False, linewidth=1))
                for angle in np.linspace(0., np.pi, 5)[:-1]:
                    x_mesh = np.array([-max(self.r)*np.cos(angle), max(self.r)*np.cos(angle)])
                    ax.plot(x_mesh, np.tan(angle)*x_mesh, color='grey', linewidth=0.5)
                ax.set_xticks([])
                ax.set_yticks([])
            ax.set_aspect(1.0)
            del z, ang
            if title == 'auto':
                label = self.vars_object.label + ': ' + self.label
                if log_scale:
                    label = 'log10 ' + label
                ax.set_title(label)
            elif title != None:
                ax.set_title(title)
        return ax, plot

# --------------------------------------------------------------------

# TIME SERIES

# a value integrated over a single sphere, e.g. Mdot through a given radius
class TSeries_SphSlice (Ops):

    def __init__ (self, vars_object, r=0.3, \
            x2min=0.,    x3min=0., \
            x2max=np.pi, x3max=2.*np.pi, \
            type='sum'):
        super().__init__(vars_object, r, x2min, x3min, r, x2max, x3max)
        self.type = type
        self.label = 'TSeries_SphSlice'
        self.times = []; self.values = []
    
    def integrate (self):
        if self.type == 'avg': # radial profiles from data by averaging
            # mask values if requested
            if self.ignore_nan:
                self.val = ignore_nan(self.val, [self.val,])[0]
            while len(self.val.shape) > 0:
                self.val = np.mean(self.val, axis=0)
        elif self.type == 'sum': # radial profiles from data by averaging
            # mask values if requested
            if self.ignore_nan:
                self.val = ignore_nan(self.val, [self.val,])[0]
            while len(self.val.shape) > 0:
                self.val = np.sum(self.val, axis=0)
        elif self.type == 'int_dV': # volume-weighted avg
            # calculate cell volumes
            dV = del_volume(self.data)
            if len(self.val.shape) > len(dV.shape): #multidim quantities
                dV = np.repeat(np.array([dV,]), self.val.shape[-1], axis=0)
                dV = np.moveaxis(dV,0,-1)
            # mask values if requested
            if self.ignore_nan:
                self.val, dV = ignore_nan(self.val, [self.val, dV])
            # perform average
            self.val *= dV
            while len(self.val.shape) > 0:
                self.val = np.sum(self.val, axis=0)
        elif self.type == 'avg_vol': # volume-weighted avg
            # calculate cell volumes
            dV = del_volume(self.data)
            if len(self.val.shape) > len(dV.shape): #multidim quantities
                dV = np.repeat(np.array([dV,]), self.val.shape[-1], axis=0)
                dV = np.moveaxis(dV,0,-1)
            # mask values if requested
            if self.ignore_nan:
                self.val, dV = ignore_nan(self.val, [self.val, dV])
            # perform average
            self.val *= dV
            while len(self.val.shape) > 0:
                self.val = np.sum(self.val, axis=0)
                V = np.sum(V, axis=0)
            self.val /= V

    def read (self, filename, old_data=None):
        # turn filename to list if not list (for averaging)
        if isinstance(filename, str):
            filenames = [filename,]
        else: filenames = filename
        # read previous results, if available
        time_done = 0.
        if old_data != None:
            self.times = list(old_data['t'])
            self.values = list(old_data['val'])
            time_done = max(self.times)
        # read and process
        from tqdm import tqdm
        for filename in tqdm(filenames):
            try:
                timeonly = ath.athdf(filename, quantities=[], \
                    x1_min=0.3, x1_max=0.3, \
                    x2_min=1.0, x2_max=1.0, \
                    x3_min=1.0, x3_max=1.0)
                if timeonly['Time'] > time_done:
                    super().read(filename)
                    self.val = 1. * self.vars_object.quantity(self.data)
                    self.integrate()
                    self.times.append(1.*self.data['Time'])
                    self.values.append(1.*self.val)
                    del self.data, self.val # free memory
                else:
                    print(' -- file %s already processed, skipping.' % filename)
                del timeonly
            except Exception as e:
                print('Could not read file: %s' % filename)
                print(e)
        self.times = np.array(self.times)
        self.values = np.array(self.values)
        self.data = {'t':self.times, 'val':self.values}

    def plot (self, fig, ax=None, subplot_spec=111, log_scale=False, kwargs={}, vmin=None, vmax=None, title=None, boxcar_smooth=0, yunit=1.0):
        if ax == None:
            ax = fig.add_subplot(subplot_spec)
        z = self.values
        t = self.times
        if boxcar_smooth > 0:
            box = np.ones(boxcar_smooth)/boxcar_smooth
            z = np.convolve(z, box, mode='valid')
            t = np.convolve(t, box, mode='valid') # kinda stupid, but ensures correctness
        if log_scale:
            z = np.log10(z)
            if vmin != None: vmin = np.log10(vmin)
            if vmax != None: vmax = np.log10(vmax)
        plot = ax.plot(t, z/yunit, **kwargs)
        ax.grid()
        if vmin != None or vmax != None:
            ax.set_ylim(vmin, vmax)
        if title == 'auto':
            label = self.vars_object.label + ': ' + self.label
            if log_scale:
                label = 'log10 ' + label
            if boxcar_smooth > 0:
                label = label + ' (' + str(boxcar_smooth) + '-frame avg)'
            ax.set_title(label)
        elif title != None:
            ax.set_title(title)
        ax.set_xlabel('Time [sim.u.]')
        return ax, plot



# a value integrated over the whole domain
class TSeries_Total (Ops):

    def __init__ (self, vars_object, \
            x1min=0., x2min=0.,    x3min=0., \
            x1max=np.inf, x2max=np.pi, x3max=2.*np.pi, \
            type='int_dV'):
        super().__init__(vars_object, x1min, x2min, x3min, x1max, x2max, x3max)
        self.type = type
        self.label = 'TSeries_Total'
        self.times = []; self.values = []
    
    def integrate (self):
        if self.type == 'avg': # radial profiles from data by averaging
            # mask values if requested
            if self.ignore_nan:
                self.val = ignore_nan(self.val, [self.val,])[0]
            while len(self.val.shape) > 0:
                self.val = np.mean(self.val, axis=0)
        elif self.type == 'sum': # radial profiles from data by averaging
            # mask values if requested
            if self.ignore_nan:
                self.val = ignore_nan(self.val, [self.val,])[0]
            while len(self.val.shape) > 0:
                self.val = np.sum(self.val, axis=0)
        elif self.type == 'int_dV': # volume-weighted avg
            # calculate cell volumes
            dV = del_volume(self.data)
            if len(self.val.shape) > len(dV.shape): #multidim quantities
                dV = np.repeat(np.array([dV,]), self.val.shape[-1], axis=0)
                dV = np.moveaxis(dV,0,-1)
            # mask values if requested
            if self.ignore_nan:
                self.val, dV = ignore_nan(self.val, [self.val, dV])
            # perform average
            self.val *= dV
            while len(self.val.shape) > 0:
                self.val = np.sum(self.val, axis=0)
        elif self.type == 'avg_vol': # volume-weighted avg
            # calculate cell volumes
            dV = del_volume(self.data)
            if len(self.val.shape) > len(dV.shape): #multidim quantities
                dV = np.repeat(np.array([dV,]), self.val.shape[-1], axis=0)
                dV = np.moveaxis(dV,0,-1)
            # mask values if requested
            if self.ignore_nan:
                self.val, dV = ignore_nan(self.val, [self.val, dV])
            # perform average
            self.val *= dV
            while len(self.val.shape) > 0:
                self.val = np.sum(self.val, axis=0)
                V = np.sum(V, axis=0)
            self.val /= V

    def read (self, filename, old_data=None):
        # turn filename to list if not list (for averaging)
        if isinstance(filename, str):
            filenames = [filename,]
        else: filenames = filename
        # read previous results, if available
        time_done = 0.
        if old_data != None:
            self.times = list(old_data['t'])
            self.values = list(old_data['val'])
            time_done = max(self.times)
        # read and process
        from tqdm import tqdm
        for filename in tqdm(filenames):
            try:
                timeonly = ath.athdf(filename, quantities=[], \
                    x1_min=0.3, x1_max=0.3, \
                    x2_min=1.0, x2_max=1.0, \
                    x3_min=1.0, x3_max=1.0)
                if timeonly['Time'] > time_done:
                    super().read(filename)
                    self.val = 1. * self.vars_object.quantity(self.data)
                    self.integrate()
                    self.times.append(1.*self.data['Time'])
                    self.values.append(1.*self.val)
                    del self.data, self.val # free memory
                else:
                    print(' -- file %s already processed, skipping.' % filename)
                del timeonly
            except Exception as e:
                print('Could not read file: %s' % filename)
                print(e)
        self.times = np.array(self.times)
        self.values = np.array(self.values)
        self.data = {'t':self.times, 'val':self.values}

    def plot (self, fig, ax=None, subplot_spec=111, log_scale=False, kwargs={}, vmin=None, vmax=None, title=None):
        if ax == None:
            ax = fig.add_subplot(subplot_spec)
        z = self.values
        if log_scale:
            z = np.log10(z)
            if vmin != None: vmin = np.log10(vmin)
            if vmax != None: vmax = np.log10(vmax)
        plot = ax.plot(self.times, z, **kwargs)
        ax.grid()
        if vmin != None or vmax != vmax:
            ax.set_ylim(vmin, vmax)
        if title == 'auto':
            label = self.vars_object.label + ': ' + self.label
            if log_scale:
                label = 'log10 ' + label
            ax.set_title(label)
        elif title != None:
            ax.set_title(title)
        ax.set_xlabel('Time [sim.u.]')
        return ax, plot

# ------------------------------------------------------------------------
# CYLINDRICAL COORDINATES

# interpolate val in slices of phi, re-grid into radii, average over phi
def phi_slice_values (arg):
    self, phi_idx = arg
    from scipy.interpolate import RegularGridInterpolator as RGI
    result = []
    val = self.vars_object.quantity(self.data)
    val_slice = val[phi_idx,:,:]
    val_fun = RGI(points=(self.data['x2v'], self.data['x1v']), \
        values=val_slice, method='linear', bounds_error=False) # val(theta,r)
    for r in self.radii:
        points = np.array(cyl2sph(r, self.z_over_r * r)).transpose()
        result.append(np.concatenate(np.array([val_fun(point) for point in points])))
    return np.array(result)

# extracts cylindrical-vertical profiles for the given list of radii
# NOTE: assumes spherical_polar input coords
class Profile_Cyl_Vertical (Ops):

    def __init__ (self, vars_object, radii, nz=128, \
            x1min=0., x2min=0.1,    x3min=0., \
            x1max=np.inf, x2max=np.pi-0.1, x3max=2.*np.pi, \
            type='avg', nproc=1):
        x1min = max(x1min, 0.9*min(radii))
        x1max = min(x1max, 1.1*max(radii) / np.cos(max([np.abs(x-0.5*np.pi) for x in [x2min, x2max]])))
        super().__init__(vars_object, \
            x1min, x2min, x3min, \
            x1max, x2max, x3max)
        self.type = type
        self.label = 'cyl_vertical_profile'
        self.radii = radii
        self.z_over_r = np.linspace(1./np.tan(x2max), 1./np.tan(x2min), nz)
        self.values = {} # phi-ageraged, single frame
        self.val = {} # frame-averaged
        self.nproc = nproc

    def integrate (self):
        from multiprocessing import Pool

        args = [(self, i) for i in range(len(self.data['x3v']))]
        with Pool(self.nproc) as pool:
            phi_slice_values_avg = pool.map(phi_slice_values, args)
        phi_slice_values_avg = np.sum(phi_slice_values_avg, axis=0)

        self.values = {}
        for ridx in range(len(self.radii)):
            self.values[self.radii[ridx]] = phi_slice_values_avg[ridx]

        if self.type == 'avg':
            for r in self.radii:
                self.values[r] /= len(self.data['x3v'])
        elif self.type == 'sum':
            pass # just leave as a sum

    def read (self, filename):
        # turn filename to list if not list (for averaging)
        if isinstance(filename, str):
            filenames = [filename,]
        else: filenames = filename
        # read and average
        counter = 0; profile_sum = []
        for filename in filenames:
            super().read(filename)
            self.integrate()
            del self.data
            for r in self.radii:
                if r not in self.val.keys():
                    self.val[r] = cp.deepcopy(self.vars_object.post_profile(self.values[r]))
                else:
                    self.val[r] += self.vars_object.post_profile(self.values[r])
            counter += 1
            del self.values
        for r in self.radii:
            self.val[r] /= counter
        self.data = {'radii':self.radii, 'z/r':self.z_over_r, 'val':self.val}

    def read_from_3D (self, ops_obj_3D):
        from scipy.interpolate import RegularGridInterpolator as RGI
        from multiprocessing import Pool
        # read in and average if needed
        self.data['x1v'] = ops_obj_3D.r
        self.data['x2v'] = ops_obj_3D.theta
        if self.x3min != self.x3max: # avg
            self.data['x3v'] = ops_obj_3D.phi
            if self.type == 'avg':
                val_slice = np.mean(ops_obj_3D.val, axis=0)
            elif self.type == 'sum':
                val_slice = np.sum(ops_obj_3D.val, axis=0)
        else: # slice
            phi_idx = max(np.where(ops_obj_3D.phi <= self.x3min)[0])
            self.data['x3v'] = ops_obj_3D.phi[phi_idx]
            val_slice = ops_obj_3D.val[phi_idx,:,:]
        # interpolate
        val_fun = RGI(points=(self.data['x2v'], self.data['x1v']), \
        values=val_slice, method='linear', bounds_error=False) # val(theta,r)
        for r in self.radii:
            points = np.array(cyl2sph(r, self.z_over_r * r)).transpose()
            with Pool(self.nproc) as pool:
                self.val[r] = np.array(pool.map(val_fun, points))
        self.data = {'radii':self.radii, 'z/r':self.z_over_r, 'val':self.val}

    def plot (self, fig, cmap=None, color_fun=None, ax=None, subplot_spec=111, log_scale=False, kwargs={}, vmin=None, vmax=None, title=None, legend=False, yunit=1.0):
        if ax == None:
            ax = fig.add_subplot(subplot_spec)
        if log_scale:
            ax.set_yscale('log')
        rmin, rmax = min(self.radii), max(self.radii)
        if len(self.radii) > 1:
            for r in self.radii:
                if color_fun == None:
                    plot = ax.plot(self.z_over_r, self.val[r]/yunit, color=cmap((r-rmin)/(rmax-rmin)), label=('$r=%.2f$' % r), **kwargs)
                else:
                    plot = ax.plot(self.z_over_r, self.val[r]/yunit, color=color_fun(r), label=('$r=%.2f$' % r), **kwargs)
        else:
            plot = ax.plot(self.z_over_r, self.val[self.radii[0]]/yunit, label=('$r=%.1f$' % self.radii[0]), **kwargs)
        ax.grid()
        if legend:
            ax.legend()
        if vmin != None or vmax != None:
            ax.set_ylim(vmin, vmax)
        if title == 'auto':
            label = self.vars_object.label + ': ' + self.label
            if log_scale:
                label = 'log10 ' + label
            ax.set_title(label)
        elif title != None:
            ax.set_title(title)
        ax.set_xlabel('$z/r$')
        return ax, plot

# extracts a butterfly diagram (BCCn at given radii averaged azimuthally vs time)
# NOTE: assumes spherical_polar input coords
class ButterflyDiagram (Ops):

    def __init__ (self, vars_object, radii, nz=128, \
            x1min=0., x2min=0.1,    x3min=0., \
            x1max=np.inf, x2max=np.pi-0.1, x3max=2.*np.pi, \
            type='avg'):
        x1min = max(x1min, 0.9*min(radii))
        x1max = min(x1max, 1.1*max(radii) / np.cos(max([np.abs(x-0.5*np.pi) for x in [x2min, x2max]])))
        super().__init__(vars_object, \
            x1min, x2min, x3min, \
            x1max, x2max, x3max)
        self.type = type
        self.label = 'ButterflyDiagram'
        self.radii = radii
        self.z_over_r = np.linspace(1./np.tan(x2max), 1./np.tan(x2min), nz)
        self.values = {} # frames
        self.val = {} # output dict
        self.times = [] # holds all the times read so far (for restarts)

    def integrate (self):
        from scipy.interpolate import RegularGridInterpolator as RGI
        
        val = self.vars_object.quantity(self.data)

        # average over phi unless slice
        if len(val.shape) > 2:
            if self.type == 'avg':
                val = np.mean(val, axis=0)
            elif self.type == 'sum':
                val = np.sum(val, axis=0)

        val_fun = RGI(points=(self.data['x2v'], self.data['x1v']), \
            values=val, method='linear', bounds_error=False) # val(theta,r)
        for r in self.radii:
            points = np.array(cyl2sph(r, self.z_over_r * r)).transpose()
            self.values[r] = np.concatenate(np.array([val_fun(point) for point in points]))

    def read (self, filename, verbose=False, temp_save_freq=5, temp_save_path=None):
        # turn filename to list if not list (for averaging)
        if isinstance(filename, str):
            filenames = [filename,]
        else: filenames = filename
        # prepare for saves if requested
        if temp_save_path != None:
            import pickle as pkl
        # read the data
        for idx in range(len(filenames)):
            filename = filenames[idx]
            if verbose: print(' - reading %s.. ' % filename.split('/')[-1], flush=True, end='')
            self.values = {}
            super().read(filename)
            time = self.data['Time']
            if time in self.times:
                if verbose: print('already read. Skipping.', flush=True)
                continue
            self.integrate()
            del self.data
            for r in self.radii:
                if r not in self.val.keys():
                    self.val[r] = {}
                self.val[r][time] = cp.deepcopy(self.vars_object.post_profile(self.values[r]))
            self.times.append(time)
            if verbose: print('done.', flush=True)
            del self.values
            # save for restart if requested
            if temp_save_path != None and idx % temp_save_freq == 0:
                with open(temp_save_path, 'wb') as f:
                    pkl.dump(self, f)
        self.data = {'radii':self.radii, 'z/r':self.z_over_r, 'val':self.val}
        if verbose: print('All files read.', flush=True)

    # TODO: update for butterfly diagrams
    def plot (self, fig, cmap=None, color_fun=None, ax=None, subplot_spec=111, log_scale=False, kwargs={}, vmin=None, vmax=None, title=None, legend=False, yunit=1.0):
        if ax == None:
            ax = fig.add_subplot(subplot_spec)
        if log_scale:
            ax.set_yscale('log')
        rmin, rmax = min(self.radii), max(self.radii)
        if len(self.radii) > 1:
            for r in self.radii:
                if color_fun == None:
                    plot = ax.plot(self.z_over_r, self.val[r]/yunit, color=cmap((r-rmin)/(rmax-rmin)), label=('$r=%.2f$' % r), **kwargs)
                else:
                    plot = ax.plot(self.z_over_r, self.val[r]/yunit, color=color_fun(r), label=('$r=%.2f$' % r), **kwargs)
        else:
            plot = ax.plot(self.z_over_r, self.val[self.radii[0]]/yunit, label=('$r=%.1f$' % self.radii[0]), **kwargs)
        ax.grid()
        if legend:
            ax.legend()
        if vmin != None or vmax != None:
            ax.set_ylim(vmin, vmax)
        if title == 'auto':
            label = self.vars_object.label + ': ' + self.label
            if log_scale:
                label = 'log10 ' + label
            ax.set_title(label)
        elif title != None:
            ax.set_title(title)
        ax.set_xlabel('$z/r$')
        return ax, plot

# visualizes pattern speed, i.e., variable vs phi vs time
class PatternSpeed (Ops):

    def __init__ (self, vars_object, radii, \
            x1min=0., x2min=(0.5*np.pi),    x3min=0., \
            x1max=np.inf, x2max=(0.5*np.pi), x3max=2.*np.pi):
        super().__init__(vars_object, \
            x1min, x2min, x3min, \
            x1max, x2max, x3max)
        self.label = 'PatternSpeed'
        self.radii = radii
        self.phis = []
        self.values = {} # frames
        self.val = {} # output dict
        self.times = [] # holds all the times read so far (for restarts)

    def integrate (self):

        rr = self.data['x1v']
        rr_idxs = [np.max(np.where(rr<=radius)[0]) for radius in self.radii]
        
        for i in range(len(self.radii)):
            idx = rr_idxs[i]
            r = self.radii[i]
            self.values[r] = self.vars_object.quantity(self.data)[:,0,idx]

    def read (self, filename, verbose=False, temp_save_freq=5, temp_save_path=None):
        # turn filename to list if not list (for averaging)
        if isinstance(filename, str):
            filenames = [filename,]
        else: filenames = filename
        # prepare for saves if requested
        if temp_save_path != None:
            import pickle as pkl
        # read the data
        for idx in range(len(filenames)):
            filename = filenames[idx]
            if verbose: print(' - reading %s.. ' % filename.split('/')[-1], flush=True, end='')
            self.values = {}
            super().read(filename)
            time = self.data['Time']
            self.phis = self.data['x3v']

            if time in self.times:
                if verbose: print('already read. Skipping.', flush=True)
                continue
            self.integrate()
            del self.data
            for r in self.radii:
                if r not in self.val.keys():
                    self.val[r] = {}
                self.val[r][time] = cp.deepcopy(self.vars_object.post_profile(self.values[r]))
            self.times.append(time)
            if verbose: print('done.', flush=True)
            del self.values
            # save for restart if requested
            if temp_save_path != None and idx % temp_save_freq == 0:
                with open(temp_save_path, 'wb') as f:
                    pkl.dump(self, f)
        self.times = np.array(self.times)
        self.phis = np.array(self.phis)
        self.data = {'radii':self.radii, 'phis':self.phis, 'val':self.val}
        if verbose: print('All files read.', flush=True)