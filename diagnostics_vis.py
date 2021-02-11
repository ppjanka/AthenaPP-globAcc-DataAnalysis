'''
[diagnostics_vis.py]
Plot definitions. Note: spherical_polar coordinates only!
Author: Patryk Pjanka, patryk.pjanka@su.se
'''
import numpy as np

exec(open("./diagnostics_header.py").read())
import diagnostics_ops as diops
import diagnostics_vars as divars

non_vminmax = {'poloidal_profile': None, 'poloidal_slice': None, 'equatorial_slice': None, 'theta_profile': None, 'radial_profile': None}

def process_figure (fig, gs, \
                   athdf_file_left, athdf_file_right, \
                   vars_obj_left, vars_obj_right, \
                   left_vmin=non_vminmax, left_vmax=non_vminmax, \
                   right_vmin=non_vminmax, right_vmax=non_vminmax, \
                   left_label=None, right_label=None, \
                   left_log_scale=True, right_log_scale=True, \
                   left_suptype=None, right_suptype=None, \
                   left_slice_or_prof='slice', right_slice_or_prof='slice', \
                   scalar_mappables=None, all_colorbars=False, \
                   left_rmax=None, right_rmax=None, \
                   left_cmap=None, right_cmap=None, \
                   left_cbar_format=None, right_cbar_format=None, \
                   left_rprofile_fit=False, right_rprofile_fit=False, \
                   markings=True, return_1D_profiles=False, \
                   disk_x2min=0., disk_x2max=np.pi, disk_x1max=np.inf, \
                   theta_H=0.2):
    if return_1D_profiles:
        outputs = {}; outputs[left_label] = {}; outputs[right_label] = {}
    # LEFT: --------------------------------------------------
    athdf_file = athdf_file_left
    vars_obj = vars_obj_left
    if left_cmap != None:
        left_kwargs = {'cmap': left_cmap}
    else:
        left_kwargs = {}
    # poloidal profile
    kwargs = {}
    if left_suptype == 'sum': kwargs['type'] = 'sum'
    if left_rmax != None: kwargs['x1max'] = left_rmax
    ops_obj = diops.Profile_ThetaR(vars_obj, **kwargs)
    ops_obj.read(athdf_file)
    TIME = ops_obj.time
    kwargs = {'subplot_spec':gs[0,1], 'flip':True, 'crop':True, 'log_scale':left_log_scale, 'levels':20, 'vmin':left_vmin['poloidal_profile'], 'vmax':left_vmax['poloidal_profile'], 'title':'auto', 'kwargs':left_kwargs}
    ax, plot = ops_obj.plot(fig, **kwargs)
    ax.set_ylim(0.,left_rmax)
    if markings:
        for angle in [(np.pi-5.*theta_H),(np.pi+5.*theta_H)]:
            ax.vlines(angle, 0., 1., 'r', linewidth=1.5, linestyle=':')
    #colorbar
    ax = fig.add_subplot(gs[0,0])
    if left_vmin['poloidal_profile'] != None and left_vmax['poloidal_profile'] != None and left_vmax['poloidal_profile'] > left_vmin['poloidal_profile'] and scalar_mappables[0] != None:
        m = scalar_mappables[0]
        if left_log_scale:
            m.set_array(np.linspace(np.log10(left_vmin['poloidal_profile']), np.log10(left_vmax['poloidal_profile']),20))
            m.set_clim(np.log10(left_vmin['poloidal_profile']), np.log10(left_vmax['poloidal_profile']))
        else:
            m.set_array(np.linspace(left_vmin['poloidal_profile'], left_vmax['poloidal_profile'],20))
            m.set_clim(left_vmin['poloidal_profile'], left_vmax['poloidal_profile'])
    else:
        m = plot
    fig.colorbar(m, cax=ax, orientation='vertical', format=left_cbar_format)
    del ops_obj
    # poloidal slices
    kwargs = {}
    if left_rmax != None: kwargs['x1max'] = left_rmax
    ops_obj = diops.PoloidalSlice(vars_obj, **kwargs)
    ops_obj.read(athdf_file)
    kwargs = {'subplot_spec':gs[0,2], 'log_scale':left_log_scale, 'levels':20, 'vmin':left_vmin['poloidal_slice'], 'vmax':left_vmax['poloidal_slice'], 'title':'auto', 'kwargs':left_kwargs}
    ax, plot = ops_obj.plot(fig, **kwargs)
    ax.set_ylim(0.,left_rmax)
    if all_colorbars: 
        fig.colorbar(plot, cax=ax, orientation='horizontal')
    del ops_obj
    kwargs = {}
    if left_rmax != None: kwargs['x1max'] = left_rmax
    ops_obj = diops.PoloidalSlice(vars_obj, intersect=np.pi, **kwargs)
    ops_obj.read(athdf_file)
    kwargs = {'ax':ax, 'subplot_spec':gs[0,2], 'flip':True, 'log_scale':left_log_scale, 'levels':20, 'vmin':left_vmin['poloidal_slice'], 'vmax':left_vmax['poloidal_slice'], 'title':'auto', 'kwargs':left_kwargs}
    ax, plot = ops_obj.plot(fig, **kwargs)
    ax.set_ylim(0.,left_rmax)
    if all_colorbars: 
        fig.colorbar(plot, cax=ax, orientation='horizontal')
    if markings:
        for angle in [-5.*theta_H, 5.*theta_H, (np.pi-5.*theta_H), (np.pi+5.*theta_H)]:
            ax.vlines(angle, 0., 1., 'r', linewidth=1.5, linestyle=':')
    del ops_obj
    # equatorial slice
    kwargs = {}
    if left_rmax != None: kwargs['x1max'] = left_rmax
    ops_obj = diops.EquatorialSlice(vars_obj, **kwargs)
    ops_obj.read(athdf_file)
    kwargs = {'subplot_spec':gs[1:3,2], 'log_scale':left_log_scale, 'levels':20, 'vmin':left_vmin['equatorial_slice'], 'vmax':left_vmax['equatorial_slice'], 'title':'auto', 'kwargs':left_kwargs}
    ax, plot = ops_obj.plot(fig, **kwargs)
    ax.set_ylim(0.,left_rmax)
    # time data
    fig.suptitle("Time = %.2f sim.u. = %.2f $P_{orb}$" % (TIME, TIME / orb_period))
    del ops_obj
    # colorbar
    ax = fig.add_subplot(gs[3,2])
    if left_log_scale:
        ax.set_title('log10 '+left_label)
    else:
        ax.set_title(left_label)
    if left_vmin['equatorial_slice'] != None and left_vmax['equatorial_slice'] != None and left_vmax['equatorial_slice'] > left_vmin['equatorial_slice'] and scalar_mappables[1] != None:
        m = scalar_mappables[1]
        if left_log_scale:
            m.set_array(np.linspace(np.log10(left_vmin['equatorial_slice']), np.log10(left_vmax['equatorial_slice']),20))
            m.set_clim(np.log10(left_vmin['equatorial_slice']), np.log10(left_vmax['equatorial_slice']))
        else:
            m.set_array(np.linspace(left_vmin['equatorial_slice'], left_vmax['equatorial_slice'],20))
            m.set_clim(left_vmin['equatorial_slice'], left_vmax['equatorial_slice'])
    else:
        m = plot
    fig.colorbar(m, cax=ax, orientation='horizontal', format=left_cbar_format)
    # theta profile
    kwargs = {}
    if left_suptype == 'sum': kwargs['type'] = 'sum'
    ops_obj = diops.Profile_Theta(vars_obj, x1max=disk_x1max, **kwargs)
    ops_obj.read(athdf_file)
    kwargs = {'subplot_spec':gs[1,0:2], 'log_scale':left_log_scale, 'vertical':True, 'title':'auto', 'vmin':left_vmin['theta_profile'], 'vmax':left_vmax['theta_profile'], 'theta_H':theta_H}
    ax, plot = ops_obj.plot(fig, **kwargs)
    ax.invert_xaxis()
    if return_1D_profiles:
        outputs[left_label]['theta'] = ops_obj.data
    del ops_obj
    # radial profile
    kwargs = {}
    if left_suptype == 'sum': kwargs['type'] = 'sum'
    ops_obj = diops.RadialProfile(vars_obj, x2min=disk_x2min, x2max=disk_x2max, **kwargs)
    ops_obj.read(athdf_file)
    kwargs = {'subplot_spec':gs[2:4,0:2], 'log_scale':left_log_scale, 'title':'auto', 'vmin':left_vmin['radial_profile'], 'vmax':left_vmax['radial_profile'], 'fit':left_rprofile_fit}
    ax, plot = ops_obj.plot(fig, **kwargs)
    ax.invert_xaxis()
    if left_rmax != None:
        ax.axvline(left_rmax, color='r', linewidth=1.5, linestyle=':')
    if return_1D_profiles:
        outputs[left_label]['r'] = ops_obj.data
    del ops_obj
    del vars_obj

    # RIGHT: --------------------------------------------------
    athdf_file = athdf_file_right
    vars_obj = vars_obj_right
    if right_cmap != None:
        right_kwargs = {'cmap': right_cmap}
    else:
        right_kwargs = {}
    # poloidal profile
    kwargs = {}
    if right_suptype == 'sum': kwargs['type'] = 'sum'
    if right_rmax != None: kwargs['x1max'] = right_rmax
    ops_obj = diops.Profile_ThetaR(vars_obj, **kwargs)
    ops_obj.read(athdf_file)
    kwargs = {'subplot_spec':gs[0,4], 'flip':False, 'crop':True, 'log_scale':right_log_scale, 'levels':20, 'vmin':right_vmin['poloidal_profile'], 'vmax':right_vmax['poloidal_profile'], 'title':'auto', 'kwargs':right_kwargs}
    ax, plot = ops_obj.plot(fig, **kwargs)
    ax.set_ylim(0.,right_rmax)
    if markings:
        for angle in [-5.*theta_H, 5.*theta_H]:
            ax.vlines(angle, 0., 1., 'r', linewidth=1.5, linestyle=':')
    #colorbar
    ax = fig.add_subplot(gs[0,5])
    if right_vmin['poloidal_profile'] != None and right_vmax['poloidal_profile'] != None and right_vmax['poloidal_profile'] > right_vmin['poloidal_profile'] and scalar_mappables[2] != None:
        m = scalar_mappables[2]
        if right_log_scale:
            m.set_array(np.linspace(np.log10(right_vmin['poloidal_profile']), np.log10(right_vmax['poloidal_profile']),20))
            m.set_clim(np.log10(right_vmin['poloidal_profile']), np.log10(right_vmax['poloidal_profile']))
        else:
            m.set_array(np.linspace(right_vmin['poloidal_profile'], right_vmax['poloidal_profile'],20))
            m.set_clim(right_vmin['poloidal_profile'], right_vmax['poloidal_profile'])
    else:
        m = plot
    cbar = fig.colorbar(m, cax=ax, orientation='vertical', format=right_cbar_format)
    del ops_obj
    # poloidal slices
    kwargs = {}
    if right_rmax != None: kwargs['x1max'] = right_rmax
    ops_obj = diops.PoloidalSlice(vars_obj, **kwargs)
    ops_obj.read(athdf_file)
    kwargs = {'subplot_spec':gs[0,3], 'log_scale':right_log_scale, 'levels':20, 'vmin':right_vmin['poloidal_slice'], 'vmax':right_vmax['poloidal_slice'], 'title':'auto', 'kwargs':right_kwargs}
    ax, plot = ops_obj.plot(fig, **kwargs)
    ax.set_ylim(0.,right_rmax)
    if all_colorbars: 
        fig.colorbar(plot, cax=ax, orientation='horizontal')
    del ops_obj
    kwargs = {}
    if right_rmax != None: kwargs['x1max'] = right_rmax
    ops_obj = diops.PoloidalSlice(vars_obj, intersect=np.pi, **kwargs)
    ops_obj.read(athdf_file)
    kwargs = {'ax':ax, 'subplot_spec':gs[0,3], 'flip':True, 'log_scale':right_log_scale, 'levels':20, 'vmin':right_vmin['poloidal_slice'], 'vmax':right_vmax['poloidal_slice'], 'title':'auto', 'kwargs':right_kwargs}
    ax, plot = ops_obj.plot(fig, **kwargs)
    ax.set_ylim(0.,right_rmax)
    if all_colorbars: 
        fig.colorbar(plot, cax=ax, orientation='horizontal')
    if markings:
        for angle in [-5.*theta_H, 5.*theta_H, (np.pi-5.*theta_H), (np.pi+5.*theta_H)]:
            ax.vlines(angle, 0., 1., 'r', linewidth=1.5, linestyle=':')
    del ops_obj
    # equatorial slice
    kwargs = {}
    if right_rmax != None: kwargs['x1max'] = right_rmax
    ops_obj = diops.EquatorialSlice(vars_obj, **kwargs)
    ops_obj.read(athdf_file)
    kwargs = {'subplot_spec':gs[1:3,3], 'log_scale':right_log_scale, 'levels':20, 'vmin':right_vmin['equatorial_slice'], 'vmax':right_vmax['equatorial_slice'], 'title':'auto', 'kwargs':right_kwargs}
    ax, plot = ops_obj.plot(fig, **kwargs)
    ax.set_ylim(0.,right_rmax)
    # colorbar
    ax = fig.add_subplot(gs[3,3])
    if right_log_scale:
        ax.set_title('log10 '+right_label)
    else:
        ax.set_title(right_label)
    if right_vmin['equatorial_slice'] != None and right_vmax['equatorial_slice'] != None and right_vmax['equatorial_slice'] > right_vmin['equatorial_slice'] and scalar_mappables[3] != None:
        m = scalar_mappables[3]
        if right_log_scale:
            m.set_array(np.linspace(np.log10(right_vmin['equatorial_slice']), np.log10(right_vmax['equatorial_slice']),20))
            m.set_clim(np.log10(right_vmin['equatorial_slice']), np.log10(right_vmax['equatorial_slice']))
        else:
            m.set_array(np.linspace(right_vmin['equatorial_slice'], right_vmax['equatorial_slice'],20))
            m.set_clim(right_vmin['equatorial_slice'], right_vmax['equatorial_slice'])
    else:
        m = plot
    cbar = fig.colorbar(m, cax=ax, orientation='horizontal', format=right_cbar_format)#, boundaries=ticks)
    del ops_obj
    # theta profile
    kwargs = {}
    if right_suptype == 'sum': kwargs['type'] = 'sum'
    ops_obj = diops.Profile_Theta(vars_obj, x1max=disk_x1max, **kwargs)
    ops_obj.read(athdf_file)
    kwargs = {'subplot_spec':gs[1,4:6], 'log_scale':right_log_scale, 'vertical':True, 'title':'auto', 'vmin':right_vmin['theta_profile'], 'vmax':right_vmax['theta_profile'], 'theta_H':theta_H}
    ax, plot = ops_obj.plot(fig, **kwargs)
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    if return_1D_profiles:
        outputs[right_label]['theta'] = ops_obj.data
    del ops_obj
    # radial profile
    kwargs = {}
    if right_suptype == 'sum': kwargs['type'] = 'sum'
    ops_obj = diops.RadialProfile(vars_obj, x2min=disk_x2min, x2max=disk_x2max, **kwargs)
    ops_obj.read(athdf_file)
    kwargs = {'subplot_spec':gs[2:4,4:6], 'log_scale':right_log_scale, 'title':'auto', 'vmin':right_vmin['radial_profile'], 'vmax':right_vmax['radial_profile'], 'fit':right_rprofile_fit}
    ax, plot = ops_obj.plot(fig, **kwargs)
    if right_rmax != None:
        ax.axvline(right_rmax, color='r', linewidth=1.5, linestyle=':')
    if return_1D_profiles:
        outputs[right_label]['r'] = ops_obj.data
    del ops_obj
    del vars_obj
    #          --------------------------------------------------
    if return_1D_profiles:
        return TIME, outputs

def get_1D_profiles (athdf_file, vars_obj, label, outputs, suptype=None, \
                   disk_x2min=0., disk_x2max=np.pi, disk_x1max=np.inf):

    outputs[label] = {}

    # theta profile
    kwargs = {}
    if suptype == 'sum': kwargs['type'] = 'sum'
    ops_obj = diops.Profile_Theta(vars_obj, x1max=disk_x1max, **kwargs)
    ops_obj.read(athdf_file)
    TIME = ops_obj.time
    outputs[label]['theta'] = ops_obj.data
    del ops_obj
    # radial profile
    kwargs = {}
    if suptype == 'sum': kwargs['type'] = 'sum'
    ops_obj = diops.RadialProfile(vars_obj, x2min=disk_x2min, x2max=disk_x2max, **kwargs)
    ops_obj.read(athdf_file)
    outputs[label]['r'] = ops_obj.data
    del ops_obj
    del vars_obj

    return TIME, outputs