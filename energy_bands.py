import numpy as np
import matplotlib.pyplot as plt
import os

from sys import version_info
if version_info[0] < 3:
    import pyface.qt
from mayavi import mlab
from mlab_latex_out import mlab_imshow_latex


# Colors for plots c0 for band 0, c1 for band 1 etc.
c0 = (0.0,0.0,0.0)
c1 = (1.0,0.0,0.0)
c2 = (0.0,1.0,0.0)
c3 = (0.0, 0.0, 1.0)   


def show_plot_3D_omega(filename,transparency = False):
    '''
    Shows 3D plot of all bands. Can be saved using the if needed.
    '''
    dirname = 'dispersion_data'
    if not os.path.exists(dirname):
        print('No data generated yet.')
    else:
        if not os.path.exists(dirname+'/'+filename):
            print('"{}/{}" not found.'.format(dirname,filename))
        else:
            data = np.loadtxt(dirname+'/'+filename)
            k1_values = np.unique(data[:,0])
            k2_values = np.unique(data[:,1])
            len_k1 = len(k1_values)
            len_k2 = len(k2_values)
            energy_b0 = data[:,2].reshape((len_k1,len_k2))
            energy_b1 = data[:,3].reshape((len_k1,len_k2))
            energy_b2 = data[:,4].reshape((len_k1,len_k2))
            energy_b3 = data[:,5].reshape((len_k1,len_k2))        

            fig = mlab.figure()
            #ax  = fig.gca(projection='3d')

            band0 = mlab.surf(k1_values, k2_values, energy_b0,color=c0)
            band1 = mlab.surf(k1_values, k2_values, energy_b1,color=c1)
            band2 = mlab.surf(k1_values, k2_values, energy_b2,color=c2)
            band3 = mlab.surf(k1_values, k2_values, energy_b3,color=c3)

            energy_min = min(np.unique([energy_b0,energy_b1,energy_b2,energy_b3]))
            energy_max = max(np.unique([energy_b0,energy_b1,energy_b2,energy_b3]))
            ax_ranges = [np.min(k1_values), np.max(k1_values), np.min(k2_values), np.max(k2_values), energy_min, energy_max]
            ax_scale = [1.0, 1.0, 0.4]
            ax_extent = ax_ranges * np.repeat(ax_scale, 2)

            band0.actor.actor.scale = ax_scale
            band1.actor.actor.scale = ax_scale
            band2.actor.actor.scale = ax_scale
            band3.actor.actor.scale = ax_scale

            latex_obj = mlab_imshow_latex('\omega')

            mlab.view(azimuth=225, elevation=70, distance=15, focalpoint=[0., 0., 0.])
            mlab.outline(band0, color=(.7, .7, .7), extent=ax_extent)
            mlab.axes(band0, color=(.7, .7, .7), extent=ax_extent,
                      ranges=ax_ranges,
                      xlabel='k1', ylabel='k2', zlabel=latex_obj)

            if transparency:    
                band0.actor.property.opacity = 0.5
                band1.actor.property.opacity = 0.5
                band2.actor.property.opacity = 0.5
                band3.actor.property.opacity = 0.5
                
                fig.scene.renderer.use_depth_peeling = 1

            mlab.show()


def plot_2D_omega_from_3D_data(filename,k1,k2):
    '''
    Plots and saves 2D plot of all bands vs k1/k2 for one value of k2/k1
    '''
    dirname = 'dispersion_data'
    if not os.path.exists(dirname):
        print('No data generated yet.')
    else:
        if not os.path.exists(dirname+'/'+filename):
            print('"{}/{}" not found.'.format(dirname,filename))
        else:
            data = np.loadtxt(dirname+'/'+filename)
            k1_values = np.unique(data[:,0])
            k2_values = np.unique(data[:,1])
            len_k1 = len(k1_values)
            len_k2 = len(k2_values)
            energy_b0 = data[:,2].reshape((len_k1,len_k2))
            energy_b1 = data[:,3].reshape((len_k1,len_k2))
            energy_b2 = data[:,4].reshape((len_k1,len_k2))
            energy_b3 = data[:,5].reshape((len_k1,len_k2))       

            if type(k1) == str:
                x_label = 'k1'
                k_x_axis = k1_values
                

                k_param_label = 'k2'
                k_param_idx = np.abs(k2_values-k2).argmin()
                k_param = k2_values[k_param_idx]
                
                plot_title = "2D plot at k2 = {:.2f}".format(k_param)

                energy_b0_2d = energy_b0[:,k_param_idx]
                energy_b1_2d = energy_b1[:,k_param_idx]
                energy_b2_2d = energy_b2[:,k_param_idx]
                energy_b3_2d = energy_b3[:,k_param_idx]
                
            else:
                x_label = 'k2'
                k_x_axis = k2_values
                
                k_param_label = 'k1'
                k_param_idx = np.abs(k1_values-k1).argmin()
                k_param = k1_values[k_param_idx]
                
                plot_title = '2D plot at k1 = {:.2f}'.format(k_param)
                
                energy_b0_2d = energy_b0[k_param_idx,:]
                energy_b1_2d = energy_b1[k_param_idx,:]
                energy_b2_2d = energy_b2[k_param_idx,:]
                energy_b3_2d = energy_b3[k_param_idx,:]

            fig = plt.figure()

            band0 = plt.plot(k_x_axis, energy_b0_2d, label = "omega_0", color=c0)
            band1 = plt.plot(k_x_axis, energy_b1_2d, label = "omega_1", color=c1)
            band2 = plt.plot(k_x_axis, energy_b2_2d, label = "omega_2", color=c2)
            band3 = plt.plot(k_x_axis, energy_b3_2d, label = "omega_3", color=c3)

            energy_min = min(np.unique([energy_b0,energy_b1,energy_b2,energy_b3]))
            energy_max = max(np.unique([energy_b0,energy_b1,energy_b2,energy_b3]))
            ax_ranges = [np.min(k1_values), np.max(k1_values), np.min(k2_values), np.max(k2_values), energy_min, energy_max]

            plt.xlabel(x_label)
            plt.ylabel("$\omega$")
            plt.title(plot_title)

            plt.savefig(dirname+'/'+filename.replace('.txt','_2D_vs_{}_at_{}_{:.2f}.png'.format(x_label,k_param_label,k_param)))
            plt.close()
            plt.clf()