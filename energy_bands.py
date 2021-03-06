import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import ntpath

from sys import version_info
if version_info[0] < 3:
    import pyface.qt
from mayavi import mlab
from mlab_latex_out import mlab_imshow_latex

#TODO: integrate LaTeX-formulas in axis- and plotlabels

# Colors for plots c0 for band 0, c1 for band 1 etc.
c0 = (0.0,0.0,0.0)
c1 = (1.0,0.0,0.0)
c2 = (0.0,1.0,0.0)
c3 = (0.0, 0.0, 1.0)   



def show_plot_3D_omega(filename,transparency = False,energy_min=None,
                                                     energy_max=None):
    '''
    Shows 3D plot of all bands. Can be saved using the if needed.
    '''
    exec_file_path = sys.argv[0].replace('.py','')
    exec_filename = ntpath.basename(exec_file_path)
    exec_dir_path = exec_file_path.replace(exec_filename,'')
    dirname = exec_dir_path+'dispersion_data'
    
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

            if energy_max == None:
                energy_max = max(np.unique([energy_b0,energy_b1,energy_b2,energy_b3]))
            if energy_min == None:
                energy_min = min(np.unique([energy_b0,energy_b1,energy_b2,energy_b3]))

            ax_ranges = [np.min(k1_values), np.max(k1_values), np.min(k2_values), np.max(k2_values), energy_min, energy_max]
            ax_scale = [1.0, 1.0, 0.4]
            ax_extent = ax_ranges * np.repeat(ax_scale, 2)

            band0.actor.actor.scale = ax_scale
            band1.actor.actor.scale = ax_scale
            band2.actor.actor.scale = ax_scale
            band3.actor.actor.scale = ax_scale

            mlab.view(azimuth=225, elevation=70, distance=15, focalpoint=[0., 0., 0.])
            mlab.outline(band0, color=(.7, .7, .7), extent=ax_extent)
            mlab.axes(band0, color=(.7, .7, .7), extent=ax_extent,
                      ranges=ax_ranges,
                      xlabel='wave vector k1', ylabel='wave vector k2', zlabel='Energy [meV]')

            if transparency:    
                band0.actor.property.opacity = 0.5
                band1.actor.property.opacity = 0.5
                band2.actor.property.opacity = 0.5
                band3.actor.property.opacity = 0.5
                
                fig.scene.renderer.use_depth_peeling = 1

            mlab.show()


def plot_2D_omega_from_3D_data(filename,k1,k2,energy_min=None,
                                              energy_max=None):
    '''
    Plots and saves 2D plot of all bands vs k1/k2 for one value of k2/k1
    '''
    exec_file_path = sys.argv[0].replace('.py','')
    exec_filename = ntpath.basename(exec_file_path)
    exec_dir_path = exec_file_path.replace(exec_filename,'')
    dirname = exec_dir_path+'dispersion_data'
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
                k_x = 'k1'
                k_x_label = 'wave vector $k_1$'
                k_x_axis = k1_values
                

                k_param = 'k2'
                k_param_value = k2
                k_param_idx = np.abs(k2_values-k2).argmin()
                k_param_value = k2_values[k_param_idx]
                
                plot_title = "2D plot at $k_2$ = {:.2f}".format(k_param_value)

                energy_b0_2d = energy_b0[:,k_param_idx]
                energy_b1_2d = energy_b1[:,k_param_idx]
                energy_b2_2d = energy_b2[:,k_param_idx]
                energy_b3_2d = energy_b3[:,k_param_idx]
                
            else:
                k_x = 'k2'
                k_x_label = 'wave vector $k_2$'
                k_x_axis = k2_values
                
                k_param = 'k1'
                k_param_value = k1
                k_param_idx = np.abs(k1_values-k1).argmin()
                k_param_value = k1_values[k_param_idx]
                
                plot_title = '2D plot at $k_1$ = {:.2f}'.format(k_param_value)
                
                energy_b0_2d = energy_b0[k_param_idx,:]
                energy_b1_2d = energy_b1[k_param_idx,:]
                energy_b2_2d = energy_b2[k_param_idx,:]
                energy_b3_2d = energy_b3[k_param_idx,:]

            fig = plt.figure()

            if energy_max == None:
                energy_max = max(np.unique([energy_b0,energy_b1,energy_b2,energy_b3]))
            if energy_min == None:
                energy_min = min(np.unique([energy_b0,energy_b1,energy_b2,energy_b3]))

            plt.xlim(np.min(k_x_axis),np.max(k_x_axis))
            plt.ylim(energy_min,energy_max)

            band0 = plt.plot(k_x_axis, energy_b0_2d, label = "Band 0", color=c0)
            band1 = plt.plot(k_x_axis, energy_b1_2d, label = "Band 1", color=c1)
            band2 = plt.plot(k_x_axis, energy_b2_2d, label = "Band 2", color=c2)
            band3 = plt.plot(k_x_axis, energy_b3_2d, label = "Band 3", color=c3)

            plt.xlabel(k_x_label)
            plt.ylabel("Energy [meV]")
            plt.legend()
            plt.title(plot_title)
            plt.savefig(dirname+'/'+filename.replace('.txt','_2D_vs_{}_at_{}_{:.2f}.png'.format(k_x,k_param,k_param_value)))
            plt.close()
            plt.clf()


def plot_2D_omega_path_k1_k2_space(filename,energy_min=None, energy_max=None):
    '''
    Plot path from Gamma (0,0) to X (pi,0) to M (pi,pi) to Gamma(0,0)
    '''
    exec_file_path = sys.argv[0].replace('.py','')
    exec_filename = ntpath.basename(exec_file_path)
    exec_dir_path = exec_file_path.replace(exec_filename,'')
    dirname = exec_dir_path+'dispersion_data'
    
    if not os.path.exists(dirname):
        print('No data generated yet.')
    else:
        if not os.path.exists(dirname+'/'+filename):
            print('"{}/{}" not found.'.format(dirname,filename))
        else:
            data = np.loadtxt(dirname+'/'+filename)
            k1_values = data[:,0]
            k2_values = data[:,1]

            k1_un_values = np.unique(k1_values)
            k2_un_values = np.unique(k2_values)

            len_k1 = len(k1_un_values)
            len_k2 = len(k2_un_values)

            energy_b0 = data[:,2].reshape((len_k1,len_k2))
            energy_b1 = data[:,3].reshape((len_k1,len_k2))
            energy_b2 = data[:,4].reshape((len_k1,len_k2))
            energy_b3 = data[:,5].reshape((len_k1,len_k2))            

            idx_start_pt = 0
            idx_end_pt = len_k1-1


            k_path_gamma_X = k1_un_values - k2_un_values[idx_start_pt]
            k_path_X_M = k1_un_values[idx_end_pt] - k2_un_values
            k_path_M_gamma = np.flip(k1_un_values) - np.flip(k2_un_values)

            k_path_tot = np.concatenate((k_path_gamma_X, k_path_X_M, k_path_M_gamma))

            energy_b0_path_gamma_X = energy_b0[:,idx_start_pt]
            energy_b0_path_X_M = energy_b0[idx_end_pt,:]

            energy_b0_path_M_gamma = np.zeros(len_k1)
            for i in np.arange(len_k1):
                energy_b0_path_M_gamma[i] = energy_b0[i,i]

            energy_b0_path_tot = np.concatenate((energy_b0_path_gamma_X,
                                                 energy_b0_path_X_M,
                                                 np.flip(energy_b0_path_M_gamma)))

            energy_b1_path_gamma_X = energy_b1[:,idx_start_pt]
            energy_b1_path_X_M = energy_b1[idx_end_pt,:]

            energy_b1_path_M_gamma = np.zeros(len_k1)
            for i in np.arange(len_k1):
                energy_b1_path_M_gamma[i] = energy_b1[i,i]

            energy_b1_path_tot = np.concatenate((energy_b1_path_gamma_X,
                                                 energy_b1_path_X_M,
                                                 np.flip(energy_b1_path_M_gamma)))

            energy_b2_path_gamma_X = energy_b2[:,idx_start_pt]
            energy_b2_path_X_M = energy_b2[idx_end_pt,:]

            energy_b2_path_M_gamma = np.zeros(len_k1)
            for i in np.arange(len_k1):
                energy_b2_path_M_gamma[i] = energy_b2[i,i]

            energy_b2_path_tot = np.concatenate((energy_b2_path_gamma_X,
                                                 energy_b2_path_X_M,
                                                 np.flip(energy_b2_path_M_gamma)))

            energy_b3_path_gamma_X = energy_b3[:,idx_start_pt]
            energy_b3_path_X_M = energy_b3[idx_end_pt,:]

            energy_b3_path_M_gamma = np.zeros(len_k1)
            for i in np.arange(len_k1):
                energy_b3_path_M_gamma[i] = energy_b3[i,i]

            energy_b3_path_tot = np.concatenate((energy_b3_path_gamma_X,
                                                 energy_b3_path_X_M,
                                                 np.flip(energy_b3_path_M_gamma)))

            fig, ax = plt.subplots(1,1)
            x = np.arange(len(k_path_tot))
            ax.set_xticks(x)

            print(x)
            xxxlabels = ['' for kk in k_path_tot]
            xxxlabels[0] = "$\Gamma$"
            xxxlabels[-1] = "$\Gamma$"
            xxxlabels[len_k1] = "X"
            xxxlabels[2*len_k1] = "M"

            print(k_path_tot[len_k1])
            print(k_path_tot[2*len_k1-1])


            print(xxxlabels)
            ax.set_xticklabels(xxxlabels)

            if energy_max == None:
                energy_max = max(np.unique([energy_b0,energy_b1,energy_b2,energy_b3]))
            if energy_min == None:
                energy_min = min(np.unique([energy_b0,energy_b1,energy_b2,energy_b3]))

            plt.ylim(energy_min,energy_max)

            band0 = ax.plot(x, energy_b0_path_tot, label = "Band 0", color=c0)
            band1 = ax.plot(x, energy_b1_path_tot, label = "Band 1", color=c1)
            band2 = ax.plot(x, energy_b2_path_tot, label = "Band 2", color=c2)
            band3 = ax.plot(x, energy_b3_path_tot, label = "Band 3", color=c3)

            plt.xlabel("Path in k1-k2-space")
            plt.ylabel("Energy [meV]")
            plt.legend()
            fig.canvas.draw()
            plt.savefig(dirname+'/'+filename.replace('.txt','_2D_path_k1-k2-space.png'))
            plt.close()
            plt.clf()