#=================================================================================
# 
# Plot vertical localization functions
#
# Original code written by Chris Reidel
#
# Modified by Lou Wicker Nov 22, 2024
#
# Lib dependencies:  numpy / matplotlib / scipy / os / OptionParser / ambiance
#------------------
#
# Note: ambiance library is available from either 'pip' and 'conda'
#
# Note:  vertical localizations are specified in "half-width" or "DART" format.
#
# Note:  To see how tp run the code, type "python plot_vloc.py -h"
#
#=================================================================================

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, NullFormatter, ScalarFormatter)
from ambiance import Atmosphere
from scipy import interpolate
import sys
import os
from optparse import OptionParser

heights    = np.arange(0e3, 30e3+100,100)[1:]
atmosphere = Atmosphere(heights)
pressures  = atmosphere.pressure/100.
fpress     = interpolate.interp1d(pressures*100.,heights)

cutoff = 1.0   # this might be useful to change....

#---------------------------------------------------------------------------------------------------
#
def compute_localization(coord='ScaleHght', ob_level_press=700., vert_norm=0.4):

    def gc_lt1(r):
        return (((-0.25*r +0.5)*r + 0.625)*r - (5.0/3.0))*r**2 + 1.0
        
    def gc_lt2(r):
        return (((((r/12.0) - 0.5)*r +0.625)*r + (5.0/3.0))*r -5.0)*r + 4.0 - 2.0 / (3.0 * r)
       
    vert_coor = coord    
    obs_level_press = ob_level_press*100.
    
#     print(' =======================================')
    
#     print(' Input observation pressure level'.format(obs_level_press))
#     
#     print(' Coordinate chosen:  {0}'.format(coord))
#     
#     print(' Localization Length:  {0}\n'.format(vert_norm))

    if vert_coor == 'ScaleHght':
      print(' Log(p) of Ob:  {0}'.format(np.log(obs_level_press)))
      location_diff = (-np.log(pressures*100) - -np.log(obs_level_press))/vert_norm
      
    elif vert_coor == 'Height':
      obs_level_hght = fpress(obs_level_press)
      print('Height of Ob:  {0}'.format(obs_level_hght))
      location_diff = (heights - obs_level_hght)/vert_norm
      
    elif vert_coor == 'Pressure':
      location_diff = (pressures*100 - obs_level_press)/(vert_norm*100)
            
    local = np.zeros((location_diff.shape[0]))
    
    location_diff  = np.abs( location_diff )
    
    r = location_diff / cutoff
    
    local = np.where(r < 2.0*cutoff, gc_lt2(r), local)

    local = np.where(r <     cutoff, gc_lt1(r), local)
            
    return local

#---------------------------------------------------------------------------------------------------
# Main function defined to return correct sys.exit() calls
#

if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option("-c", "--coord", dest="coord", type="string", default= "ScaleHght", \
                                  help="type of vertical localizaton:  Height, Pressure, ScaleHeight [default]")
    parser.add_option("-l", "--loc", dest="loc",   type="float",  default=0.4,  \
                                  help = "Vertical Localization length in Z, P, or ScaleHeight")

    parser.add_option("-o", "--ob_p", dest="ob_level_press",   type="float",  default=700.,  \
                                  help = "Pressure level (mb) of ob location")

    (options, args) = parser.parse_args()
        
    print(' =======================================================')
    print(' Input:  Type of Vertical Localization: {0}'.format(options.coord))
    print(' Input:  Half-width of Vertical Localization: {0}'.format(options.loc))
    print(' Input:  Pressure level of observation: {0} mb\n'.format(options.ob_level_press))

    local = compute_localization(coord=options.coord, vert_norm=options.loc, ob_level_press=options.ob_level_press)

    print(local.shape)
        
    ticks = []
    start = 1000.
    while start > 10:
      ticks.append(int(start))
      start = start - 50.
    ticks.append(int(10))

    ticks = np.array([10,50,100,150,175,200,225,250,275,300,325,350,375,400,450,500,600,700,800,900,1000])[::-1]

    f = interpolate.interp1d(pressures,heights/1000.)
    ynew = f(ticks[:-1])
    plt.rc('font', weight='bold')
    fig,ax = plt.subplots(1,1,figsize=(8,9))
    ax.fill_betweenx(pressures,local,color='gray',zorder=10)

    local2 = np.where(local > 0.5, local, 0.0)
    ax.fill_betweenx(pressures,local2,color='#539ecd',zorder=10)

    ax.plot(local2,pressures,'-k',zorder=10)

    ###
    ax2 = ax.twinx()
    #ax.grid(True,zorder=1,axis='y')
    ax2.grid(True,axis='both',linewidth=0.5,linestyle='dashed',zorder=1.0)
#   ax2.plot(local,pressures,'-k',alpha=1.0,zorder=10)
#   ax2.fill_betweenx(pressures,local,color='#539ecd',zorder=10)
    ax2.set_yscale("log")
    ax2.yaxis.set_major_formatter(ScalarFormatter())
    ax2.yaxis.set_minor_formatter(NullFormatter())
    ax2.set_yticks(ticks[:-1])
    ax2.set_yticklabels(np.round(ynew,2))
    ax2.set_ylim(pressures[0], 50)
    ####
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_yticks(ticks)
    ax.set_yticklabels(ticks)
    ax.set_ylim(pressures[0], 50)
    ax.set_ylabel('Pressure (hPa)',weight='bold')
    ax2.set_ylabel('Height (KM)',weight='bold')
    ax.set_xlim(-0.008,1.02)
    ax.set_xticks(np.arange(0,1+0.1,0.1))
    
    if options.coord == 'ScaleHght':
      ax.set_title('Scale Height Vertical Localization Standard Atmosphere\nCutoff: {0} - vert_normal: {1}'.format(cutoff,options.loc),weight='bold')
    elif options.coord == 'Pressure':
      ax.set_title('Pressure Vertical Localization Standard Atmosphere\nCutoff: {0} - vert_normal: {1}'.format(cutoff,options.loc),weight='bold')
    elif options.coord == 'Height':
      ax.set_title('Height Vertical Localization Standard Atmosphere\nCutoff: {0} - vert_normal: {1}'.format(cutoff,options.loc),weight='bold')
    
    print(' Output:  Saving PNG file:  VLoc_{0}.png'.format(options.coord))
    plt.savefig('VLoc_{0}.png'.format(options.coord),dpi=500,bbox_inches='tight')
    
    os.system('open VLoc_%s.png' % options.coord)
    
    
    
    
    
    #fig = plt.figure()
    #plt.plot(dist*6378.1,local,'-k')
    #plt.fill_between(dist*6378.1,local, color='#539ecd')
    #plt.plot(np.array([dist[480],dist[480]])*6378.1,np.array([0,local[480]]),'--k')
    #plt.plot(np.array([dist[160],dist[160]])*6378.1,np.array([0,local[160]]),'--k')
    #plt.grid(True)
    #plt.ylim(-0.008,1.02)
    #plt.ylabel('Localization Factor',name='Calibri',size=12,weight='bold')
    #plt.xlabel('Distance (Km)',name='Calibri',size=12,weight='bold')
    #plt.savefig('horizontal_localization.eps',dpi=500,bbox_inches='tight')


    #fig = plt.figure()
    #plt.plot(local,dist*80.,'-k')
    #plt.fill_between(local,dist*80., color='#539ecd')
    ##plt.plot(dist*6.25,local,'-k')
    ##plt.fill_between(dist*6.25,local,where=np.abs(dist*6.25)<1., color='#539ecd')
    #plt.plot(np.array([0,local[480]]),np.array([dist[480],dist[480]])*80.,'--k')
    #plt.plot(np.array([0,local[160]]),np.array([dist[160],dist[160]])*80.,'--k')
    #plt.ylabel('Distance (Km)',name='Calibri',size=12,weight='bold')
    #plt.xlabel('Localization Factor',name='Calibri',size=12,weight='bold')
    #plt.xlim(-0.008,1.02)
    #
    #plt.grid(True)
    #plt.savefig('vertical_localization.eps',dpi=500,bbox_inches='tight')


    #local_vert.py
    #Displaying local_vert.py.
