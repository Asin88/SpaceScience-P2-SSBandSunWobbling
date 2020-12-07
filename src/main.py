# -*- coding: utf-8 -*-
"""
Created on Sun Dec  6 20:25:00 2020

@author: AAKRITI

Space Science Project 2
"""

# === Goal: Explore the movement of Sun around Solar System Barycenter (SSB) ===

"""
According to Kepler's laws, the Earth revolves in an elliptical orbit with the 
Sun at one of it focii. However, applying Newton's mechanics related to centre 
of mass, it is seen that the centre of mass of the total solar system does not
necessarily coincide with the centre of the Sun. Even though the Sun accounts
for nearly 90% of the total mass of the solar system, the mass of the other 
bodies such as planets, moons, comets, and asteroids are enough to 
offset the centre of mass of the solar system. This centre of mass is known as
the Solar System Barycentre (SSB).

Consequently, the bodies of the solar system including the Sun orbit the SSB,
and not the centre of the Sun. 

In this project, we shall see how far away from the Sun the SSB is. We shall 
analyse whether the distance is significant. We shall also plot the movement of
the Sun around the SSB over a certain time period selected randomly. 

We will be using the NASA toolkit SPICE in Python3 using the wrapper library
spiceypy. The plotting will be done using the matplotlib library. 

In order to perform the necessary computations, we shall need some NAIF kernels.
We have created a subfolder called _kernels under /data/external folder to
hold all kernels needed in this and further projects. We are now also adding a 
meta file 'SPICE_kernels_metafile_P2.txt' under /references/ that holds the 
relative paths of all kernels needed in this project. This file will be 
referenced just once in the program to load all necessary kernels at once. 
"""

# **Import modules**
#%%import datetime as dt
import spiceypy
import math
import numpy as np
from matplotlib import pyplot as plt
#%%

# === Main Function ===
"""
To compute the position vectors of SSB wrt Sun, we shall use the spkgps 
function for positional information. The corresponding kernel has been 
downloaded and its path has been added to the kernel metafile.

Function details found at 
https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkgps_c.html 

Function description: Compute the geometric position of a target body relative 
to an observing body.

Detailed_Input
 
targ        is the standard NAIF ID code for a target body. 
 
et          is the epoch (ephemeris time) at which the position 
           of the target body is to be computed. 
 
ref         is the name of the reference frame to 
           which the vectors returned by the routine should 
           be rotated. This may be any frame supported by 
           the CSPICE subroutine sxform_c. 
 
obs         is the standard NAIF ID code for an observing body. 
 
Detailed_Output
 
pos         contains the position of the target 
           body, relative to the observing body. This vector is 
           rotated into the specified reference frame. Units 
           are always km. 
 
lt          is the one-way light time from the observing body 
           to the geometric position of the target body at the 
           specified epoch. 
               
The target object is SSB (NAIF ID = 0) and the observer is the Sun 
(NAIF ID = 10). NAIF IDs are referenced from 
https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html
The reference frame is the standard “Ecliptic Plane” of our planet set in the 
year 2000. 
"""
"""
The radius of the Sun can be found from the bodvcd function.

Abstract
     Fetch from the kernel pool the double precision values 
     of an item associated with a body, where the body is 
     specified by an integer ID code.
     
Brief_I/O

VARIABLE  I/O  DESCRIPTION
--------  ---  --------------------------------------------------
BODYID     I   Body ID code.
ITEM       I   Item for which values are desired. ('RADII',
                'NUT_PREC_ANGLES', etc. )
MAXN       I   Maximum number of values that may be returned.
DIM        O   Number of values returned.
VALUES     O   Values.

Detailed_Input

BODYID     is the NAIF integer ID code for a body of interest.
            For example, if the body is the earth, the code is
            399.

ITEM       is the item to be returned. Together, the NAIF ID
            code of the body and the item name combine to form a
            kernel variable name, e.g.,

                  'BODY599_RADII'     
                  'BODY401_POLE_RA' 

            The values associated with the kernel variable having
            the name constructed as shown are sought.  Below
            we'll take the shortcut of calling this kernel variable
            the "requested kernel variable."

            Note that ITEM *is* case-sensitive.  This attribute
            is inherited from the case-sensitivity of kernel
            variable names.

MAXN       is the maximum number of values that may be returned.
            The output array VALUES must be declared with size at
            least MAXN.  It's an error to supply an output array
            that is too small to hold all of the values associated
            with the requested kernel variable.

Detailed_Output

DIM        is the number of values returned; this is always the
            number of values associated with the requested kernel
            variable unless an error has been signaled.

VALUES     is the array of values associated with the requested
            kernel variable.  If VALUES is too small to hold all
            of the values associated with the kernel variable, the
            returned values of DIM and VALUES are undefined.
                
We extract the Sun radii (x, y, z components of the Sun ellipsoid) and 
use the only x component as the radius is same along all coordinates.
"""

def main():
    #%%
    #Unload all kernels
    kernel_meta = r"E:\Data Science Projects\Space Science\SpaceScience-P2-SSBandSunWobbling\references\SPICE_kernels_metafile_P2.txt"
    spiceypy.unload( kernel_meta )
    #%%
    #%%
    # === Load the necessary kernels from meta file ===
    spiceypy.kclear()
    spiceypy.furnsh(kernel_meta)
    spiceypy.furnsh(r'E:\Data Science Projects\Space Science\SpaceScience-P2-SSBandSunWobbling\data\external\_kernels\lsk\naif0012.tls')
    #%%
    #%%
    #List loaded kernels
    count = spiceypy.ktotal( 'ALL' ) 
    for i in range(0, count):
        [ file, type, source, handle] = spiceypy.kdata(i, 'ALL');
        print( 'File   {0}'.format(file) )
        print( 'Type   {0}'.format(type) )
        print( 'Source {0}\n'.format(source) )
    #%%
    
    # === Compute position of SSB wrt Sun ===
    
    #First, we compute the position vectors of the SSB wrt the Sun at a
    #randomly chosen time. Let's choose the time as my birthday (24 April 1988,
    #9AM). 
    #%%
    #We shall first convert this time from UTC to ET using the utc2et function.
    init_time_utc = dt.datetime(year=1988,month=4,day=24,hour=9,minute=0,second=0)
    #Convert datetime object to str
    init_time_utc_str = init_time_utc.strftime('%Y-%m-%dT%H:%M:%S')
    #Convert UTC to ET
    init_time_et = spiceypy.utc2et(init_time_utc_str)
    #%%
    #%%
    #Compute the position of SSB wrt Sun
    SSB_position_init, _ = spiceypy.spkgps(targ=0,et=init_time_et,\
                                            ref="ECLIPJ2000",obs=10)
    #Report position
    print(f"On {init_time_utc_str}: \n\nPosition of SSB wrt Sun: \n \
          X-axis = {SSB_position_init[0]} km \n \
          Y-axis = {SSB_position_init[1]} km \n \
          Z-axis = {SSB_position_init[2]} km")
    
    #Compute distance of SSB from Sun in km
    SSB_sun_distance = math.sqrt(SSB_position_init[0]**2.0 \
                                 + SSB_position_init[1]**2.0 \
                                 + SSB_position_init[2]**2.0)
    #Convert km to AU
    SSB_sun_distance_AU = spiceypy.convrt(SSB_sun_distance, 'km', 'AU')
    #Report value
    print(f"Distance between SSB and Sun in km: {SSB_sun_distance} km")
    print(f"Distance between SSB and Sun in AU: {SSB_sun_distance_AU} AU")
    #%%
    # === Plot the distance between the Sun and SSB ===
    
    #To visualise the SSB wrt Sun, we need to first scale the distance to the
    #radius of the Sun using bodvcd function. 
    #%%
    _, radii_sun = spiceypy.bodvcd(bodyid=10, item='RADII', maxn=3)
    radius_sun = radii_sun[0]
    
    #Scale distance between SSB and Sun to radius of Sun
    SSB_position_init_scaled = SSB_position_init/radius_sun
    #%%
    #Plot using matplotlib
    # We only plot the x, y components (view on the ecliptic plane)
    #%%
    SSB_position_init_scaled_xy = SSB_position_init_scaled[0:2]
    print(SSB_position_init_scaled_xy)
    # Set a dark background... since... space is dark
    plt.style.use('dark_background')
    # Create a figure and ax
    fig, ax = plt.subplots(figsize=(12, 8))
    # Create a yellow circle that represents the Sun, add it to the ax
    sun_circ = plt.Circle((0.0, 0.0), 1.0, color='yellow', alpha=0.8)
    ax.add_artist(sun_circ)
    # Plot the centre of the Sun
    ax.plot(0, 0, 'o', color='magenta')
    # Plot the SSB movement
    ax.plot(SSB_position_init_scaled_xy[0], \
            SSB_position_init_scaled_xy[1], \
            'o', color='blue')
    
    # Set some parameters for the plot, set an equal ratio, set a grid, and set
    # the x and y limits
    ax.set_aspect('equal')
    ax.grid(True, linestyle='dashed', alpha=0.5)
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    # Some labelling
    ax.set_xlabel('X in Sun-Radius')
    ax.set_ylabel('Y in Sun-Radius')
    #Annotate the points
    plt.annotate("Centre of Sun", (0,0.1))
    SSB_label = "SSB ({:.2f},{:.2f})".format(SSB_position_init_scaled_xy[0],SSB_position_init_scaled_xy[1])
    plt.annotate(SSB_label, (SSB_position_init_scaled_xy[0], SSB_position_init_scaled_xy[1]+0.1))
    plt.title(f"Position of SSB wrt Centre of Sun on {init_time_utc}")
    #Set font sizes
    small_font_size = 8
    med_font_size = 10
    big_font_size = 12
    plt.rc('font', size=small_font_size)
    plt.rc('axes', titlesize=med_font_size)     
    plt.rc('axes', labelsize=med_font_size)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=small_font_size)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=small_font_size)    # fontsize of the tick labels
    plt.rc('figure', titlesize=big_font_size)  # fontsize of the figure title
    # Saving the figure in high quality
    plt.savefig(r'E:\Data Science Projects\Space Science\SpaceScience-P2-SSBandSunWobbling\reports\figures\SSB_WRT_SUN.png', dpi=300)
    plt.show()
    #%%

# === Run the Script ===
    
if __name__ == "__main__":
    #%%
    #Output file
    outfile = open(r"E:\Data Science Projects\Space Science\SpaceScience-P2-SSBandSunWobbling\reports\outfile.txt","w+")
    #Call the main function
    main()
    #Close output file
    outfile.close()
    #%%