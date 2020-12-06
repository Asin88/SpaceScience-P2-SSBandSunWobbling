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

# === Main Function ===

def main():
    
    
# === Run the Script ===
    
if __name__ == "__main__":
    #Output file
    outfile = open("E:\Data Science Projects\Space Science\SpaceScience-P2-SSBandSunWobbling\reports\outfile.txt","w+")
    #Call the main function
    main()
    #Close output file
    outfile.close()

