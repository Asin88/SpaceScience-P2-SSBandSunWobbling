# SpaceScience-P2-SSBandSunWobbling

Space Science Project 2
**Goal: Explore the relation between Sun and the Solar System Barycentre (SSB).**


According to Kepler's laws, the Earth revolves in an elliptical orbit with the 
Sun at one of its focii. However, applying Newton's mechanics related to centre 
of mass, it is seen that the centre of mass of the total solar system does not
necessarily coincide with the centre of the Sun. Even though the Sun accounts
for nearly 90% of the total mass of the solar system, the mass of the other 
bodies such as planets, moons, comets, and asteroids are enough to 
offset the centre of mass of the solar system. This centre of mass is known as
the Solar System Barycentre (SSB).

In this project, we shall see how far away from the Sun the SSB is. We shall 
analyse whether the distance is significant. We shall also plot the 3D movement
of the Sun around the SSB over a certain time period selected randomly. The 
current program is the second version of the previous program. Here, we shall
redo the visualization in 3d coordinates. We will be using the NASA toolkit SPICE 
in Python3 using the wrapper library spiceypy.

First, we compute the position vectors of the SSB wrt the Sun at a
randomly chosen time. Let's choose the time as my birthday (24 April 1988,
9AM). We use the spkgps function in spiceypy for this purpose. The calculated
value is more than 600,000 km. In other words, the SSB is more than 6 lakh km 
away from the centre of the Sun on 24th April 1988 at 9AM. 

That might seem to be a huge distance. But when it comes to astronomical distances,
context is necessary to make it intuitive. So let's compare the value with the 
Sun's dimensions. The radius of the Sun is considered to be approximately 
696,000 km. Thus, on the selected day, although the SSB was not at the centre 
of the Sun, it was not situated outside the Sun either. Thus, when considered 
at the scale of the solar system, the difference in positions of Sun and SSB is
that significant. 

The figure below shows the SSB with respect to centre of Sun on the selected 
day. For an intuitive visualization, we have scaled the distances to the 
radius of the Sun.

![2D SSB wrt Sun on 24April1988 9AM]( https://github.com/Asin88/SpaceScience-P2-SSBandSunWobbling/blob/main/reports/figures/SSB_WRT_SUN_INIT.png?raw=true)

Let's look at it in 3d as well. 

![3D SSB wrt Sun on 24April1988 9AM](https://github.com/Asin88/SpaceScience-P2-SSBandSunWobbling/blob/main/reports/figures/SSB_WRT_SUN_INIT_3D.png?raw=true)

So we calculated position of SSB on one particular moment in time and the difference
from position of Sun's centre was noticed. What about other days? Since the position 
of planets around the Sun keeps changing, the mass distribution of the solar system 
also varies. This results in shifts in SSB. Let's have a look at how much does the 
SSB shift over a time period of, say, 30 years. Again using spkgps and matplotlib,
we plot the daily position of SSB wrt centre of Sun after scaling the distances to 
radius of Sun. 

![2D SSB wrt Sun over 30 years](https://github.com/Asin88/SpaceScience-P2-SSBandSunWobbling/blob/main/reports/figures/SSB_WRT_SUN_30_YEARS.png?raw=true)

Here is a 3D rendering of the same figure.

![3D SSB wrt Sun over 30 years](https://github.com/Asin88/SpaceScience-P2-SSBandSunWobbling/blob/main/reports/figures/SSB_WRT_SUN_30_YEARS_3D.png?raw=true)

Isn't it dramatic! The changes in position of SSB is amazing, especially considering 
the fact that it is effect of the mass of the planets and bodies other than the Sun. 
It is observed that the SSB stays outside the Sun about 53% of the time n these 30 
years. This means that in 30 years since my birth, half the time the centre of the 
solar system was not the Sun! That's a total of more than 15 years (almost half my 
lifetime). 

In fact, the closest the SSB came to the centre of the Sun was almost exactly two 
years after the initial time (just one day short). It was only about 44,000 km away. 
The SSB has moving in a sort of spiral trajectory on or near the Sun. 9 years since 
the initial time, it was the farthest away (more than 1 lakh km) in 30 years. 

The importance of this observation is revealed when we consider the fact that 
everything in the solar system revolves around the overall centre of mass. This is 
an extension of the two-body problem where two bodies affected by each other's
gravity revolve around the common centre of mass. Our project has shown that the SSB
is not necessarily the centre of the Sun and may even be situated outside the Sun. 
Thus, technically the earth revolves around not the Sun but the SSB. This also 
means that even the Sun revolves around the SSB. Half the time, when SSB is inside 
the Sun, the Sun was sort of wobbling while spinning. The other half of the time, it
was actually orbiting the SSB. 

The radius of the orbit is small in an astronomical sense, which is why it is often 
ignored in simple calculations. But when performing scientific observations or 
experiments which are impacted by the position of the Sun, even minor changes might 
be significant. 

Moreover, this wobbling movement has been used to detect the presence of planets around 
other stars as well. When a star has a planetary system, it wobbles just like our Sun.
The wobbling is detected through its light spectrum that shifts between blue and red 
colours due to the minor differences in its distance from us. 

Let us visualise the movement of the Sun by making SSB as the fixed reference 
instead of the Sun. The below figure plots the daily position of the centre of the 
Sun wrt SSB over the same 30 years. 

![3D Sun wrt SSB over 30 years](https://github.com/Asin88/SpaceScience-P2-SSBandSunWobbling/blob/main/reports/figures/SUN_WRT_SSB_30_YEARS_3D.png?raw=true)

This project has brought some astounding observations to our view. It has taken us a 
step deeper into the science behind celestial movements. In the next update, we shall
increase our time period to include the past, present and future positions of the Sun. 
We shall also attempt a visualization of the wobbling Sun. Later we will see how much
The Earth's orbit is affected by changing our reference point from the Sun to SSB. 
