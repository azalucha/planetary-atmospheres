# planetary-atmospheres

This code represents the decade of work I did developing subroutines for the MITgcm (http://mitgcm.org) to model planetary atmospheres.  The code is written in Fortran 77 and is meant to be plugged in to the MITgcm dynamical core (last version that I ran it with was checkpoint65x).  Simply place these subroutines in your "code" directory.  This code is freely distributed AS IS.  I might have been in the middle of editing a subroutine when I ceased researching planetary atmospheres, so use with caution.  Some modules have been deliberately excluded, namely, radcodes 9, 10, and 13 (see AMZVARS.h for a description of the variables I added) because they are not mine to distribute.  Radcode 9 (NASA Ames Mars atmosphere radiation code) is available at https://spacescience.arc.nasa.gov/mars-climate-modeling-group/code_form.html.  Radcodes 10 and 13 (Pluto) is available from Darrell Strobel (JHU) and/or Xun Zhu (APL) upon request.  Citations for the remaining radiation codes are as follows:

2: Zalucha, A. M., R. A. Plumb, and R. J. Wilson, "An Analysis of the Effect of Topography on the Martian Hadley Cells," Journal of the Atmospheric Sciences 67, 673, 2010 (doi: 10.1175/2009JAS3130.1).

4, 8: Zalucha, A. M., “Demonstration of a GCM for Mars, GJ 1214b, Pluto, and Triton,” Comparative Climatology of Terrestrial Planets, 2012.

6 (includes Pluto and Triton): Zalucha, A. M. and T. I. Michaels, "A 3D General Circulation Model for Pluto and Triton with Fixed Volatile Abundance and Simplified Surface Forcing," Icarus 223, 819, 2013 (arXiv:1211.0009).

11, 12: I used these in a proposal, but they were never published, so use your best judgement (cite this repository?)
