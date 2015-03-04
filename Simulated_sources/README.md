#Simulated Source Code

The beginnings of python code to generate sources in the SpIES mosaics to calculate completeness.
The code is documented, so have a look at the file before running, however, you will need:

1. A mosaic.fits file (which you will need to adjust the name in the python file)
 
2. The empty version of the position.ascii file (which you will need to clear before running as well)

 - Important note: The first line in the position.ascii file should ALWAYS be as follows::
   # RA DEC MAG FLUX  Any different will result in an odd output

 

What does the code do right now?

-At the moment, the code will:

  -Generate random centroid positions

  -Generate random amplitude for the gaussian

  -Evaluate the 2D Gaussian function out to a specified radius

  -Calculate the total flux in the Gaussian (summing the values from the previous step)

  -Convert that to AB mags

  -Write a file containing the random RA DEC MAG and FLUX values

  -Output a new file which is the sum of the input mosaic file and the matrix created with the gaussian sources(all other cells are 0)
  
