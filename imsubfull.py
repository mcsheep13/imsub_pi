# ProperImage Code Iteration 2 (add time and background adjustments)

# Imports here
from copy import copy
import numpy as np
from astropy.io import fits
import properimage.single_image as si
from properimage.operations import subtract
import os

# Move to working files directory
os.chdir("D:\\Aligned_Frames\\12June2012\\Aligned")

# Intake files and process (only once for template, loop for all others)
imt = fits.open('template.fits')
expt = 120	# Hard coded exposure time
bgval = np.median(imt[0].data[1000:1100,1070:1170])	# Median of blank space for background
imtemp = (imt[0].data-bgval)/expt
imt[0].data = imtemp
imt.writeto('adj_template.fits')

# Image file processing
n = 0
# for n in range (695):
ind=f'{n:04}'
imt = fits.open('Aligned_12June{}.fits'.format(ind))
expt = 15	# Exposure time coded in
bgval = np.median(imt[0].data[1000:1100,1070:1170])	# Median of blank space for background
imtemp = (imt[0].data-bgval)/expt
imt[0].data = imtemp
imt.writeto('adj_Aligned{}.fits'.format(ind))


# Subtraction time

# Move to correct directory
os.chdir("D:\\Aligned_Frames\\12June2012\\Aligned\\propim_results\\adj_results")

# Set template path
ref_path = "D:\\Aligned_Frames\\12June2012\\Aligned\\adj_template.fits"

# Function here
n = 0
# for n in range(695):
# Set new image path
ind=f'{n:04}'
new_path = "D:\\Aligned_Frames\\12June2012\\Aligned\\adj_Aligned{}.fits".format(ind)

# Process image subtraction
result = subtract(
    ref=ref_path,
    new=new_path,
    smooth_psf=False,
    fitted_psf=True,
    align=False,
    iterative=False,
    beta=False,
    shift=False
)

# Extract data
D = np.real(result[0])
P = np.real(result[1])
Scorr = np.real(result[2])
mask = np.real(result[3])

# Export data to files
im = fits.open(new_path)
primary = fits.PrimaryHDU(header=im[0].header)
file1 = fits.ImageHDU(D)
file2 = fits.ImageHDU(P)
file3 = fits.ImageHDU(Scorr)
file4 = mask
struct1 = fits.HDUList([primary,file1])
struct2 = fits.HDUList([primary,file2])
struct3 = fits.HDUList([primary,file3])
struct1.writeto('adj_D{}.fits'.format(ind))
struct2.writeto('adj_P{}.fits'.format(ind))
struct3.writeto('adj_Scorr{}.fits'.format(ind))
np.savetxt('adj_mask{}.txt'.format(ind), np.where(file4)[0])

# Confirm that code ran successfully
print('Processing finished')