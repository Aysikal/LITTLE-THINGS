from astropy.io import fits
import numpy as np
from astropy.io import fits
import os
import matplotlib.pyplot as plt
from astropy.stats import SigmaClip
from photutils.background import Background2D,SExtractorBackground
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from scipy import ndimage
from scipy.ndimage import gaussian_filter
plt.gray()
galaxy_name = 'DDO69'
filter = "V"
#entered by hand ------------------------------------------------------------------------------------------------------------------------------------------------
center_v = [430 , 500]
center_b = [430 , 504]
center_u = [380 , 530]
center = center_v
box_size = 200
window_size = (40, 40)
B_exp = 2400
U_exp = 1800
V_exp = 1200
exp = V_exp

# import files----------------------------------------------------------------------------------------------------------------------------------------------------
light_file = fits.open(r"C:/Users\AYSAN\Desktop/project/Galaxy\Data\DDO69\d69_V.fits")
light = light_file[0].data

starless_file = fits.open(r"C:/Users\AYSAN\Desktop/project/Galaxy\starless fits\background first\starless_DDO69_V_background_subtracted.fit")
starless = starless_file[0].data

# create background------------------------------------------------------------------------------------------------------------------------------------------------
sigma_clip = SigmaClip(sigma=3.0)
bkg_estimator = SExtractorBackground()
bkg = Background2D(light, (250 , 300) , filter_size=(3, 3),
                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
data = light - bkg.background
min_value = int(np.min(data))
Acceptable_i = []
for i in range(0,-min_value):
    newdata = data + i
    num_negative_values = np.sum(newdata < 0)
    ratio = num_negative_values / newdata.size
    if ratio < 0.005:
        Acceptable_i.append(i)
    else:
        i = i+1
i = np.min(Acceptable_i)

#background correction--------------------------------------------------------------------------------------------------------------------------------------------

corrected_light = light - bkg.background + i
print(i)

plt.imshow(bkg.background, origin='lower', cmap='Greys_r',
           interpolation='nearest')
cbar = plt.colorbar()
cbar.set_label('pixel value')
plt.title("DDO69 background (V) , boxsize = 250x300")
plt.show()
norm = ImageNormalize(vmin=0., stretch=LogStretch())
image1 = light
image2 = corrected_light
norm = ImageNormalize(vmin=0., stretch=LogStretch())
fig, axs = plt.subplots(1, 2, figsize=(10, 5))

# Display the images
im1 = axs[0].imshow(image1, origin = "lower" , aspect='auto' , norm = norm)
im2 = axs[1].imshow(image2, origin = "lower" , aspect='auto' , norm = norm)
axs[0].set_title('Light')
axs[1].set_title('Background subtracted')
fig.suptitle('DDO69 V-flter (scale = log)')
# Remove the space between the two images
plt.subplots_adjust(wspace=0.08)
# Create an axis for the colorbar on the right side of axs[1].
divider = make_axes_locatable(axs[1])
cax = divider.append_axes("right", size="5%", pad=0.1)
# Create a colorbar
cbar = fig.colorbar(im1, cax=cax)
cbar.set_label('log(pixel value)')
# Show the plot
plt.show()

#export background corrected:

output_filename = 'DDO69_V_background_subtracted.fits'
# Create a PrimaryHDU (header/data unit) from your array
primary_hdu = fits.PrimaryHDU(corrected_light)
# Create an HDUList and append the PrimaryHDU
hdul = fits.HDUList([primary_hdu])
# Write the HDUList to the FITS file
hdul.writeto(output_filename, overwrite=True)

#starless----------------------------------------------------------------------------------------------------------------------------------------------------------

image3 = starless
image2 = corrected_light
fig, axs = plt.subplots(1, 2, figsize=(10, 5))
# Display the images
im3 = axs[1].imshow(image3, origin = "lower" , aspect='auto' , norm = norm)
im2 = axs[0].imshow(image2, origin = "lower" , aspect='auto' , norm = norm)
axs[0].set_title('Background subtracted')
axs[1].set_title('starless image')
fig.suptitle('DDO69 V-filter (scale = log)')
# Remove the space between the two images
plt.subplots_adjust(wspace=0.08)
# Create an axis for the colorbar on the right side of axs[1].
divider = make_axes_locatable(axs[1])
cax = divider.append_axes("right", size="5%", pad=0.1)
# Create a colorbar
cbar = fig.colorbar(im1, cax=cax)
cbar.set_label('log(pixel value)')
# Show the plot
plt.show()


#center of mass ---------------------------------------------------------------------------------------------------------------------------------------------------
norm = ImageNormalize(vmin=0., stretch=LogStretch())
image_center_of_mass = ndimage.center_of_mass(starless)
#coordinates
x, y = image_center_of_mass[1], image_center_of_mass[0]
# Create a figure and axes
fig, ax = plt.subplots()
# Display  image
ax.imshow(starless,norm=norm,origin="lower")
# Mark the point with a red circle
circle = plt.Circle((x, y), radius=5, fill=False, color='red')
ax.add_patch(circle)
plt.show() 

pixel_scale = 1.134 #(arcsec)
starless[starless <= 0] = 1
#starless---------------------------------------------------------------------------------------------------------------------------------------------------------
norm = ImageNormalize(vmin=0., stretch=LogStretch())
plt.imshow(starless, norm = norm)
plt.title("starless %s %s"%(galaxy_name,filter))
plt.show()
#starless magnitude table------------------------------------------------------------------------------------------------------------------------------------------
flux = (starless/(exp*((pixel_scale)**2)))
magnitude_table = -2.5 * np.log10(flux) + 25
plt.imshow(magnitude_table, origin = "lower")
plt.title("magnitude table %s %s"%(galaxy_name,filter))
plt.colorbar()
plt.show()

# Slice the array
galaxy_box = starless[center[1]-box_size : center[1]+box_size , center[0]-box_size : center[0]+box_size]
plt.imshow(galaxy_box, origin = "lower")
plt.title("galaxy box %s %s"%(galaxy_name,filter))
plt.show()

#center of mass ---------------------------------------------------------------------------------------------------------------------------------------------------
norm = ImageNormalize(vmin=0., stretch=LogStretch())
image_center_of_mass = ndimage.center_of_mass(galaxy_box)
#coordinates
x, y = image_center_of_mass[1], image_center_of_mass[0]
# Create a figure and axes
fig, ax = plt.subplots()
# Display  image
ax.imshow(galaxy_box,norm=norm,origin="lower")
plt.title("center of mass %s %s"%(galaxy_name,filter))
# Mark the point with a red circle
circle = plt.Circle((x, y), radius=3, fill=False, color='red')
ax.add_patch(circle)
plt.show() 

#smoothing (moving average)------------------------------------------------------------------------------------------------------------------------------------------
flux = (galaxy_box/(exp*((pixel_scale)**2)))
mag_table = -2.5 * np.log10(flux) + 25
moving_averages = gaussian_filter(mag_table, sigma = 5, mode='reflect')
print(moving_averages)
plt.imshow(moving_averages, origin = "lower")
plt.title("smoothed galaxy box %s %s"%(galaxy_name,filter))
plt.colorbar()
plt.show()

plt.imshow(mag_table, origin = "lower")
plt.title("magnitude table %s %s"%(galaxy_name,filter))
plt.colorbar()
plt.show()

plt.figure()
plt.imshow(moving_averages, origin='lower', cmap='gray', alpha=1)
# Draw contour lines at the levels specified
CS = plt.contour(moving_averages, levels=[26,26.5,27,27.5,28,28.5,29,29.5,30,30.5,31,31.5,32])
plt.clabel(CS, inline=True, fontsize=7)
plt.title("contour lines %s %s"%(galaxy_name,filter))
plt.colorbar(CS)
# Show the plot
plt.show()







