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

plt.gray()


light_file = fits.open(r"C:/Users\AYSAN\Desktop/project/Galaxy\Data\DDO69\d69_V.fits")
light = light_file[0].data


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

starless_file = fits.open(r"C:/Users\AYSAN\Desktop/project/Galaxy\starless fits\background first\starless_DDO69_V_background_subtracted.fit")
starless = starless_file[0].data

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


output_filename = 'DDO69_V_background_subtracted.fits'

# Create a PrimaryHDU (header/data unit) from your array
primary_hdu = fits.PrimaryHDU(corrected_light)

# Create an HDUList and append the PrimaryHDU
hdul = fits.HDUList([primary_hdu])

# Write the HDUList to the FITS file
hdul.writeto(output_filename, overwrite=True)


