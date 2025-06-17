#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: gregz
@author: mayad
"""

import argparse as ap
import glob
import numpy as np
import os.path as op
import sys
import warnings
from astropy.convolution import convolve, Gaussian1DKernel
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.stats import mad_std
from astropy.table import Table
from astropy.time import Time
from distutils.dir_util import mkpath
from utils import setup_logging
from astropy.stats import biweight_location as biweight
from scipy.interpolate import interp1d
from scipy.ndimage import percentile_filter
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Turn off annoying warnings (even though some deserve attention)
warnings.filterwarnings("ignore")

parser = ap.ArgumentParser(add_help=True)

parser.add_argument("folder", help='''Input folder''', type=str)

parser.add_argument('date', type=str,
                     help='''string name for date, ex: 20250613''')

parser.add_argument('outfolder', type=str,
                     help='''name of the output file''')

parser.add_argument('-n', '--name', type=str,
                    help='''Name of the science target''', default=None)

parser.add_argument("-ra", "--reduce_all",
                    help='''Reduce all files in folder''',
                    action="count", default=0)

parser.add_argument("-bl", "--bias_label",
                    help='''The objet name for bias files''',
                   type=str, default='bias')

parser.add_argument("-al", "--arc_label",
                    help='''The objet name for arc files''',
                   type=str, default='arc')

parser.add_argument("-dl", "--dark_label",
                    help='''The objet name for arc files''',
                   type=str, default='dark')

parser.add_argument("-fl", "--flat_label",
                    help='''The objet name for dome flat files''',
                   type=str, default='flat')

parser.add_argument("-tfl", "--twilight_flat_label",
                    help='''The objet name for twilight flat files''',
                   type=str, default='twi')

argv = None
args = parser.parse_args(args=argv)
folder = args.folder
outfolder = args.outfolder
# Make output folder if it doesn't exist
mkpath(outfolder)

log = setup_logging('virusp_reductions')

sns.set_context('talk')
sns.set_style('ticks')

plt.rcParams["font.family"] = "Times New Roman"

def get_script_path():
    '''
    Get script path, aka, where does Antigen live?
    '''
    return op.dirname(op.realpath(sys.argv[0]))

def identify_sky_pixels(sky, per=50, size=50):
    """
    Identifies sky pixels by applying a percentile filter and sigma-clipping.

    Parameters:
        sky (array-like): Input sky intensity values.
        per (int, optional): Percentile value for the filter. Default is 50 (median).
        size (int, optional): Size of the filter window. Default is 50.

    Returns:
        tuple: A boolean mask array indicating sky pixels and the filtered continuum array.
    """
    # Apply a percentile filter to smooth the sky data and estimate the continuum
    cont = percentile_filter(sky, per, size=size)

    try:
        # Apply sigma-clipping to identify outliers (sky pixels)
        # Use MAD-based standard deviation for robust statistics
        mask = sigma_clip(sky - cont, masked=True, maxiters=None,
                          stdfunc=mad_std, sigma=5)
    except:
        # Fallback for older versions of sigma_clip where maxiters was iters
        mask = sigma_clip(sky - cont, iters=None, stdfunc=mad_std, sigma=5)
    
    # Return the mask (True for sky pixels) and the filtered continuum
    return mask.mask, cont


def get_skymask(sky, per=50, size=50, niter=3):
    """
    Iteratively identifies and masks sky pixels in an input array by applying
    a percentile filter and sigma-clipping.

    Parameters:
        sky (array-like): Input sky intensity values.
        per (int, optional): Percentile value for the filter. Default is 50 (median).
        size (int, optional): Size of the filter window. Default is 50.
        niter (int, optional): Number of iterations for refining the sky mask. Default is 3.

    Returns:
        tuple: A boolean mask array indicating sky pixels and the final filtered continuum array.
    """
    # Keep a copy of the original sky data for final comparison
    sky_orig = sky * 1.0

    # Iteratively refine the sky mask
    for i in np.arange(niter):
        # Apply a percentile filter to estimate the continuum
        cont = percentile_filter(sky, per, size=size)

        try:
            # Apply sigma-clipping to identify sky pixels using robust statistics
            mask = sigma_clip(sky - cont, masked=True, maxiters=None,
                              stdfunc=mad_std, sigma=5, sigma_lower=500)
        except:
            # Fallback for older versions of sigma_clip
            mask = sigma_clip(sky - cont, iters=None, stdfunc=mad_std,
                              sigma=5, sigma_lower=500)

        # Update the sky values for masked pixels using the continuum
        sky[mask.mask] = cont[mask.mask]

    # Perform a final sigma-clipping pass using the original sky data
    try:
        mask = sigma_clip(sky_orig - cont, masked=True, maxiters=None,
                          stdfunc=mad_std, sigma=5, sigma_lower=500)
    except:
        mask = sigma_clip(sky_orig - cont, iters=None, stdfunc=mad_std,
                          sigma=5, sigma_lower=500)

    # Return the final mask and continuum array
    return mask.mask, cont


def rectify(scispectra, errspectra, wave_all, def_wave):
    """
    Rectifies scientific and error spectra by interpolating them onto a common wavelength grid.

    Parameters:
        scispectra (2D array): Array of scientific spectra to be rectified.
        errspectra (2D array): Corresponding error spectra for each scientific spectrum.
        wave_all (2D array): Wavelength grids corresponding to each input spectrum.
        def_wave (1D array): Target wavelength grid for interpolation.

    Returns:
        tuple:
            - scirect (2D array): Rectified scientific spectra on the target wavelength grid.
            - errorrect (2D array): Rectified error spectra on the target wavelength grid.
    """
    # Initialize arrays to store rectified scientific spectra and errors
    scirect = np.zeros((scispectra.shape[0], len(def_wave)))
    errorrect = np.zeros((scispectra.shape[0], len(def_wave)))

    # Placeholder for wavelength grid indices (used in optional refinement)
    indices1 = np.ones(scirect.shape, dtype=int)

    # Loop through each spectrum to interpolate onto the target wavelength grid
    for i in np.arange(scispectra.shape[0]):
        # Compute wavelength bin sizes for flux normalization
        dw = np.diff(wave_all[i])
        dw = np.hstack([dw[0], dw])  # Ensure length matches the wavelength array

        # Interpolate the scientific spectrum, normalizing by wavelength bin size
        scirect[i] = np.interp(def_wave, wave_all[i], scispectra[i] / dw,
                               left=np.nan, right=np.nan)

        # Interpolate the error spectrum, normalizing by wavelength bin size
        errorrect[i] = np.interp(def_wave, wave_all[i], errspectra[i] / dw,
                                 left=np.nan, right=np.nan)

        # Store indices for possible further refinement (optional block below)
        indices1[i] = np.searchsorted(wave_all[i], def_wave) + i * 1032

    # Optional: Weighted error interpolation for more accurate error propagation
    # Uncomment the following block if needed for error refinement:
    #
    # x_var = (def_wave[np.newaxis, :] * np.ones((scirect.shape[0], 1))).ravel()
    # x_fix = wave_all.ravel()
    # indices1 = indices1.ravel()
    # indices2 = indices1 - 1
    # indices2[indices2 < 0] = 0
    # indices1[indices1 >= len(x_fix)] = len(x_fix) - 1
    #
    # distances1 = np.abs(x_fix[indices1] - x_var)
    # distances2 = np.abs(x_fix[indices2] - x_var)
    # total_distance = distances1 + distances2
    # weight1 = distances1 / total_distance
    # weight2 = distances2 / total_distance
    # errorrect = (weight2**2 * errspectra.ravel()[indices1]**2 +
    #              weight1**2 * errspectra.ravel()[indices2]**2)
    # errorrect = np.sqrt(errorrect)
    # errorrect = np.reshape(errorrect, scirect.shape)
    # errorrect[np.isnan(scirect)] = np.nan

    # Return the rectified scientific and error spectra
    return scirect, errorrect


def get_fiber_to_fiber(spectrum, n_chunks=100):
    """
    Computes the fiber-to-fiber correction by normalizing each fiber's spectrum
    to the average spectrum across all fibers, then smooths the correction factor
    using interpolation.

    Parameters:
        spectrum (2D array): Array of spectra from multiple fibers (fibers x wavelength).
        n_chunks (int, optional): Number of chunks to split the wavelength range into 
                                  for smoothing. Default is 100.

    Returns:
        tuple:
            - initial_ftf (2D array): Initial fiber-to-fiber correction factors.
            - ftf (2D array): Smoothed fiber-to-fiber correction factors.
    """
    # Compute the average spectrum across all fibers using a robust biweight statistic
    average = biweight(spectrum, axis=0, ignore_nan=True)

    # Calculate the initial fiber-to-fiber correction by dividing each fiber by the average spectrum
    initial_ftf = spectrum / average[np.newaxis, :]

    # Create a wavelength grid and divide it into chunks for smoothing
    X = np.arange(spectrum.shape[1])
    x = np.array([np.mean(chunk) for chunk in np.array_split(X, n_chunks)])

    # Initialize the smoothed correction array
    ftf = spectrum * 0.

    # Loop through each fiber to compute the smoothed correction factor
    for i in np.arange(len(spectrum)):
        # Compute the biweight statistic for each chunk of the initial correction factor
        y = np.array([biweight(chunk, ignore_nan=True) for chunk in np.array_split(initial_ftf[i], n_chunks)])
        
        # Select valid (finite) values for interpolation
        sel = np.isfinite(y)
        if sel.sum() == 0.:
            continue
        # Interpolate the correction factor using quadratic interpolation
        I = interp1d(x[sel], y[sel], kind='quadratic', bounds_error=False, fill_value='extrapolate')

        # Apply the interpolation to the full wavelength range
        ftf[i] = I(X)

    # Return both the initial and smoothed fiber-to-fiber correction factors
    return initial_ftf, ftf
    

def get_wavelength(spectrum, trace, good, xref, lines, use_kernel=True, limit=100):
    """
    Computes the wavelength solution for each fiber in a spectrograph based on trace and spectral data.

    Args:
        spectrum (ndarray): 2D array of spectra, each row corresponding to a fiber.
        trace (ndarray): 2D array with trace positions for each fiber.
        good (ndarray): Boolean array indicating which fibers have valid data.
        use_kernel (bool): Whether to apply kernel smoothing when identifying arc lines. Default is True.
        limit (float): Limit on how far to search for matching arc lines. Default is 100.

    Returns:
        tuple: (wavelength, res, X, W)
            - wavelength (ndarray): Wavelength solution for each fiber.
            - res (ndarray): Residuals from the biweight mean calculation.
            - X (ndarray): Adjusted positions in trace space for arc lines.
            - W (ndarray): Arc line positions for each fiber.
    """
    
    # Initialize wavelength array and starting fiber position
    init_fiber = fiberref
    wavelength = np.zeros_like(spectrum)
    loc = xref.copy()

    # W will store arc line positions for each fiber
    W = np.zeros((trace.shape[0], len(lines)))
    mask, cont = identify_sky_pixels(spectrum[init_fiber], per=5)  # Identify sky lines
    y = spectrum[init_fiber] - cont  # Subtract continuum
    W[init_fiber] = get_arclines_fiber(y, loc, limit=limit, use_kernel=use_kernel)
    
    # Process fibers before the reference fiber in reverse order
    for i in np.arange(init_fiber)[::-1]:
        mask, cont = identify_sky_pixels(spectrum[i], per=5)  # Identify sky lines
        y = spectrum[i] - cont  # Subtract continuum
        if good[i]:  # Only process if the fiber is marked as good
            loc = get_arclines_fiber(y, loc, limit=limit,
                                     use_kernel=use_kernel)
            W[i] = loc

    # Reset location and process fibers after the reference fiber
    loc = xref.copy()
    for i in np.arange(init_fiber + 1, spectrum.shape[0]):
        mask, cont = identify_sky_pixels(spectrum[i])
        y = spectrum[i] - cont
        if good[i]:
            loc = get_arclines_fiber(y, loc, limit=limit, 
                                     use_kernel=use_kernel)
            W[i] = loc
    
    # Initialize X (adjusted trace positions) and residuals array
    X = np.zeros_like(W)
    xall = np.arange(trace.shape[1])
    res = np.zeros(W.shape[1])
    
    # Interpolate missing values and fit polynomial to each arc line
    for i in range(W.shape[1]):
        x = np.zeros(W.shape[0])
        bad = np.where(~good)[0]
        gind = np.where(good)[0]

        # Fill in missing arc line positions for bad fibers
        for b in bad:
            W[b, i] = W[gind[np.argmin(np.abs(b - gind))], i]

        # Interpolate positions in trace space
        for j in range(W.shape[0]):
            x[j] = np.interp(W[j, i], xall, trace[j])

        # Fit a 4th-order polynomial to arc line positions
        sel = (W[:, i] > 0) * np.isfinite(x)
        
        X[:, i] = np.polyval(np.polyfit(x[sel], W[sel, i], 4), x)

        # Compute residuals using biweight mean
        res[i] = mad_std(X[:, i] - W[:, i], ignore_nan=True)
    fits.HDUList([fits.PrimaryHDU(X), fits.ImageHDU(lines)]).writeto('test.fits', overwrite=True)
    goodlines = res < 0.2
    # Compute final wavelength solution for each fiber
    for j in range(W.shape[0]):
        if good[j]:
            wavelength[j] = np.polyval(np.polyfit(X[j][goodlines], 
                                                  lines[goodlines], 5), xall)
    
        
    # Plot wavelength solution for inspection
    plot_wavelength(lines, W, wavelength)
    
    return wavelength, res, X, W


def get_arclines_fiber(spectrum, init_loc=None, limit=100, use_kernel=True):
    """
    Identifies arc line positions in a given spectrum by detecting peaks.

    Args:
        spectrum (ndarray): 1D array representing the spectrum of a fiber.
        init_loc (ndarray, optional): Initial guess locations for arc lines. Default is None.
        limit (float): Minimum peak value to consider a valid arc line. Default is 1000.
        use_kernel (bool): Whether to apply a box kernel convolution to smooth the spectrum. Default is True.

    Returns:
        ndarray: Array of arc line positions (pixel indices) in the spectrum.
    """

    # Apply box kernel convolution to smooth the spectrum if use_kernel is True
    if use_kernel:
        B = Gaussian1DKernel(1.0)
        y1 = convolve(spectrum, B)
    else:
        y1 = spectrum.copy()

    # Identify peaks in the spectrum by finding zero-crossings in the first derivative
    diff_array = y1[1:] - y1[:-1]
    loc = np.where((diff_array[:-1] > 0) & (diff_array[1:] < 0))[0]

    # Filter peaks based on the limit threshold
    peaks = y1[loc + 1]
    loc = loc[peaks > limit] + 1
    peaks = y1[loc]

    # Helper function to refine peak positions using quadratic interpolation
    def get_trace_chunk(flat, XN):
        YM = np.arange(flat.shape[0])
        inds = np.zeros((3, len(XN)))
        inds[0] = XN - 1
        inds[1] = XN
        inds[2] = XN + 1
        inds = inds.astype(int)

        # Quadratic interpolation to refine peak positions
        Trace = (YM[inds[1]] - (flat[inds[2]] - flat[inds[0]]) /
                 (2. * (flat[inds[2]] - 2. * flat[inds[1]] + flat[inds[0]])))
        return Trace

    # Refine peak positions using quadratic interpolation
    loc = get_trace_chunk(y1, loc)

    # Match refined peak positions with initial guess locations, if provided
    if init_loc is not None:
        final_loc = []
        for i in init_loc:
            final_loc.append(loc[np.argmin(np.abs(np.array(loc) - i))])
        loc = final_loc

    return loc


def get_spectra(array_flt, array_trace, npix=5):
    """
    Extract spectra by dividing the flat field and averaging the central pixels.

    Parameters
    ----------
    array_flt : 2D numpy array
        Twilight image.
    array_trace : 2D numpy array
        Trace for each fiber.
    npix : int, optional
        Number of pixels to extract around the trace center. Default is 5.

    Returns
    -------
    spec : 2D numpy array
        Extracted and rectified spectrum for each fiber.
    """

    # Initialize the output spectrum array
    spec = np.zeros((array_trace.shape[0], array_trace.shape[1]))

    # Get the number of rows in the flat field image
    N = array_flt.shape[0]

    # Create an array of x-axis pixel indices
    x = np.arange(array_flt.shape[1])

    # Calculate the lower and upper bounds for pixel extraction
    LB = int((npix + 1) / 2)  # Lower bound
    HB = -LB + npix + 1       # Upper bound

    # Iterate through each fiber
    for fiber in np.arange(array_trace.shape[0]):

        # Skip fibers with trace positions too close to the image edges
        if np.round(array_trace[fiber]).min() < LB:
            continue
        if np.round(array_trace[fiber]).max() >= (N - LB):
            continue

        # Convert trace positions to integer indices
        indv = np.round(array_trace[fiber]).astype(int)

        # Iterate through pixels around the trace center
        for j in np.arange(-LB, HB):

            # Calculate weight for the lower boundary pixel
            if j == -LB:
                w = indv + j + 1 - (array_trace[fiber] - npix / 2.)

            # Calculate weight for the upper boundary pixel
            elif j == HB - 1:
                w = (npix / 2. + array_trace[fiber]) - (indv + j)

            # Assign weight 1 for central pixels
            else:
                w = 1.

            # Add the weighted pixel values to the spectrum
            spec[fiber] += array_flt[indv + j, x] * w

    # Normalize the spectrum by the number of extracted pixels
    return spec / npix


def get_spectra_error(array_flt, array_trace, npix=5):
    '''
    Extract spectra by dividing the flat field and averaging the central
    two pixels
    
    Parameters
    ----------
    array_flt : 2d numpy array
        Twilight image
    array_trace : 2d numpy array
        Trace for each fiber
    npix : int, optional
        Number of pixels for averaging (default is 5)
    
    Returns
    -------
    twi_spectrum : 2d numpy array
        Rectified twilight spectrum for each fiber  
    '''
    
    # Initialize spectrum array to store extracted spectra
    spec = np.zeros((array_trace.shape[0], array_trace.shape[1]))

    # Get number of rows in the flat field image
    N = array_flt.shape[0]

    # Create an array of x-coordinates for the flat field image
    x = np.arange(array_flt.shape[1])

    # Calculate bounds for pixel averaging
    LB = int((npix + 1) / 2)
    HB = -LB + npix + 1

    # Iterate over each fiber to extract its spectrum
    for fiber in np.arange(array_trace.shape[0]):
        # Skip fibers with traces too close to image edges
        if np.round(array_trace[fiber]).min() < LB:
            continue
        if np.round(array_trace[fiber]).max() >= (N - LB):
            continue

        # Convert trace positions to integer indices
        indv = np.round(array_trace[fiber]).astype(int)

        # Loop over neighboring pixels for averaging
        for j in np.arange(-LB, HB):
            if j == -LB:
                # Calculate weight for lower boundary pixels
                w = indv + j + 1 - (array_trace[fiber] - npix / 2.)
            elif j == HB - 1:
                # Calculate weight for upper boundary pixels
                w = (npix / 2. + array_trace[fiber]) - (indv + j)
            else:
                # Set weight to 1 for central pixels
                w = 1.

            # Accumulate weighted sum of squared values from the flat field
            spec[fiber] += array_flt[indv + j, x] ** 2 * w

    # Return the root mean square error normalized by npix
    return np.sqrt(spec) / npix


def get_spectra_chi2(array_flt, array_sci, array_err, array_trace, npix=5):
    '''
    Extract spectra by dividing the flat field and averaging the central
    two pixels
    
    Parameters
    ----------
    array_flt : 2d numpy array
        Twilight image
    array_sci : 2d numpy array
        Science image
    array_err : 2d numpy array
        Error estimate for each pixel
    array_trace : 2d numpy array
        Trace for each fiber
    npix : int, optional
        Number of pixels for averaging (default is 5)
    
    Returns
    -------
    spec : 2d numpy array
        Chi-squared spectra for each fiber  
    '''

    # Initialize spectrum array to hold chi-squared values
    spec = np.zeros((array_trace.shape[0], array_trace.shape[1]))

    # Get the number of rows in the flat field image
    N = array_flt.shape[0]

    # Create an array of x-coordinates for the images
    x = np.arange(array_flt.shape[1])

    # Calculate bounds for pixel averaging
    LB = int((npix + 1) / 2)
    HB = -LB + npix + 1

    # Iterate over each fiber to extract its chi-squared spectrum
    for fiber in np.arange(array_trace.shape[0]):
        # Initialize a chi-squared array with shape (npix+1, 3, len(x))
        chi2 = np.zeros((npix + 1, 3, len(x)))

        # Skip fibers with traces too close to the image edges
        if np.round(array_trace[fiber]).min() < LB:
            continue
        if np.round(array_trace[fiber]).max() >= (N - LB):
            continue

        # Convert trace positions to integer indices
        indv = np.round(array_trace[fiber]).astype(int)

        # Loop over neighboring pixels for averaging
        for j in np.arange(-LB, HB):
            # Calculate weights for boundary pixels
            if j == -LB:
                w = indv + j + 1 - (array_trace[fiber] - npix / 2.)
            elif j == HB - 1:
                w = (npix / 2. + array_trace[fiber]) - (indv + j)
            else:
                # Use a weight of 1 for central pixels
                w = 1.

            # Apply weights to science, flat field, and error images
            chi2[j + LB, 0] = array_sci[indv + j, x] * w
            chi2[j + LB, 1] = array_flt[indv + j, x] * w
            chi2[j + LB, 2] = array_err[indv + j, x] * w

        # Compute the normalization factor for the flux
        norm = chi2[:, 0].sum(axis=0) / chi2[:, 1].sum(axis=0)

        # Calculate the chi-squared numerator: (data - model)^2
        num = (chi2[:, 0] - chi2[:, 1] * norm[np.newaxis, :]) ** 2

        # Calculate the denominator: (error + regularization term)^2
        denom = (chi2[:, 2] + 0.01 * chi2[:, 0].sum(axis=0)[np.newaxis, :]) ** 2

        # Compute the chi-squared value for each fiber
        spec[fiber] = 1. / (1. + npix) * np.sum(num / denom, axis=0)

    # Return the final chi-squared spectrum array
    return spec


def get_trace(twilight, ref):
    """
    Extract fiber traces from a twilight flat field image.

    Parameters
    ----------
    twilight : 2d numpy array
        Twilight flat field image used to determine fiber locations.

    Returns
    -------
    trace : 2d numpy array
        The calculated trace positions for each fiber across the image.
    good : 1d numpy array (boolean)
        Boolean mask indicating which fibers are valid (non-missing).
    """

    # Determine the number of valid (good) fibers
    N1 = np.isfinite(ref['px']).sum()
    good = np.where(np.isfinite(ref['px']))[0]  # Indices of good fibers

    # Helper function to calculate trace positions for a chunk of the image
    def get_trace_chunk(flat, XN):
        # YM represents the y-axis pixel coordinates
        YM = np.arange(flat.shape[0])

        # Create a 3-row array for XN-1, XN, and XN+1 indices
        inds = np.zeros((3, len(XN)))
        inds[0] = XN - 1.
        inds[1] = XN + 0.
        inds[2] = XN + 1.
        inds = np.array(inds, dtype=int)

        # Calculate the trace using a second-order derivative method
        Trace = (YM[inds[1]] - (flat[inds[2]] - flat[inds[0]]) /
                 (2. * (flat[inds[2]] - 2. * flat[inds[1]] + flat[inds[0]])))
        return Trace

    # Assign the input image to a variable
    image = twilight

    # Determine the number of chunks based on whether the image is binned
    N = 80

    # Split the x-axis into chunks and calculate the mean x-position for each chunk
    xchunks = np.array([np.mean(x) 
                        for x in np.array_split(np.arange(image.shape[1]), N)])

    # Split the image into vertical chunks
    chunks = np.array_split(image, N, axis=1)

    # Calculate the mean flat field for each chunk
    flats = [np.mean(chunk, axis=1) for chunk in chunks]

    # Initialize an array to hold the trace positions for each fiber
    Trace = np.zeros((len(ref), len(chunks)))

    # Initialize a counter and a list to store peak positions
    k = 0
    P = []

    # Iterate over each chunk to calculate the fiber traces
    for flat, x in zip(flats, xchunks):
        # Calculate the difference between adjacent pixels
        diff_array = flat[1:] - flat[:-1]

        # Identify peaks by finding zero-crossings in the difference array
        loc = np.where((diff_array[:-1] > 0.) & (diff_array[1:] < 0.))[0]
        loc = loc[loc > 2]  # Ignore peaks near the image edges

        # Filter out weak peaks
        peaks = flat[loc + 1]
        loc = loc[peaks > 0.3 * np.median(peaks)] + 1

        # Store the detected peak positions
        P.append(loc)

        # Get the trace positions for the detected peaks
        trace = get_trace_chunk(flat, loc)

        # Initialize an array to hold the trace for this chunk
        T = np.zeros((len(ref)))

        # If the number of detected peaks exceeds the number of good fibers, trim the excess
        if len(trace) > N1:
            trace = trace[-N1:]

        # If the number of detected peaks matches the number of good fibers
        if len(trace) == N1:
            T[good] = trace
            # Interpolate missing fibers based on nearby good fibers
            for missing in np.where(np.isnan(ref['px']))[0]:
                gind = np.argmin(np.abs(missing - good))
                T[missing] = T[good[gind]] + ref['px'][missing] - ref['px'][good[gind]]

        # If the number of detected peaks matches the total number of fibers
        if len(trace) == len(ref):
            T = trace

        # Store the calculated trace for this chunk
        Trace[:, k] = T
        k += 1

    # Fit a 7th-order polynomial to smooth the traces across the x-axis
    x = np.arange(twilight.shape[1])
    trace = np.zeros((Trace.shape[0], twilight.shape[1]))
    for i in np.arange(Trace.shape[0]):
        sel = Trace[i, :] != 0.
        trace[i] = np.polyval(np.polyfit(xchunks[sel], Trace[i, sel], 7), x)

    # Return the final trace array and the good fiber mask
    return trace, np.isfinite(ref['px']), Trace, xchunks


def plot_wavelength(lines, W, wavelength):
    '''
    Plots the residuals of the wavelength solution using a violin plot.
    
    Parameters
    ----------
    lines : 1D ndarray
        Expected wavelengths of arc lamp lines.
    W : 2D ndarray
        Measured positions of the arc lines for different fibers.
    wavelength : 2D ndarray
        Wavelength solution for the fibers.
    
    Returns
    -------
    None.
    '''
    
    # Prepare data for seaborn violin plot
    data = []
    num_fibers = W.shape[0]
    for i, line in enumerate(lines):
        residuals = []
        for fiber in range(2, num_fibers): # skip the first two broken fibers
            X = np.arange(wavelength.shape[1])
            pred = np.interp(W[fiber, i], X, wavelength[fiber])
            if np.all(wavelength[fiber] != 0.):
                residuals.append(line - pred)
        for res in residuals:
            data.append({'Wavelength': line, 'Residual': res})
    
    df = pd.DataFrame(data)
    
    # Create the violin plot
    plt.figure(figsize=(14, 7))
    plt.gca().set_position([0.15, 0.19, 0.75, 0.71])
    sns.violinplot(x='Wavelength', y='Residual', data=df, palette='coolwarm',
                   inner=None, saturation=0.8)
    # Customize plot appearance
    plt.xlabel('Arc Line Wavelength ($\mathrm{\AA}$)')
    plt.ylabel(r'Measured - Expected ($\mathrm{\AA}$)')
    plt.title('Wavelength Solution Residuals')
    plt.xticks(rotation=45)

    # Save the plot as a PNG file with the given name
    plt.savefig(op.join(outfolder, 'wavelength_measures.png'))

def plot_trace(full_trace, trace, x, orders=[5, 130, 230]):

    '''
    Plots the residuals of the trace correction and saves the figure.
    
    Parameters
    ----------
    trace_cor : 2D ndarray
        The array of trace correction residuals to be plotted.
    name : str
        A name or identifier used for saving the plot.

    Returns
    -------
    None.
    
    '''

    X = np.arange(full_trace.shape[1])
    # Create a figure with specified size
    plt.figure(figsize=(8, 7))
    colors = plt.get_cmap('Set2')(np.linspace(0, 1, len(orders)))
    for order, color in zip(orders, colors):
        mean_trace = np.mean(full_trace[order])
        plt.scatter(x, trace[order] - mean_trace, color='k', edgecolor='k',
                    s=30,)
        plt.scatter(x, trace[order] - mean_trace, color=color, edgecolor='k',
                    s=20, alpha=0.5)
        plt.plot(X, full_trace[order] - mean_trace, color=color, lw=1,
                 label='Fiber: %i' % (order+1))

    plt.legend()

    # Adjust the appearance of the ticks on both axes
    ax = plt.gca()
    ax.tick_params(axis='both', which='both', direction='in', zorder=3)
    ax.tick_params(axis='y', which='both', left=True, right=True)
    ax.tick_params(axis='x', which='both', bottom=True, top=True)
    ax.tick_params(axis='both', which='major', length=8, width=2)
    ax.tick_params(axis='both', which='minor', length=4, width=1)
    ax.minorticks_on() 

    # Label the axes
    plt.xlabel('Column')
    plt.ylabel('Trace - Mean(Trace)')

    # Save the plot as a PNG file with the given name
    plt.savefig(op.join(outfolder, 'trace_measures.png'))
    

def prep_image(image, channel):
    '''
    Orient the images from blue to red (left to right)
    Fibers are oriented to match configuration files
    
    Parameters
    ----------
    image : 2d numpy array
        fits image
    amp : str
        Amplifier for the fits image
    ampname : str
        Amplifier name is the location of the amplifier
    
    Returns
    -------
    image : 2d numpy array
        Oriented fits image correcting for what amplifier it comes from
        These flips are unique to the VIRUS/LRS2 amplifiers
    '''
    
    overscan_length = 32
    bias_value = biweight(image[:, -(overscan_length-2):])
    image = image[:, :-overscan_length] - bias_value

    if channel == "b":
        image[:] = image[::-1, :]
    if channel == "g":
        image[:] = image[:, ::-1]
    if channel == "r":
        image[:] = image[::-1, :]
    if channel == "d":
        image[:] = image[:, ::-1]
    # Overscan subtraction

    
    return image


def base_reduction(data, masterbias, channel):
    """
    Perform basic image reduction by applying bias subtraction, gain correction, 
    and calculating the error estimate.

    Parameters
    ----------
    data : 2d numpy array
        Raw input image to be reduced.
    masterbias : 2d numpy array
        Master bias frame to be subtracted from the image.

    Returns
    -------
    image : 2d numpy array
        Reduced image with bias subtracted and gain applied.
    E : 2d numpy array
        Error estimate for each pixel, including read noise and photon noise.
    """

    # Preprocess the raw image (e.g., background subtraction, padding)
    image = prep_image(data, channel)

    # Subtract the master bias from the image
    image[:] -= masterbias

    # Apply gain correction to convert counts to electrons
    image[:] *= gain

    # Calculate the error estimate (read noise + photon noise)
    E = np.sqrt(rdnoise**2 + np.where(image > 0., image, 0.))

    # Return the reduced image and the error estimate
    return image, E


def subtract_sky(spectra, good):
    """
    Subtract the sky background from spectra by identifying sky fibers 
    and performing a biweight calculation.

    Parameters
    ----------
    spectra : 2d numpy array
        The input spectra data where rows represent fibers and columns represent wavelengths.
    good : 1d numpy array of bools
        Boolean mask indicating which fibers are good (non-sky).

    Returns
    -------
    2d numpy array
        Spectra with the sky background subtracted for each fiber.
    """
    
    # Get the number of fibers and number of wavelength bins
    nfibs, N = spectra.shape

    # Define range for biweight calculation (middle third of the data)
    n1 = int(1. / 3. * N)
    n2 = int(2. / 3. * N)

    # Calculate the biweight of spectra over the middle third of each fiber's data
    y = biweight(spectra[:, n1:n2], axis=1, ignore_nan=True)

    # Identify sky pixels based on the biweighted data and apply a mask
    mask, cont = identify_sky_pixels(y[good], size=15)

    # Create a mask for fibers that are not good and are sky fibers
    m1 = ~good
    m1[good] = mask
    skyfibers = ~m1

    # Compute the biweighted sky spectrum based on sky fibers
    init_sky = biweight(spectra[skyfibers], axis=0, ignore_nan=True)

    # Subtract the sky spectrum from the original spectra
    return spectra - init_sky[np.newaxis, :]


def get_pca_sky_residuals(data, ncomponents=5):
    """
    Perform PCA on the input data to extract the principal components and 
    reconstruct the data using a specified number of components.

    Parameters
    ----------
    data : 2d numpy array
        Input data where rows represent samples (e.g., spectra) and columns represent features (e.g., wavelengths).
    ncomponents : int, optional
        The number of principal components to retain for reconstruction. Default is 5.

    Returns
    -------
    pca : PCA object
        The fitted PCA model.
    A : 2d numpy array
        The reconstructed data using the first `ncomponents` principal components.
    """

    # Initialize PCA with the specified number of components
    pca = PCA(n_components=ncomponents)

    # Fit the PCA model and transform the data into the principal components space
    H = pca.fit_transform(data)

    # Reconstruct the data using the principal components
    A = np.dot(H, pca.components_)

    # Return the PCA model and the reconstructed data
    return pca, A


def get_residual_map(data, pca):
    """
    Compute the residual map by subtracting a PCA-based model from the input data.
    
    Parameters
    ----------
    data : 2d numpy array
        Input data where rows represent samples (e.g., spectra) and columns represent features.
    pca : PCA object
        Fitted PCA model used for reconstruction of the data.
    
    Returns
    -------
    res : 2d numpy array
        The residual map, which is the difference between the input data and the PCA model.
    """

    # Initialize the residual map with zeros
    res = data * 0.

    # Loop over each column (feature) in the input data
    for i in np.arange(data.shape[1]):

        # Compute the absolute deviation from the median for each column
        A = np.abs(data[:, i] - np.nanmedian(data[:, i]))

        # Identify the "good" data points (those within 3 times the median deviation)
        good = A < (3. * np.nanmedian(A))

        # Select finite values and good points for the model fitting
        sel = np.isfinite(data[:, i]) * good

        # Compute the PCA coefficients for the selected data points
        coeff = np.dot(data[sel, i], pca.components_.T[sel])

        # Reconstruct the model based on the computed coefficients
        model = np.dot(coeff, pca.components_)

        # Store the residual (difference) for this feature in the residual map
        res[:, i] = model

    # Return the computed residual map
    return res



def get_arc_pca(arcskysub, good, mask, components=15):
    """
    Perform PCA on the arc-sky-subtracted data with preprocessing to remove 
    bad data points and mask the non-relevant pixels.

    Parameters
    ----------
    arcskysub : 2d numpy array
        The arc-sky-subtracted data to be processed.
    good : 1d numpy array of bools
        Boolean mask indicating which fibers are good (non-bad).
    mask : 1d numpy array of bools
        Mask indicating the relevant pixels (e.g., fiber locations).
    components : int, optional
        Number of PCA components to retain. Default is 15.

    Returns
    -------
    pca : PCA object
        The fitted PCA model.
    """

    # Initialize X as the arc-sky-subtracted data
    X = arcskysub

    # Set values outside the mask to 0 and apply the "good" fiber mask
    X[:, ~mask] = 0.
    X[~good] = 0.

    # Transpose X to have samples in rows and features in columns
    X = X.swapaxes(0, 1)

    # Calculate the mean and standard deviation along the rows (fiber axis)
    M = np.nanmean(X, axis=1)
    Var = np.nanstd(X, axis=1)

    # Prevent division by zero by setting variance of zero to 1
    Var[Var == 0.] = 1.

    # Normalize the data by subtracting the mean and dividing by the variance
    X = (X - M[:, np.newaxis]) / Var[:, np.newaxis]

    # Replace NaN values in the normalized data with 0
    X[np.isnan(X)] = 0.

    # Perform PCA on the normalized data
    pca, A = get_pca_sky_residuals(X, ncomponents=components)

    # Return the fitted PCA model
    return pca


def get_continuum(skysub, masksky, nbins=50):
    """
    Compute the continuum by interpolating the biweighted median of 
    sky-subtracted data, with masked values excluded.

    Parameters
    ----------
    skysub : 2d numpy array
        Sky-subtracted data where each row represents a spectrum.
    masksky : 1d numpy array of bools
        Mask indicating which sky pixels should be excluded in the continuum calculation.
    nbins : int, optional
        Number of bins to divide the spectrum into for the biweight calculation. Default is 50.

    Returns
    -------
    bigcont : 2d numpy array
        The continuum for each spectrum in the sky-subtracted data.
    """

    # Initialize the output array for the continuum with zeros
    bigcont = skysub * 0.

    # Loop over each row (spectrum) in the sky-subtracted data
    for j in np.arange(skysub.shape[0]):
        # Copy the current spectrum and mask the sky pixels
        y = skysub[j] * 1.
        y[masksky] = np.nan

        # Divide the spectrum into bins and calculate the mean of each bin
        x = np.array([np.mean(chunk) for chunk in np.array_split(np.arange(len(y)), nbins)])

        # Calculate the biweighted median for each bin
        z = np.array([biweight(chunk, ignore_nan=True) for chunk in np.array_split(y, nbins)])

        # Select bins with finite values for interpolation
        sel = np.isfinite(z)

        # If there are enough valid bins, perform quadratic interpolation
        if sel.sum() > 5:
            I = interp1d(x[sel], z[sel], kind='quadratic', bounds_error=False, 
                         fill_value='extrapolate')

            # Store the interpolated continuum for the current spectrum
            bigcont[j] = I(np.arange(len(y)))

    # Return the computed continuum for all spectra
    return bigcont

    
def reduce(fn, biastime_list, masterbias_list, flttime_list,
           trace_list, wave_time, wave_list, ftf_list, channel, 
           pca=None):
    """
    Reduce the raw data by performing a series of processing steps, 
    including bias subtraction, flat-fielding, sky subtraction, 
    and PCA-based residuals analysis.

    Parameters
    ----------
    fn : str
        The filename of the FITS file containing the data.
    biastime_list : list
        List of times associated with bias frames.
    masterbias_list : list
        List of master bias frames for bias correction.
    flttime_list : list
        List of times associated with flat field frames.
    trace_list : list
        List of fiber trace data.
    wave_time : list
        List of times associated with wavelength calibration.
    wave_list : list
        List of wavelength calibration data.
    ftf_list : list
        List of flat-field corrections.
    pca : PCA object, optional
        A pre-fitted PCA model for residual map analysis. Default is None.

    Returns
    -------
    pca : PCA object
        The fitted PCA model, returned if `pca` is None.
    continuum : 1d numpy array
        The computed continuum for the spectrum.
    """

    # Open the FITS file and extract the observation time
    f = fits.open(fn)
    t = Time(f[0].header['DATE-OBS'] + 'T' + f[0].header['UT'])
    mtime = t.mjd

    # Select appropriate master bias frame based on observation time
    masterbias = masterbias_list[get_cal_index(mtime, biastime_list)]

    # Perform basic image reduction (bias subtraction, gain adjustment)
    image, E = base_reduction(f[0].data, masterbias, channel)

    # Get the fiber trace and selection mask for the current observation
    trace, good = trace_list[get_cal_index(mtime, flttime_list)]

    # Extract spectra from the image using the trace data
    spec = get_spectra(image, trace)

    # Calculate the spectrum error using the flat-field and error image
    specerr = get_spectra_error(E, trace)

    # Compute the chi-square of the spectrum to identify bad pixels
    chi2 = get_spectra_chi2(masterflt - masterbias, image, E, trace)
    badpix = chi2 > 20.  # Pixels with chi2 > 20 are considered bad
    specerr[badpix] = np.nan
    spec[badpix] = np.nan

    # Retrieve the wavelength calibration data for the current observation
    wavelength = wave_list[get_cal_index(mtime, wave_time)]

    # Rectify the spectrum and error based on the wavelength
    specrect, errrect = rectify(spec, specerr, wavelength, def_wave)

    # Apply flat-field correction
    ftf = ftf_list[get_cal_index(mtime, flttime_list)]
    specrect[:] /= (ftf * f[0].header['EXPTIME'])
    errrect[:] /= (ftf * f[0].header['EXPTIME'])

    # Generate a sky mask and the continuum for sky subtraction
    skymask, cont = get_skymask(biweight(specrect, axis=0, ignore_nan=True), size=25)

    # Subtract the sky from the spectrum
    skysubrect = subtract_sky(specrect, good)

    # If PCA is not provided, compute it from the sky-subtracted data
    if pca is None:
        pca = get_arc_pca(skysubrect, good, skymask, components=pca_comp)
        return pca

    # Adjust the sky mask and compute the continuum
    skymask[1:] += skymask[:-1]
    skymask[:-1] += skymask[1:]
    cont1 = get_continuum(skysubrect, skymask, nbins=50)

    # Compute the residuals by subtracting the continuum
    Z = skysubrect - cont1
    res = get_residual_map(Z, pca)

    # Mask out residuals where sky mask is not valid
    res[:, ~skymask] = 0.0

    # Write the final reduced data to a new FITS file
    write_fits(skysubrect - res, skysubrect, specrect, errrect, f[0].header)

    # Return the biweighted spectrum and continuum
    return biweight(specrect, axis=0,ignore_nan=True), cont


def write_fits(skysubrect_adv, skysubrect, specrect, errorrect, header):
    """
    Writes the sky-subtracted, rectified spectra and error data to a FITS file, 
    preserving the header information and adding necessary meta-information.

    Parameters
    ----------
    skysubrect_adv : 2D numpy array
        The advanced sky-subtracted spectrum.
    skysubrect : 2D numpy array
        The basic sky-subtracted spectrum.
    specrect : 2D numpy array
        The rectified spectrum.
    errorrect : 2D numpy array
        The error associated with the rectified spectrum.
    header : FITS header
        The header information to be preserved in the output FITS file.
    """
    
    hdulist = []  # List to store HDU objects for the FITS file

    # Loop through the data arrays and create HDUs for each
    for image, ftp in zip([skysubrect_adv, skysubrect, specrect, errorrect], 
                          [fits.PrimaryHDU, fits.ImageHDU, fits.ImageHDU, fits.ImageHDU]):
        
        # Create an HDU object from each image, setting it to 'float32' type
        hdu = ftp(np.array(image, dtype='float32'))
        
        # Remove any conflicting CD matrix elements first
        for key in ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'CDELT1', 'CDELT2']:
            if key in hdu.header:
                del hdu.header[key]
        
        # Define your wavelength solution
        wavelength_step = def_wave[1] - def_wave[0]  # Compute wavelength step
        
        # Set WCS parameters correctly
        hdu.header['CRVAL1'] = def_wave[0]  # First wavelength (Angstroms)
        hdu.header['CRPIX1'] = 1  # First pixel corresponds to first wavelength
        hdu.header['CD1_1'] = wavelength_step  # Set CD1_1 to match wavelength step
        hdu.header['CTYPE1'] = 'WAVE'  # Spectral axis label
        
        # Set fiber axis metadata
        hdu.header['CRVAL2'] = 1  # Reference value for fiber index
        hdu.header['CRPIX2'] = 1  # First pixel for fiber axis
        hdu.header['CD2_2'] = 1  # Step of 1 fiber per index
        hdu.header['CTYPE2'] = 'FIBER'  # Labeling fiber axis

        # Copy relevant keys from the input header, avoiding duplicates
        for key in header.keys():
            if key in hdu.header:
                continue
            if ('CCDSEC' in key) or ('DATASEC' in key):  # Exclude CCDSEC and DATASEC keys
                continue
            if ('BSCALE' in key) or ('BZERO' in key):  # Exclude BSCALE and BZERO keys
                continue
            try:
                hdu.header[key] = header[key]  # Copy header data to the new HDU
            except:
                continue

        # Format the output filename using observation date and object name
        t = Time(header['DATE-OBS'] + 'T' + header['UT'])
        objname = '_'.join(header['OBJECT'].split())
        iname = '_'.join([objname, t.strftime('%Y%m%dT%H%M%S'), 'multi'])  # Generate filename

        # Append the HDU to the list
        hdulist.append(hdu)

    # Write the HDU list to the output file, overwriting if necessary
    fits.HDUList(hdulist).writeto(op.join(outfolder, iname + '.fits'), overwrite=True)


def make_mastercal_list(filenames, breakind, channel):
    """
    Creates a list of master calibration images and corresponding times 
    by splitting the input list of filenames at given indices.

    Parameters
    ----------
    filenames : list of str
        List of FITS file paths containing calibration data.
    breakind : list of int
        List of indices to split the filenames into different chunks.

    Returns
    -------
    masters : list of 2D numpy arrays
        List of median calibration images for each chunk.
    times : list of float
        List of mean observation times (MJD) corresponding to each chunk.
    """

    # Define break points for splitting the filenames into chunks
    breakind1 = np.hstack([0, breakind])  # Start indices for chunks
    breakind2 = np.hstack([breakind, len(filenames)+1])  # End indices for chunks

    masters = []  # List to store median calibration images
    times = []    # List to store mean observation times

    # Iterate over the chunks defined by breakind1 and breakind2
    for bk1, bk2 in zip(breakind1, breakind2):
        # Collect and preprocess frames within the current chunk
        frames = [prep_image(fits.open(f)[0].data, channel) 
                  for cnt, f in enumerate(filenames) 
                  if ((cnt > bk1) * (cnt < bk2))]  # Only include frames in the current chunk
        
        # Extract observation times (MJD) for frames in the current chunk
        t = [Time(fits.open(f)[0].header['DATE-OBS']).mjd 
             for cnt, f in enumerate(filenames) 
             if ((cnt > bk1) * (cnt < bk2))]
        

        # Append the median frame and the mean time for the current chunk
        masters.append(np.nanmedian(frames, axis=0))
        times.append(np.mean(t))

    return masters, times


def get_cal_index(mtime, time_list):
    """
    Finds the index of the closest calibration time to a given observation time.

    Parameters
    ----------
    mtime : float
        The observation time (MJD) to compare against the calibration times.
    time_list : list of float
        A list of calibration times (MJD).

    Returns
    -------
    int
        The index of the closest calibration time in `time_list`.
    """

    # Find the index of the closest time in time_list to mtime by minimizing the absolute time difference
    return np.argmin(np.abs(mtime - np.array(time_list)))


def get_filenames(gnames, typelist, names):
    """
    Finds filenames that match a list of keywords by checking if any of the keywords 
    appear in the associated types.

    Parameters
    ----------
    gnames : list of str
        List of filenames to check.
    typelist : list of str
        List of types or descriptions corresponding to the filenames.
    names : list of str
        List of keywords to search for within the types/descriptions.

    Returns
    -------
    np.ndarray
        Array of filenames from `gnames` that match any of the keywords in `names`.
    """

    matches = []  # List to store matched filenames
    # Iterate through each filename and its corresponding type
    for gn, tp in zip(gnames, typelist):
        # Check if any keyword appears in the type (case-insensitive)
        for name in names:
            if name.lower() in str(tp).lower():
                matches.append(gn)  # Append matching filename to the list
    return np.array(matches)  # Return matched filenames as a numpy array


# =============================================================================
# Get Folder and Filenames
# =============================================================================
ROOT_DATA_PATH = args.folder
date = args.date
allfilenames = sorted(glob.glob(op.join(ROOT_DATA_PATH, 'VIRUS2', date, 
                                     '*', '*', '*.fits')))
unit_list = [fn.split('_')[-4] for fn in allfilenames]
units = np.unique(unit_list)
observations = np.unique([fn.split('_')[-6]  for fn in allfilenames])
DIRNAME = get_script_path()

fiber_radius = 2.483 / 2.
pca_comp = 15


# =============================================================================
# Get Line List
# =============================================================================
for unit in units:
    channel = unit[-1].lower()
    if channel == 'b':
        def_wave = np.linspace(3700., 4630., 2064)
        continue
    if channel == 'g':
        def_wave = np.linspace(4610., 5925., 2064)
        line_list = Table.read(op.join(DIRNAME, 'line_list', 
                                        'virus2_green.txt'), format='ascii')
        limit = 20
        gain = 2.017
        rdnoise = 3.09
    if channel == 'r':
        def_wave = np.linspace(5900., 7610., 2064)
        continue
    if channel == 'd':
        def_wave = np.linspace(7590., 9300., 2064)
        continue

    filenames = [fn for fn, un in zip(allfilenames, unit_list) if un == unit]

    fiberref = 130
    lines = np.array(line_list['wavelength'])
    xref = np.array(line_list['col'])
    use_kernel = True
    
    
    # =============================================================================
    # Make a list of the objects for each filename, ignore those without 'OBJECT' 
    # in the header
    # =============================================================================
    typelist = []
    gnames = []
    timelist = []
    for fn in filenames:
        try:
            obj = fn.split('_')[-5] # fits.open(f)[0].header['OBJECT']
            dateobs = fn.split('_')[-2].split('T')[0]
            dateobs = '%s-%s-%s' % (dateobs[:4], dateobs[4:6], dateobs[6:])
            datestring = dateobs # fits.open(f)[0].header['DATE-OBS']
        except:
            continue
        typelist.append(obj)
        gnames.append(fn)
        timelist.append(Time(datestring))
    
    # =============================================================================
    # Get the bias filenames, domeflat filenames, and arc lamp filenames
    # =============================================================================
    log.info('Sorting Files')
    bnames = [args.bias_label]
    dnames = [args.dark_label]
    anames = [args.arc_label]
    tnames = [args.twilight_flat_label]
    dfnames = [args.flat_label]
    snames = ['feige', 'bd']
    bias_filenames = get_filenames(gnames, typelist, bnames)
    twiflt_filenames = get_filenames(gnames, typelist, tnames)
    domeflt_filenames = get_filenames(gnames, typelist, dfnames)
    arc_filenames = get_filenames(gnames, typelist, anames)
    std_filenames = get_filenames(gnames, typelist, snames)
    if args.reduce_all:
        gna = []
        for gn in gnames:
            if gn in bias_filenames:
                continue
            if gn in arc_filenames:
                continue
            if gn in twiflt_filenames:
                continue
            if gn in domeflt_filenames:
                continue
            gna.append(gn)
        sci_filenames = np.array(gna)
    else:
        if args.name is not None:
            sci_filenames = get_filenames(gnames, typelist, [args.name])
        else:
            sci_filenames = []
    
    if len(twiflt_filenames) > 0:
        flt_filenames = twiflt_filenames
    else:
        flt_filenames = domeflt_filenames
    
    # =============================================================================
    # Use the file numbers for connecting blocks of observations
    # =============================================================================
    biasnum = [int(op.basename(op.dirname(op.dirname(fn)))) for fn in bias_filenames]
    fltnum = [int(op.basename(op.dirname(op.dirname(fn)))) for fn in flt_filenames]
    arcnum = [int(op.basename(op.dirname(op.dirname(fn)))) for fn in arc_filenames]
    bias_breakind = np.where(np.diff(biasnum) > 1)[0]
    flt_breakind = np.where(np.diff(fltnum) > 1)[0]
    arc_breakind = np.where(np.diff(arcnum) > 1)[0]
    
    # Load reference fiber locations from a predefined file
    ref = Table.read(op.join(DIRNAME, 'IFUcen', 'IFUcen_VIRUS2_D3G.txt'), format='ascii')  
    ref.reverse()
    # =============================================================================
    # Make a master bias, master dome flat, and master arc for the first set of OBS
    # =============================================================================
    log.info('Making master bias frames')
    masterbias_list, biastime_list = make_mastercal_list(bias_filenames,
                                                         bias_breakind, channel)
    
    log.info('Making master flat frames')
    
    masterflt_list, flttime_list = make_mastercal_list(flt_filenames,
                                                       flt_breakind, channel)
    
    log.info('Making master arc frames')
    masterarc_list, arctime_list = make_mastercal_list(arc_filenames,
                                                       arc_breakind, channel)
    
    
    # =============================================================================
    # Get trace from the dome flat
    # =============================================================================
    trace_list = []
    fltspec = []
    log.info('Getting trace for each master flat')
    for masterflt, mtime in zip(masterflt_list, flttime_list):
        masterbias = masterbias_list[get_cal_index(mtime, biastime_list)]
        trace, good, Tchunk, xchunk = get_trace(masterflt-masterbias, ref)
        plot_trace(trace, Tchunk, xchunk)
        trace_list.append([trace, good])
        domeflat_spec = get_spectra(masterflt-masterbias, trace)
        domeflat_error = 0. * domeflat_spec
        fltspec.append([domeflat_spec, domeflat_error])
    
    # =============================================================================
    # Get wavelength from arc lamps
    # =============================================================================
    wave_list = []
    wave_time = []
    bk1 = np.hstack([0, arc_breakind+1])
    log.info('Getting wavelength for each master arc')
    import traceback
    for masterarc, mtime, bk in zip(masterarc_list, arctime_list, bk1):
        masterbias = masterbias_list[get_cal_index(mtime, biastime_list)]
        trace, good = trace_list[get_cal_index(mtime, flttime_list)]
        lamp_spec = get_spectra(masterarc-masterbias, trace)
        fits.PrimaryHDU(lamp_spec).writeto('test.fits',overwrite=True)
        try:
            wavelength, res, X, W = get_wavelength(lamp_spec, trace, good, 
                                                   xref, lines, limit=limit, 
                                                   use_kernel=use_kernel)
        except:
            log.warning('Could not get wavelength solution for masterarc')
            log.warning('First file of failed masterarc included: %s' %
                        (arc_filenames[bk]))
            traceback.print_exc()
            continue
        wave_list.append(wavelength)
        wave_time.append(mtime)
    
    # =============================================================================
    # Rectify domeflat spectra and get fiber to fiber
    # =============================================================================
    ftf_list = []
    log.info('Getting fiber to fiber for each master domeFlat')
    for fltsp, mtime in zip(fltspec, flttime_list):
        wavelength = wave_list[get_cal_index(mtime, wave_time)]        
        domeflat_spec, domeflat_error = fltsp
        domeflat_rect, domeflat_error_rect = rectify(domeflat_spec, domeflat_error,
                                                     wavelength, def_wave)
        ftf, ftf_smooth = get_fiber_to_fiber(domeflat_rect)
        ftf_list.append(ftf)
        
    
    pca = reduce(arc_filenames[0], biastime_list, masterbias_list, flttime_list,
                 trace_list, wave_time, wave_list, ftf_list, channel, pca=None)
    for fn in sci_filenames:
        log.info('Reducing: %s' % fn)
        sky, cont = reduce(fn, biastime_list, masterbias_list, flttime_list,
                           trace_list, wave_time, wave_list, ftf_list, 
                           channel, pca=pca)
