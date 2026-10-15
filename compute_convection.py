import numpy as np
import matplotlib.pyplot as plt


def measure_convective_plume(Vr_patch, patch_physical_width):
    """
    Measures the characteristic size of convective plumes using 2D spatial autocorrelation.

    Parameters:
    Vr_patch (2D numpy array): A square patch of radial velocities extracted from the envelope.
    patch_physical_width (float): The physical width of this patch in Solar Radii (R_sun).

    Returns:
    plume_diameter (float): The measured characteristic size of the plumes in R_sun.
    target_cell_size (float): The ideal Arepo cell diameter (D / 10) in R_sun.
    """

    # 1. Subtract the mean (we only want the velocity fluctuations)
    # Replace any NaNs with 0 so the FFT doesn't crash
    dV = np.nan_to_num(Vr_patch)
    dV = dV - np.mean(dV)

    # 2. Compute 2D Autocorrelation via FFT (Wiener-Khinchin theorem)
    F = np.fft.fft2(dV)
    power_spectrum = np.abs(F) ** 2
    autocorr_2d = np.fft.ifft2(power_spectrum).real

    # Shift the zero-frequency (perfect correlation lag=0) to the center
    autocorr_2d = np.fft.fftshift(autocorr_2d)

    # Normalize so the center peak is exactly 1.0
    autocorr_2d /= np.max(autocorr_2d)

    # 3. Radially Average the 2D Autocorrelation
    # Create grids of X and Y pixel indices
    y, x = np.indices(autocorr_2d.shape)
    center_y, center_x = np.array(autocorr_2d.shape) // 2

    # Calculate the integer distance (radius in pixels) of each pixel from the center
    r_pixels = np.sqrt((x - center_x) ** 2 + (y - center_y) ** 2).astype(int)

    # Fast way to average values in concentric rings using bincount
    ring_sums = np.bincount(r_pixels.ravel(), autocorr_2d.ravel())
    ring_counts = np.bincount(r_pixels.ravel())
    radial_profile = ring_sums / ring_counts

    # 4. Find the first Zero-Crossing
    # The point where correlation drops below 0 is the plume radius
    zero_crossings = np.where(np.diff(np.sign(radial_profile)))[0]

    if len(zero_crossings) == 0:
        print("Error: No zero-crossing found. Try a larger patch of the envelope.")
        return None, None

    plume_radius_pixels = zero_crossings[0]

    # 5. Convert back to physical units
    pixel_size_rsun = patch_physical_width / Vr_patch.shape[0]
    plume_radius_rsun = plume_radius_pixels * pixel_size_rsun
    plume_diameter_rsun = plume_radius_rsun * 2

    # The "Goldilocks" Arepo resolution (10 cells across the plume)
    target_cell_size = plume_diameter_rsun / 10.0

    # --- Optional: Plot the result to verify visually ---
    plt.figure(figsize=(8, 4))
    r_physical = np.arange(len(radial_profile)) * pixel_size_rsun
    plt.plot(r_physical, radial_profile, lw=2)
    plt.axhline(0, color='black', linestyle='--')
    plt.axvline(plume_radius_rsun, color='red', linestyle='--', label=f'Radius = {plume_radius_rsun:.1f} R_sun')
    plt.title("Spatial Autocorrelation of Convection")
    plt.xlabel("Distance ($R_\odot$)")
    plt.ylabel("Correlation")
    plt.legend()
    plt.xlim(0, plume_diameter_rsun * 1.5)
    plt.show()

    print(f"Measured Plume Diameter: {plume_diameter_rsun:.2f} R_sun")
    print(f"Ideal Arepo Cell Diameter (D/10): {target_cell_size:.2f} R_sun")

    return plume_diameter_rsun, target_cell_size