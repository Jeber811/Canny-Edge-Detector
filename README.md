# Canny-Edge-Detector (marrh.c)

Multi-stage edge detection algorithm for PGM formatted image files. Involves:

Gaussian function calculation, non-maxima suppression of magnitudes perpendicular to the angle of the gradient, automatic HI and LO pixel value generation, and hysteresis thresholding calculation of the HI and LO pixel values

Takes sigma and percent values as user-inputted parameters.

# Sobel-Edge-Detector (Sobel.c)

Gradient-based edge detection algorithm for PGM formatted image files. Involves:

Computation of horizontal and vertical intensity gradients using fixed Sobel convolution masks, calculation of gradient magnitudes from the X and Y responses, normalization of all magnitude values to the 0–255 range, and binary thresholding using user-provided HI and LO pixel values.

Takes HI and LO threshold values as user-inputted parameters.
