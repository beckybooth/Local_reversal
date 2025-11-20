# A three-dimensional model for the reversal in the local large-scale interstellar magnetic field

### The reversal plane

The reversal plane is positioned between the Sun and Sagittarius Arm and tilted towards the Sun above the Galactic midplane and away below. 
Above the reversal, the large-scale magnetic field is counterclockwise (red arrows) as viewed from the North Galactic Pole. Below the reversal, the magnetic field is clockwise (blue arrows).
The positive x-axis points towards the Galactic centre. 

<img src = "https://github.com/beckybooth/Local_reversal/blob/main/plane_geometry1.png">

### Parameters tha describe the reversal plane geometry are: 

$\ell_n$ - the horizontal tilt angle of the plane, the horizontal angle between the positive x-axis and the plane normal

$b_n$ - the vertical tilt angle of the plane, the vertical angle between the positive x-axis and the plane normal

$x_0$ = the distance between the SUn and the point on the positive x-axis that intersects the plane

<img src = "https://github.com/beckybooth/Local_reversal/blob/main/plane_geometry2.png">

### Functions in local_reversal.py

<hr>

#### reversal_intersect(l, b, elln, bn, x0)
  
Returns distance (kpc) to the reversal plane along line of sight (l,b)

$$R_{p} = \frac{A x_0}{A \cos b \cos\ell + B \cos b \sin \ell + C \sin b}$$

where

$A   = \cos b_n \ \cos \ell_n$

$B   = \cos b_n \ \sin \ell_n$

$C   = \sin b_n$

<hr>

#### B_line_of_sight(l, b, a, alpha, p, cosb, sinb, cosbeta, sinbeta)

Returns the line of sight component of magnetic field unit vector (assumes rhat is directed away)

Clockwise field: a = +1

Counterclockwise field a = -1

$$B_\parallel = -a \ \cos b \ \cos \beta \ \sin(\ell + p - \alpha) + \sin b \ \sin\beta$$

<hr>

### Simulated Faraday depth (no approximation)

#### phi_sim_full(l, b,  R,  betaCW, betaCCW, x0, pitch = 11.5, elln = 168.5, bn = -60)

Returns the simulated Faraday depth for LOS $(\ell, b)$ out to path length R using scipy.integrate.quad for numerical integration

For now we assume $n_e\|B| = 0.04\\ cm^{-3}\\mu\textrm{G}$

For a LOS where $R_p > R\textrm{ or }R_p < 0$ (does not intersect the reversal):

$$\phi_{\textrm{no-rev}}(r) = 0.812  \int_{r}^{0} n_e\ B_{\parallel CW} \ dr^\prime$$

For a LOS where $R > R_p$ (does intersect the reversal):

$$\phi_{\textrm{rev}}(r) = 0.812 { \int_{r}^{R_p} n_e\, B_{\parallel CCW} \ dr^\prime +\int_{R_p}^0 n_e\ B_{\parallel CW} \, dr^\prime\}$$

<hr>

### Simulated Faraday depth approximations, assuming alpha = 180 degrees (runs faster)

#### M1screen(l, b,  R,  betaCW, betaCCW, x0, neB, pitch = 11.5, elln = 168.5, bn = -60)

Returns the simulated M1 for a screen (assuming alpha = 180 degrees and uniform ne|B|)

For a LOS where $R_p > R\textrm{ or }R_p < 0$ (does not intersect the reversal):

$$ M1_{\textrm{screen}_{\textrm{no-rev}}} = \epsilon \ \zeta \ R $$

For a LOS where $R > R_p$ (does intersect the reversal):

$$ M1_{\textrm{screen}_{\textrm{rev}}} =  \epsilon \ [\zeta \ R_p + \eta \ (R-R_p) \ ] $$

where

$\epsilon   = - 0.812 \ n_e |B|$

$\zeta   =\cos(b) \ \cos(\beta_{CW}) \ \sin(\ell + p) + \sin(b) \ \sin(\beta_{CW})$

$\eta   =-\cos(b) \ \cos(\beta_{CCW}) \ \sin(\ell + p) + \sin(b) \ \sin(\beta_{CCW}) $

#### M1slab(l, b,  R,  betaCW, betaCCW, x0, neB, pitch = 11.5, elln = 168.5, bn = -60)

Returns the simulated M1 for a slab (assuming alpha = 180 degrees and uniform ne|B|)

For a LOS where $R_p > R\textrm{ or }R_p < 0$ (does not intersect the reversal):

   $$ M1_{{\textrm{slab}}_{\textrm{no-rev}}}  = \frac{1}{2} \ \epsilon \ \zeta \ R $$

For a LOS where $R > R_p$ (does intersect the reversal):

$$M1_{\textrm{slab}_{\textrm{rev}}} = \frac{\epsilon \ [\zeta \ (R_pR - \frac{1}{2}R_p^2) + \eta\ (\frac{1}{2}R^2 + \frac{1}{2} \ R_p^2 - R_p R) \ ]}{R}$$

where

$\epsilon   = - 0.812 \ n_e |B|$

$\zeta   =\cos(b) \ \cos(\beta_{CW}) \ \sin(\ell + p) + \sin(b) \ \sin(\beta_{CW})$

$\eta   =-\cos(b) \ \cos(\beta_{CCW}) \ \sin(\ell + p) + \sin(b) \ \sin(\beta_{CCW}) $

<hr>

### Functions to demonstrate reversal plane geometry

#### get_reversal_plane(x0, elln = 168.5, bn = -60, dmax = 1, numpoints = 100)

returns the x, y, and z coordinates of the reversal plane (array, aray, array)

####  get_reversal_plane_long_lat(x0, elln = 168.5, bn = -60, dmax = 1, numpoints = 100)

returns the r(kpc), longitude (deg), latitude (deg), coordinates of the reversal plane (array, aray, array)
