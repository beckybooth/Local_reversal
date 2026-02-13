import numpy as np
import healpy as hp
from scipy.integrate import quad

#########################################################################################################

rGC = 8.15 # distance to Galactic centre

#########################################################################################################

def long_lat(x, y, z):
    '''
    args: 
    float x: heliocentric x
    float y: heliocentric y
    float z: heliocentric z

    returns:
    float: heliocentric r
    float: Galactic longitude in degrees
    float: Galactic latitude in degress
    '''

    r = np.sqrt(x*x + y*y + z*z)
    
    long = np.degrees(np.arctan2(y,x))
    lat = np.degrees(np.arctan2(z,r))
    
    return long, lat

def xyz(r, l, b):
    '''
    args:
    float r: heliocentric r
    float l: Galactic longitude in radians
    float b: Galactic latitude in radians
    
    returns:
    float: heliocentric x
    float: heliocentric y
    float: heliocentric z
    '''

    x = r*np.cos(l)*np.cos(b)
    y = r*np.sin(l)*np.cos(b)
    z = r*np.sin(b)

    return x,y,z

########################################################################################################
# Reversal geometry functions

def get_ABC(elln, bn):
    '''
    args:
    float elln: radians
    float bn: radians

    returns:
    float: A, B, C
    '''

    A = np.cos(bn)*np.cos(elln)
    B = np.cos(bn)*np.sin(elln)
    C = np.sin(bn)

    return A, B, C

def reversal_intersect(l, b, elln, bn, x0):
    '''
    args:
    float l: Galactic longitude in radians
    float b: Galactic latitude in radians
    float elln: reversal plane horizontal angle in radians
    float bn: reversal plane vertical angle in radians
    float x0: distance to plane along x-axis in kpc

    returns:
    float: distance (kpc) to intersection with along line of sight (l,b)
    '''

    A, B, C = get_ABC(elln, bn)

    x = np.cos(l)*np.cos(b)
    y = np.sin(l)*np.cos(b)
    z = np.sin(b)
    
    return (A*x0) / ( A*x + B*y + C*z )

#########################################################################################################

# simulated Faraday depth (no approximation)

def B_line_of_sight(l, b, a, alpha, p, cosb, sinb, cosbeta, sinbeta):
    '''
    args:
    float l: Galactic longitude in radians
    float b: Galactic latitude in radians
    int a: +1 for a CW field and -1 for a CCW field
    float alpha: Galactocentric azimuthal angle in radians
    float p: magnetic field pitch in radians
    cosb: the cosine of Galactic latitude (precalculated to make integration more efficient)
    sinb: the sine of Galactic latitude (precalculated to make integration more efficient)
    cosbeta: the cosine of magnetic field vertical tilt (precalculated to make integration more efficient)
    sinbeta: the sine of magnetic field vertical tilt (precalculated to make integration more efficient)
    
    returns:
    float: line of sight component of magnetic field unit vector (assumes rhat is directed away)
    '''

    return -a*cosb*cosbeta*np.sin(l + p - alpha) + sinb*sinbeta

def integrand(r, a, l, b, p, cosl, sinl, cosb, sinb, cosbeta, sinbeta, neB):
    '''
    args:
    float r: heliocentric r
    int a: +1 for a CW field and -1 for a CCW field
    float l: Galactic longitude in radians
    float b: Galactic latitude in radians
    float p: magnetic field pitch in radians
    cosl: the cosine of Galactic longitude (precalculated to make integration more efficient)
    sinl: the sine of Galactic longitude (precalculated to make integration more efficient)
    cosb: the cosine of Galactic latitude (precalculated to make integration more efficient)
    sinb: the sine of Galactic latitude (precalculated to make integration more efficient)
    cosbeta: the cosine of magnetic field vertical tilt (precalculated to make integration more efficient)
    sinbeta: the sine of magnetic field vertical tilt (precalculated to make integration more efficient)

    returns:
    float: The integrand for the Faraday depth numerical integration 0.812*neB*Blos
    '''

    # calculate Galactocentric azimuthal angle
    xgc = r*cosb*cosl - rGC
    ygc = r*cosb*sinl
    zgc = r*sinb
    rgc = np.sqrt(xgc*xgc + ygc*ygc + zgc*zgc)
    alpha = np.arctan2(ygc,xgc)

    # get  line of sight component of magnetic field unit vector
    BLOS = B_line_of_sight(l, b, a, alpha, p, cosb, sinb, cosbeta, sinbeta)

    return  BLOS*neB

def phi_sim_full(l, b,  R,  betaCW, betaCCW, x0, pitch = 11.5, elln = 168.5, bn = -60, neB = 0.04):
    '''
    args:
    float l: Galactic longitude in degrees
    float b: Galactic latitude in degrees
    float R: path length in kpc
    float betaCW: tilt angle of the CW field (unreversed field) in degrees
    float betaCCw: tilt angle of the CCW field (reversed field) in degrees
    float x0: distance to plane along x-axis in kpc
    float p: magnetic field pitch in degrees (default = 11.5)
    float elln: reversal plane horizontal angle in degrees (default = 168.5)
    float bn: reversal plane vertical angle in degrees (default = -60)

    returns:
    float: the simulated Faraday depth out to path length R using scipy.integrate.quad for numerical integration
    '''

    # convert to radians
    elln = np.radians(elln)
    bn = np.radians(bn)
    l = np.radians(l)
    b = np.radians(b)
    pitch = np.radians(pitch)

    # precalculate to make integration more efficent
    sinb = np.sin(b)
    cosb = np.cos(b)
    cosl = np.cos(l)
    sinl = np.sin(l)


    # distance to the reversal along the current LOS
    Rp = reversal_intersect(l, b, elln, bn, x0)
    
    # do integration
    if Rp > R or Rp < 0:
        # if the LOS does NOT intersect the reversal

        sinbeta = np.sin(np.radians(betaCW))
        cosbeta = np.cos(np.radians(betaCW))
        a = 1 # CW

        return quad(integrand, R, 0 , args = (a, l, b, pitch, cosl, sinl, cosb, sinb, cosbeta, sinbeta, neB))[0] * 0.812*1000 
        # (multiply by 1000 to convert kpc to pc)

    else:
        
        # if the LOS DOES intersect the reversal
   
        sinbeta = np.sin(np.radians(betaCW))
        cosbeta = np.cos(np.radians(betaCW))
        a = 1 # CW

        front = quad(integrand, Rp, 0 , args = (a, l, b, pitch, cosl, sinl, cosb, sinb, cosbeta, sinbeta, neB))[0] * 0.812*1000

        sinbeta = np.sin(np.radians(betaCCW))
        cosbeta = np.cos(np.radians(betaCCW))
        a = -1 # CCW

        back = quad(integrand, R, Rp , args = (a, l, b, pitch,  cosl, sinl, cosb, sinb, cosbeta, sinbeta, neB))[0] * 0.812*1000
        
        #print(a, front, back)

        return (front + back) 
    

#########################################################################################################

# simulated Faraday depth approximations, assuming alpha = 180 degrees (runs faster)


def M1screen(l, b,  R,  betaCW, betaCCW, x0, neB, pitch = 11.5, elln = 168.5, bn = -60):
    '''
    args:
    float l: Galactic longitude in degrees
    float b: Galactic latitude in degrees
    float R: path length in kpc
    float betaCW: tilt angle of the CW field (unreversed field) in degrees
    float betaCCw: tilt angle of the CCW field (reversed field) in degrees
    float x0: distance to plane along x-axis in kpc
    float neB: electron density times magnetic field strength free parameter
    float p: magnetic field pitch in degrees (default = 11.5)
    float elln: reversal plane horizontal angle in degrees (default = 168.5)
    float bn: reversal plane vertical angle in degrees (default = -60)

    returns:
    float: the simulated M1 for a screen (assuming alpha = 180 degrees and uniform ne|B|)
    '''
     
    # convert to radians
    elln = np.radians(elln)
    bn = np.radians(bn)
    l = np.radians(l)
    b = np.radians(b)
    pitch = np.radians(pitch)
    betaCW = np.radians(betaCW)
    betaCCW = np.radians(betaCCW)

    # distance to the reversal along the current LOS
    Rp = reversal_intersect(l, b, elln, bn, x0)

    epsilon = -0.812*neB
    zeta = np.cos(b)*np.cos(betaCW)*np.sin(l+pitch) + np.sin(b)*np.sin(betaCW)
    eta = -np.cos(b)*np.cos(betaCCW)*np.sin(l+pitch) + np.sin(b)*np.sin(betaCCW)

    if Rp > R or Rp < 0:
    # if the LOS does NOT intersect the reversal

        return epsilon * zeta * R * 1000
        # (multiply by 1000 to convert kpc to pc)
    else:
    # if the LOS DOES intersect the reversal    

        return epsilon*( zeta*Rp + eta*(R - Rp)) * 1000 

def M1slab(l, b,  R,  betaCW, betaCCW, x0, neB, pitch = 11.5, elln = 168.5, bn = -60):
    '''
    args:
    float l: Galactic longitude in degrees
    float b: Galactic latitude in degrees
    float R: path length in kpc
    float betaCW: tilt angle of the CW field (unreversed field) in degrees
    float betaCCw: tilt angle of the CCW field (reversed field) in degrees
    float x0: distance to plane along x-axis in kpc
    float neB: electron density times magnetic field strength free parameter
    float p: magnetic field pitch in degrees (default = 11.5)
    float elln: reversal plane horizontal angle in degrees (default = 168.5)
    float bn: reversal plane vertical angle in degrees (default = -60)

    returns:
    float: the simulated M1 for a slab (assuming alpha = 180 degrees and uniform ne|B|)
    '''

    # convert to radians
    elln = np.radians(elln)
    bn = np.radians(bn)
    l = np.radians(l)
    b = np.radians(b)
    pitch = np.radians(pitch)
    betaCW = np.radians(betaCW)
    betaCCW = np.radians(betaCCW)

    # distance to the reversal along the current LOS
    Rp = reversal_intersect(l, b, elln, bn, x0)

    epsilon = -0.812*neB
    zeta = np.cos(b)*np.cos(betaCW)*np.sin(l+pitch) + np.sin(b)*np.sin(betaCW)
    eta = -np.cos(b)*np.cos(betaCCW)*np.sin(l+pitch) + np.sin(b)*np.sin(betaCCW)

    if Rp > R or Rp < 0:
    # if the LOS does NOT intersect the reversal

        return 0.5 *epsilon * zeta * R *1000
        # (multiply by 1000 to convert kpc to pc)

    else:
    # if the LOS DOES intersect the reversal  

        term1 = zeta * (Rp*R - 0.5*Rp*Rp)
        term2 = eta * (0.5*R*R + 0.5*Rp*Rp - Rp*R)
        
        return (epsilon * (term1 + term2))/R *1000
    
#########################################################################################################
# extra functions for visualizing the reversal plane
    

def get_reversal_plane(x0, elln = 168.5, bn = -60, dmax = 1, numpoints = 100):
    '''
    Returns the position the reversal plane in Cartesian coordinates (x, y, z) kpc

    args:
    float x0: distance to plane along x-axis in kpc
    float elln: reversal plane horizontal angle in degrees (default = 168.5)
    float bn: reversal plane vertical angle in degrees (default = -60)
    float dmax: the maximum x and y distance along the plane (default = 1 kpc)
    int numpoints: the number of data points along each axis (default = 100)

    returns:
    array, array, array
    the x, y, and z coordinates of the reversal plane
    '''


# Returns the position the reversal plane in Cartesian coordinates (x, y, z) kpc

    elln = np.radians(elln)
    bn = np.radians(bn)

    A, B, C = get_ABC(elln, bn)

    x = np.linspace(-dmax, dmax, numpoints)
    y = np.linspace(-dmax, dmax, numpoints)
    xx, yy = np.meshgrid(x, y)
    x = xx.flatten()
    y = yy.flatten()

    z = (-A*(x-x0) - B*y)/C

    return x, y, z


def get_reversal_plane_long_lat(x0, elln = 168.5, bn = -60, dmax = 1, numpoints = 100):
    '''
    Returns the position the reversal plane in Galactic coordinates (r, l, b)

    args:
    float x0: distance to plane along x-axis in kpc
    float elln: reversal plane horizontal angle in degrees (default = 168.5)
    float bn: reversal plane vertical angle in degrees (default = -60)
    float dmax: the maximum x and y distance along the plane (default = 1 kpc)
    int numpoints: the number of data points along each axis (default = 100)

    returns:
    array, array, array
    the r(kpc), longitude (deg), latitude (deg), coordinates of the reversal plane
    '''

    x, y, z = get_reversal_plane(x0, elln, bn, dmax, numpoints)
    r = np.sqrt(x**2 + y**2 + z**2)
    long, lat = long_lat(x, y, z)

    return r, long, lat