import numpy as np
import phd_toolbox.constants as cst

# ORBIT EFFECT
@staticmethod
def orbit_shift(I_L, r_w, phi_w, Brho, mu0=cst.mu0_H_m):
    '''
    This functions computes the orbit shift induced by a wire in both planes.
    Taken from https://indico.cern.ch/event/456856/contributions/1968793/attachments/1196177/1740660/BBLR_LYON_2015.pdf

    Args:
        I_L: float or array, current in the wire [A.m]
        r_w: float or array, beam-wire distance, regardless the plane [m]
        phi_w: float or array, angle between the beam and the wires, zero being horizontal [rad]
        Brho: float, magnetic rigidity
        mu0: vacuum permeability [H/m] 
    
    Returns:
        [dx,dy]: lst of dim 2 (for x and y). Each element can be either a float or an array depending on the input.
    '''
    dx = mu0*I_L*np.cos(phi_w)/2/np.pi/Brho/r_w
    dy = mu0*I_L*np.sin(phi_w)/2/np.pi/Brho/r_w
    return [dx,dy]

# TUNE EFFECT
@staticmethod
def tune_shift(I_L, r_w, betx, bety, phi_w, Brho, mu0=cst.mu0_H_m):
    '''
    This functions computes the tune shift induced by a wire in both planes.
    Taken from https://indico.cern.ch/event/456856/contributions/1968793/attachments/1196177/1740660/BBLR_LYON_2015.pdf

    Args:
        I_L: float or array, current in the wire [A.m]
        r_w: float or array, beam-wire distance, regardless the plane [m]
        phi_w: float or array, angle between the beam and the wires, zero being horizontal [rad]
        betx,bety: floats, beta functions at the wire location
        Brho: float, magnetic rigidity
        mu0: vacuum permeability [H/m] 
    
    Returns:
        [dqx,dqy]: lst of dim 2 (for x and y). Each element can be either a float or an array depending on the input.
    '''
    dqx = -mu0*I_L*betx*np.cos(2*phi_w)/8/np.pi/np.pi/Brho/r_w/r_w
    dqy = mu0*I_L*bety*np.cos(2*phi_w)/8/np.pi/np.pi/Brho/r_w/r_w
    return [dqx,dqy]
