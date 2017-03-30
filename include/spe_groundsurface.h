#ifndef __NCPA_SPE_GROUNDSURFACE_H__
#define __NCPA_SPE_GROUNDSURFACE_H__


/*******************************************************************************************/
/* Internal integral in P_G(x): the one over theta                                         */
/* params = { r_y, z_y, r_x, theta_x, phi_x, k, R, dhdy1, dhdy2, phi_1, phi_2, phi_3 }     */
/*******************************************************************************************/
double ground_surface_integral_over_theta( double theta, void *params );

/*******************************************************************************************/
/* External integral in P_G(x): the one over r                                         */
/* params = { z_y, r_x, theta_x, phi_x, k, R, dhdy1, dhdy2, phi_1, phi_2, phi_3 }     */
/*******************************************************************************************/
double ground_surface_integral_over_r( double r, void *params );







#endif