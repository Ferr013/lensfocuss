#ifndef _fastell_h
#define _fastell_h

/*
 *
 *
 *
 * C forward function declarations.
 *
 *
 *
 *
 */

//This routine calculates the deflection due to an elliptical 
//    mass distribution, quickly and accurately.
//The parameters are position (x1in,x2in), overall factor
//     (q), power (gam) which should be between -1 and 2, axis ratio 
//     (arat) which is <=1, core radius squared (s2), and the output 
//     two-component deflection (defl).
//     The projected mass density distribution, in units of the 
//     critical density, is kappa(x1,x2)=q [u2+s2]^(-gam), where
//     u2=[x1^2+x2^2/(arat^2)].
void fastellmag_(double *x1in, double *x2in, double *q, double *gam, double *arat, double *s2, double *defl, double *magmx);
void fastelldefl_(double *x1in, double *x2in, double *q, double *gam,  double *arat, double *s2, double *defl);
void ellipphi_(double *x1in, double *x2in, double *q, double *gam, double *arat, double *s2, double *phi);


#endif //_fastell_h
