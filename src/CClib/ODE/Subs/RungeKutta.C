#include <cmath>
#include <ODE.h>


int ODE::RungeKuttaStep(double* y, double* yp, double* dy, double t,double h)
{
//
//	Runge-Kutta, 4th order,	for ODE	dy/dt = f(y,t)
//		calculate dy = y(t+dt) - y(t),	where y(t) = y and dt = h 
//		and return error status
//
	double *y1 = dy;		// y1 is intermediate storage only
		
	double h_half = h/2.0;
	double h_sixth = h/6.0;
	double t_center = t + h_half;
	
	double *ptr_k1, *ptr_k2, *ptr_k3, *ptr_k4;
	double *ptr_y0, *ptr_y0_p, *ptr_y1, *ptr_dy;
	
	int n;
	
	// Compute k1 = h/2*f(y, t)		and y1 = y + k1

	// yp = f(y,t)
	
	for(ptr_k1 = k1, ptr_y0_p = yp,
				ptr_y0 = y, ptr_y1 = y1, n = n_dim; n--; )
	{
		*ptr_k1 = (*ptr_y0_p++)*h_half;
		*ptr_y1++ = *ptr_y0++ + *ptr_k1++;
	}
	
	
	// Compute k2 = h/2*f(y+k1, t + h/2)	and y1 = y + k2
	
	if(F(k2, y1, t_center))
		return FAILED;
	
	for(ptr_k2 = k2, ptr_y1 = y1, ptr_y0 = y, n = n_dim; n--; )
	{
		*ptr_k2 *= h_half;
		*ptr_y1++ = *ptr_y0++ + *ptr_k2++;
	}
	
	
	// Compute k3 = h*f(y+k2, t + h/2)	and y1 = y + k3
	
	if(F(k3, y1, t_center))
		return FAILED;
	
	for(ptr_k3 = k3, ptr_y1 = y1, ptr_y0 = y, n = n_dim; n--; )
	{
		*ptr_k3 *= h;
		*ptr_y1++ = *ptr_y0++ + *ptr_k3++;
	}
	
	
	
	// Compute k4 = h*f(y+k3, t + h)
	// and dy = 1/6*(2*k1 + 4*k2 + 2*k3 + k4)
	
	if(F(k4, y1, t+h))
		return FAILED;
		
	for(ptr_k1=k1,ptr_k2=k2,ptr_k3=k3,ptr_k4=k4,
					ptr_dy=dy, n=n_dim; n--;)
	{
		*ptr_dy++ = (1.0/3.0)*(*ptr_k1++ + *ptr_k3++) +
				(2.0/3.0)*(*ptr_k2++) + h_sixth*(*ptr_k4++);
	}
	
	return 0;
}
