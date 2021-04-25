#include "udf.h"
#include<stdio.h>
#include<stdlib.h>
#define PI 3.141592654


DEFINE_PROPERTY(Carreau_Yasuda_model,c,t)
{
/* mu in Pa-s */
real n = 0.3568; /* power law index */
real mu_cy;
real mu_o = 0.056; /*zero shear*/
real mu_inf = 0.0035; /* inf shear*/
real lambda = 3.313; /* time constant */
real A = 2;
mu_cy = mu_inf+((mu_o-mu_inf)/pow(  (1+pow(lambda*C_STRAIN_RATE_MAG(c,t),A)),(1-n)/A));
return mu_cy;
}

DEFINE_PROPERTY(Casson_model,c,t)
{
/*Hematocrit = 0.4*/
real yield_stress = .005;
real strain = C_STRAIN_RATE_MAG(c,t);
real n_casson =0.0035;
real mu_casson;

	if (strain == 0.)
		mu_casson = n_casson;
	else
	{
		mu_casson = yield_stress/strain + pow(n_casson*yield_stress,.5)/pow(strain,0.5)+n_casson;
	}

return mu_casson;

}

DEFINE_PROPERTY(Cross_model,c,t)
{
/*Hematocrit = 0.4*/
real mu_0_cr = 0.056;
real mu_inf_cr = 0.00345;
real lambda_cr = 1.007;
real mu_cross;
real m_cr = 1.028;
mu_cross = mu_inf_cr+((mu_0_cr-mu_inf_cr)/(1+pow(lambda_cr*C_STRAIN_RATE_MAG(c,t),m_cr)));
return mu_cross;
}

DEFINE_PROPERTY(Power_model,c,t)
{
/*Hematocrit = 0.4*/
real k_pow = 0.017;
real n_pow = 0.708;
real mu_pow;
real strain = C_STRAIN_RATE_MAG(c,t);

	if (strain <= 0.05)
		mu_pow = 0.04077084195415557;
	else
	{
		mu_pow=k_pow*pow(strain,n_pow-1);
	}

return mu_pow;
}


DEFINE_PROFILE(inlet_velocity,th,i)
{
	face_t f;
	begin_f_loop(f,th)
		double t = (CURRENT_TIME*2-floor(CURRENT_TIME*2))/2; /*t is the local time within each period */

	{
		if(t <= 0.218)
			F_PROFILE(f,th,i) = 0.5*sin(4*PI*(t+0.0160236));
		else
			F_PROFILE(f,th,i) = 0.1;
	}
	end_f_loop(f,th);
}