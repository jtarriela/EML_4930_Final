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
mu_cy = mu_inf+((mu_o-mu_inf)/pow((1+pow(lambda*C_STRAIN_RATE_MAG(c,t),A)),(1-n)/A));
return mu_cy;
}

DEFINE_PROPERTY(Casson_model,c,t)
{
/*Hematocrit = 0.4*/
real mu_p = 0.00145;
real Hct = .04;
real mu_casson;
real mu_inf_casson;
real N_inf;
real strain;

strain = C_STRAIN_RATE_MAG(c,t);
mu_inf_casson = pow(0.625*Hct,0.5);
N_inf = pow(mu_p*pow((1-Hct),-0.25),0.5);


	if (strain == 0.)
		mu_casson = pow(N_inf,2);
	else
	{
		mu_casson = pow(mu_inf_casson,2)/strain+2*mu_inf_casson*N_inf/pow(strain,0.5)+pow(N_inf,2);
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

mu_pow = k_pow*pow(C_STRAIN_RATE_MAG(c,t),n_pow-1);
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