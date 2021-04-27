/* This file saved as unix __.c file */
/* Used on CIRCE cluster */

#include "udf.h"
#include<stdio.h>
#include<stdlib.h>
#define PI 3.141592654


DEFINE_PROPERTY(Carreau_Yasuda,c,t)
{
real n = 0.3568;
real mu_cy;
real mu_o = 0.056;
real mu_inf = 0.0035;
real lambda = 3.313;
real A = 2;
mu_cy = mu_inf+((mu_o-mu_inf)/pow((1+pow(lambda*C_STRAIN_RATE_MAG(c,t),A)),(1-n)/A));

return mu_cy;
}

DEFINE_PROPERTY(Casson,c,t)
{
/*Hematocrit = 0.4*/
real mu_p = 0.00145;
real Hct = .04;
real mu_casson;
real mu_inf_casson;
real N_inf;


mu_inf_casson = pow(0.625*Hct,0.5);
N_inf = pow(pow(mu_p*(1-Hct),-0.25),0.5);

mu_casson = pow(mu_inf_casson,2)/C_STRAIN_RATE_MAG(c,t)+2*mu_inf_casson*N_inf/pow(C_STRAIN_RATE_MAG(c,t),0.5)+pow(N_inf,2);

return mu_casson;
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