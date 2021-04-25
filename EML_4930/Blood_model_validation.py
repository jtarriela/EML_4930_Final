## checking blood model viscosity
import numpy as np
import math
import pint
from pint import UnitRegistry
ureg = UnitRegistry()
Q_ = ureg.Quantity
import pandas as pd
import matplotlib.pyplot as plt


class blood_model_validation():
    def __init__(self, Q_, mu_0,mu_inf,strain=None,strain_range=None):
        if strain is not None:
            self.strain = Q_(strain, '1/second')

        if strain_range is not None:
            self.strain_list = strain_range
            self.df = pd.DataFrame()
            self.df['Strain (Pascal)'] = pd.Series(self.strain_list)

        self.mu_0 = Q_(mu_0,'pascal*second')
        self.mu_inf = Q_(mu_inf, 'pascal*second')
        # super().__init__(Q_) # DONT FORGET TO COPY/PASTE THIS

        self.Parameters_Carreau = {'lambda_': Q_(3.313, 'second'),
                                   'n': Q_(0.3568, 'dimensionless'),
                                   'A': Q_(2, 'dimensionless')}

        self.Parameters_Cross = {'lambda_': Q_(1.007, 'second'),
                                 'm':Q_(1.028, 'dimensionless')}

        self.Parameters_Casson_1 = {'mu_p': Q_(0.00145, 'dimensionless'),
                                    'hematocrit':Q_(0.4, 'dimensionless')}

        self.Parameters_Casson_2 = {'yield_stress': Q_(0.005, 'pascal'),
                                    'n':Q_(0.0035, 'pascal*second')}

        self.Parameters_Power ={'n': Q_(0.708, 'dimensionless'),
                                'k': Q_(0.017, 'dimensionless')}

    def mod_carreau_yasuda(self, lambda_, n, A, strain = None):
        if strain is not None:
            self.strain = Q_(strain, '1/second')

        self.mu_cy = self.mu_inf + ((self.mu_0 - self.mu_inf))/math.pow((1+math.pow(lambda_*self.strain,A)),(1-n)/A)
        return self.mu_cy

    def mod_cross(self, lambda_, m, strain=None):
        if strain is not None:
            self.strain = Q_(strain, '1/second')

        self.mu_cr = self.mu_inf + ((self.mu_0 - self.mu_inf) / (1 + math.pow(lambda_ * self.strain, m)))
        return self.mu_cr

    def mod_casson_1(self, mu_p, hematocrit, strain = None):
        if strain is not None:
            strain = strain
        else:
            strain = self.strain * Q_('second')

        mu_inf_casson = np.sqrt(0.625*hematocrit)
        n_inf_casson = np.sqrt(mu_p*math.pow((1-hematocrit),-.25))
        self.mu_casson_1 = mu_inf_casson**2/strain + \
                         2*mu_inf_casson*n_inf_casson/np.sqrt(strain)+\
                         n_inf_casson**2

        return self.mu_casson_1

    def mod_casson_2(self,yield_stress,n, strain = None):
        if strain is not None:
            self.strain = Q_(strain, '1/second')

        self.mu_casson_2 = yield_stress/self.strain + \
                           np.sqrt(n*yield_stress)/np.sqrt(self.strain)+\
                           n

        return self.mu_casson_2

    def mod_power(self, n, k, strain):
        if strain is not None:
            self.strain = Q_(strain, '1/second')

        strain = self.strain*Q_('second')
        self.mu_power = k*math.pow(strain,(n-1))
        return self.mu_power

    def calc_all_one_point(self):

        model_cy = bm.mod_carreau_yasuda(**self.Parameters_Carreau)
        model_cross = bm.mod_cross(**self.Parameters_Cross)
        model_casson_1 = bm.mod_casson_1(**self.Parameters_Casson_1)
        model_casson_2 = bm.mod_casson_2(**self.Parameters_Casson_2)
        model_power = bm.mod_power(**self.Parameters_Power)

        print('Viscosity for Carreau Yasuda:{} \n'
              'viscosity for Cross Model: {} \n'
              'Viscosity for Casson Model 1: {} \n'
              'Viscosity for Casson Model 2 {} \n'
              'Viscosity for Power Model: {} \n'.format(model_cy,
                                                        model_cross,
                                                        model_casson_1,
                                                        model_casson_2,
                                                        model_power))

    def mu_calculations(self):

        temp = np.arange(0,len(self.strain_list))

        for i in temp:
            self.Parameters_Carreau.update({'strain': self.strain_list[i]})
            self.Parameters_Cross.update({'strain': self.strain_list[i]})
            self.Parameters_Casson_1.update({'strain': self.strain_list[i]})
            self.Parameters_Casson_2.update({'strain': self.strain_list[i]})
            self.Parameters_Power.update({'strain': self.strain_list[i]})

            self.mod_carreau_yasuda(**self.Parameters_Carreau)
            self.mod_cross(**self.Parameters_Cross)
            self.mod_casson_1(**self.Parameters_Casson_1)
            self.mod_casson_2(**self.Parameters_Casson_2)
            self.mod_power(**self.Parameters_Power)

            self.df.loc[i, 'Carreau Yasuda']=round(self.mu_cy.magnitude, 8)
            self.df.loc[i, 'Cross']=round(self.mu_cr.magnitude, 8)
            self.df.loc[i, 'Casson_1'] = round(self.mu_casson_1.magnitude, 8)
            self.df.loc[i, 'Casson_2'] = round(self.mu_casson_2.magnitude, 8)
            self.df.loc[i, 'Power'] = round(self.mu_power.magnitude, 8)

            # print(test)
        del self.Parameters_Carreau["strain"]
        del self.Parameters_Cross["strain"]
        del self.Parameters_Casson_1["strain"]
        del self.Parameters_Casson_2["strain"]
        del self.Parameters_Power["strain"]

        self.plot_()

        return self.df

    def plot_(self):
        fig1 = plt.figure()
        cy = plt.plot(self.df['Strain (Pascal)'], self.df['Carreau Yasuda'], 'go--', linewidth=1,
                      markersize=2, label="Carreau Yasuda")
        cr = plt.plot(self.df['Strain (Pascal)'], self.df['Cross'], 'bo-', linewidth=1, markersize=2,
                      label="Cross")
        ca1 = plt.plot(self.df['Strain (Pascal)'], self.df['Casson_1'], 'm--', linewidth=1, markersize=2,
                       label="Casson (1)")
        ca2 = plt.plot(self.df['Strain (Pascal)'], self.df['Casson_2'], '^k', linewidth=1, markersize=2,
                       label="Casson (2)")
        pw = plt.plot(self.df['Strain (Pascal)'], self.df['Power'], 'ro--', linewidth=1, markersize=2,
                      label="Power")

        plt.legend(loc="upper right", fontsize='medium', frameon=True, shadow=True)
        plt.title('Blood Model Viscosity vs Shear Stress')
        plt.ylabel('Viscosity (Pa*sec)')
        plt.xlabel('Shear Rate (1/sec')


        plt.show()
    #
    # def plot_est_viscosity_range(self):

if __name__ == "__main__":

    # strain pass through individual methods as unit of pressure
    # will convert to 1/second

    shear_range = np.arange(1,500)

    Parameters_BM_Class = {'Q_': Q_,
                           'strain': 450, # Pass universal strain for all models
                           'strain_range': shear_range,
                           'mu_0': 0.056,
                           'mu_inf': 0.0035}
    bm = blood_model_validation(**Parameters_BM_Class)
    predicted_mu = bm.mu_calculations()

    # Parameters_Power = {'n': Q_(0.708, 'dimensionless'),
    #                     'k': Q_(0.017, 'dimensionless'),
    #                     'strain':0.05}
    # power=bm.mod_power(**Parameters_Power)
    # print(power)

    # sim_strain_df = pd.read_csv(r'C:\Users\Joseph Tarriela\OneDrive - University of South Florida\School Year\2020_21\Spring 2021\Final Project\Post-Exports\Shear_stress_df.xlsx')

