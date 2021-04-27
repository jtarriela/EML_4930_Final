## checking blood model viscosity
import numpy as np
import math
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

        self.plot_experimental_visc_shear()

        return self.df

    def plot_experimental_visc_shear(self):
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
        plt.title('Blood Model Viscosity vs Shear Rate')
        plt.ylabel('Viscosity (Pa*sec)')
        plt.xlabel('Shear Rate (1/sec)')


        plt.show()

    def compute_sim_viscosity(self, df):
        sim_strain_df = df
        column_list = sim_strain_df.columns

        Parameters_BM_Class = {'Q_': Q_,
                               'mu_0': 0.056,
                               'mu_inf': 0.0035}
        bm = blood_model_validation(**Parameters_BM_Class)

        # sim_visc_df = sim_strain_df['Time [ s ]'].copy() ## Rename column to flowtime
        # sim_visc_df = sim_visc_df.rename("Flow Time", axis='columns')
        for i in range(len(column_list)):
            print(column_list[i])
            col_name = column_list[i]

            shear_time_list = sim_strain_df[col_name].tolist()

            if 'Casson' in str(col_name):  ##### Case Sensitive
                # print('Casson to be computed')
                for i in range(len(sim_strain_df)):
                    Parameters_Casson_1 = {'mu_p': Q_(0.00145, 'dimensionless'),
                                           'hematocrit': Q_(0.4, 'dimensionless'),
                                           'strain': shear_time_list[i]}

                    viscosity = bm.mod_casson_1(**Parameters_Casson_1)
                    sim_strain_df.loc[i, 'Casson Pa*sec'] = round(viscosity.magnitude, 8)

            elif 'CY' in str(col_name):
                # print('yay')
                for i in range(len(sim_strain_df)):
                    Parameters_Carreau = {'lambda_': Q_(3.313, 'second'),
                                          'n': Q_(0.3568, 'dimensionless'),
                                          'A': Q_(2, 'dimensionless'),
                                          'strain': shear_time_list[i]}

                    viscosity = bm.mod_carreau_yasuda(**Parameters_Carreau)
                    sim_strain_df.loc[i, 'CY Pa*sec'] = round(viscosity.magnitude, 8)


            elif 'Cross' in str(col_name):
                # print('cool')
                for i in range(len(sim_strain_df)):
                    Parameters_Cross = {'lambda_': Q_(1.007, 'second'),
                                        'm': Q_(1.028, 'dimensionless'),
                                        'strain': shear_time_list[i]}

                    viscosity = bm.mod_cross(**Parameters_Cross)
                    sim_strain_df.loc[i, 'Cross Pa*sec'] = round(viscosity.magnitude, 8)

        self.output_df = sim_strain_df
        return self.output_df

    def inlet_profile(self):
        self.inlet_profile_time = np.linspace(0, .5, 1000)
        self.inlet_profile_velocity = []
        self.inlet_profile_df = pd.DataFrame()

        for t in range(len(self.inlet_profile_time)):
            if self.inlet_profile_time[t] < 0.218:
                vel = 0.5*math.sin(4 * math.pi*(self.inlet_profile_time[t]+0.0160236))
            elif t > 0.218:
                vel = 0.1
            self.inlet_profile_velocity.append(vel)


        self.inlet_profile_df['Time (s)'] = self.inlet_profile_time
        self.inlet_profile_df['Velocity (m/s)'] = self.inlet_profile_velocity

        fig5 = plt.figure()
        plt.plot(self.inlet_profile_df['Time (s)'], self.inlet_profile_df['Velocity (m/s)'], 'c-', linewidth=2)
        plt.title('Inlet Velocity vs Time')
        plt.xlabel('Time (sec)')
        plt.ylabel('Velocity (m/sec)')
        plt.show()

        return self.inlet_profile_df

    def plot_wss(self):
        avg_wss = pd.read_csv('https://raw.githubusercontent.com/jtarriela/FDA_Blood_Pump/EML_4930/EML_4930/Data/Raw%20Data/avg_wss.csv')
        fig6 = plt.figure()
        plt.plot(avg_wss['Time [ s ]'], avg_wss['CY [ Pa ]'], 'y-', linewidth=2, label = "Carreau Yasuda")
        plt.plot(avg_wss['Time [ s ]'], avg_wss['Cassion [ Pa ]'], 'k-', linewidth=2, label = 'Cassion')
        plt.plot(avg_wss['Time [ s ]'], avg_wss['Cross [ Pa ]'], 'r-', linewidth=2, label = 'Cross')
        plt.plot(avg_wss['Time [ s ]'], avg_wss['Newtownian [ Pa ]'], 'g-', linewidth=2, label = 'Newtownian')
        plt.legend(loc="upper right", fontsize='medium', frameon=True, shadow=True)
        plt.title('Average WSS vs Time')
        plt.xlabel('Time (sec)')
        plt.ylabel('WSS (Pascal)')
        plt.show()

    def plot_mass_resid(self):
        mass_resid = pd.read_csv('https://raw.githubusercontent.com/jtarriela/FDA_Blood_Pump/EML_4930/EML_4930/Data/Raw%20Data/net_mass_flow_rate.csv')
        fig6 = plt.figure()
        plt.plot(mass_resid['Time [ s ]'], mass_resid['CY [ kg s^-1 ]'], 'y-', linewidth=2, label = "Carreau Yasuda")
        plt.plot(mass_resid['Time [ s ]'], mass_resid['Cassion [ kg s^-1 ]'], 'k-', linewidth=2, label = 'Cassion')
        plt.plot(mass_resid['Time [ s ]'], mass_resid['Cross [ kg s^-1 ]'], 'r-', linewidth=2, label = 'Cross')
        plt.plot(mass_resid['Time [ s ]'], mass_resid['Newtownian [ kg s^-1 ]'], 'g-', linewidth=2, label = 'Newtownian')
        plt.legend(loc="upper right", fontsize='medium', frameon=True, shadow=True)
        plt.title('Mass Residual vs Time')
        plt.xlabel('Time (sec)')
        plt.ylabel('Net Mass Flow Rate (kg/s)')
        plt.show()

    def plot_viscosity_time(self):
        fig2 = plt.figure()
        x_axis = self.output_df['Time [ s ]']
        plt.plot(x_axis, self.output_df['CY Pa*sec'], 'go--', linewidth=1, markersize=2, label="Carreau Yasuda")
        plt.plot(x_axis, self.output_df['Cross Pa*sec'], 'b--', linewidth=1, markersize=2, label="Cross")
        plt.plot(x_axis, self.output_df['Casson Pa*sec'], 'm--', linewidth=1, markersize=2, label="Casson (1)")
        plt.legend(loc="upper right", fontsize='medium', frameon=True, shadow=True)
        plt.title('Simulation Viscosity per Time Step')
        plt.ylabel('Viscosity (Pa*sec)')
        plt.xlabel('Time (sec)')
        plt.show()

    def plot_shear_time(self):
        fig3 = plt.figure()
        x_axis = self.output_df['Time [ s ]']
        plt.plot(x_axis, self.output_df['CY  [ s^-1 ]'], 'go-', linewidth=1, markersize=2, label="Carreau Yasuda")
        plt.plot(x_axis, self.output_df['Cross  [ s^-1 ]'], 'b-', linewidth=1, markersize=2, label="Cross")
        plt.plot(x_axis, self.output_df['Casson  [ s^-1 ]'], 'm-', linewidth=1, markersize=2, label="Casson (1)")
        plt.plot(x_axis, self.output_df['Newtownian  [ s^-1 ]'], 'c-', linewidth=1, markersize=2, label="Newtownian")
        plt.legend(loc="upper right", fontsize='medium', frameon=True, shadow=True)
        plt.title('Average Fluid Shear Rate per Time Step')
        plt.ylabel('Shear Rate (1/s)')
        plt.xlabel('Time (sec)')
        plt.show()

    def plot_viscosity_vs_shear(self):
        fig4 = plt.figure()
        x_axis = self.output_df['Time [ s ]']
        plt.plot(self.output_df['CY  [ s^-1 ]'], self.output_df['CY Pa*sec'], 'go', linewidth=1, markersize=5,
                 label="Carreau Yasuda")
        plt.plot(self.output_df['Cross  [ s^-1 ]'], self.output_df['Cross Pa*sec'], 'bo', linewidth=1, markersize=5,
                 label="Cross")
        plt.plot(self.output_df['Casson  [ s^-1 ]'], self.output_df['Casson Pa*sec'], 'mo', linewidth=1, markersize=5,
                 label="Casson (1)")
        plt.legend(loc="upper right", fontsize='medium', frameon=True, shadow=True)
        plt.title('Viscosity vs Shear Rate')
        plt.xlabel('Shear Rate (1/s)')
        plt.ylabel('Viscosity (Pa*s)')
        plt.show()



if __name__ == "__main__":
    # IN FLUENT, STRAIN RATE == SHEAR RATE
    # ALL REFERENCES TO STRAIN == SHEAR

    ##### DO NOT DELETE THIS CODE BLOCK #####
    ##### KWARG DEMO OF CLASS INPUT #####
    shear_range = np.arange(1,500)
    #
    Parameters_BM_Class = {'Q_': Q_,
                           'strain': 450, # Pass universal strain for all models
                           'strain_range': shear_range,
                           'mu_0': 0.056,
                           'mu_inf': 0.0035}
    bm = blood_model_validation(**Parameters_BM_Class)
    predicted_mu = bm.mu_calculations()
    #####--------------------------------#####

    # Parameters_Power = {'n': Q_(0.708, 'dimensionless'),
    #                     'k': Q_(0.017, 'dimensionless'),
    #                     'strain':0.05}
    # power=bm.mod_power(**Parameters_Power)
    # print(power)

    sim_strain_df = pd.read_csv('https://raw.githubusercontent.com/jtarriela/FDA_Blood_Pump/EML_4930/EML_4930/Data/Shear_stress_df.csv')

    sim_strain_visc_df = bm.compute_sim_viscosity(sim_strain_df)
    #
    # inlet_vel_df = bm.inlet_profile()
    # bm.plot_wss()
    # bm.plot_mass_resid()

    bm.plot_viscosity_time()