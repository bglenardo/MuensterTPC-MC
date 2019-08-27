import numpy as np
from scipy.special import gamma as gamma_func
import pandas as pd

m_e = 0.511 # MeV/c^2
alpha = 1./137.
hbar = 197.3 # MeV*fm (really this is hbar*c, assuming c=1)

class Geoneutrinos:

    def __init__( self ):

        self.U238_dict = dict()
        self.Th232_dict = dict()
        self.K40_dict = dict()

    def nu_beta_spectrum( self, Q, A, Z, E_nu ):
        W_e = Q - E_nu + m_e
        p_e = np.sqrt(Q - E_nu) * np.sqrt(Q - E_nu + 2.*m_e)
        return p_e * W_e * E_nu**2 * self.fermi_function(Z,A,W_e,p_e)
    
    def fermi_function( self, Z, A, W_e, p_e ):
        R = 1.2 * A**(1./3.) #fm
        eta = alpha * Z * W_e / p_e
        gamma = np.sqrt(1. - (alpha*Z)**2)
    
        return 4 * (2*p_e*R/hbar)**(2.*gamma-2.) * np.exp( np.pi * eta ) *\
            np.absolute( gamma_func( gamma + 1j * eta ) )**2 / \
            np.absolute( gamma_func( 2*gamma + 1 ) )**2


    def LoadData(self, spreadsheet):

        self.df_U238_branching = pd.read_excel('Geo_neutrinos/Th_U_K decay betas.xlsx',\
                                          sheet_name='U238 chain branching')
        self.df_Th232_branching = pd.read_excel('Geo_neutrinos/Th_U_K decay betas.xlsx',\
                                           sheet_name='Th232 chain branching')
        self.df_K40_branching = pd.read_excel('Geo_neutrinos/Th_U_K decay betas.xlsx',\
                                         sheet_name='K40 chain branching')
        
        for index,row in self.df_U238_branching.iterrows():
            print('Adding {}'.format(row['Isotope ']))
            self.U238_dict.update( {row['Isotope ']: \
                 pd.read_excel(spreadsheet,sheet_name=row['Isotope '])} )
        for index,row in self.df_Th232_branching.iterrows():
            print('Adding {}'.format(row['Isotope']))
            self.Th232_dict.update( {row['Isotope']: \
                       pd.read_excel(spreadsheet,sheet_name=row['Isotope'])} )
        for index,row in self.df_K40_branching.iterrows():
            print('Adding {}'.format(row['Isotope']))
            self.K40_dict.update( {row['Isotope']: \
                       pd.read_excel(spreadsheet,sheet_name=row['Isotope'])} )    
    


    def U238_spectrum( self, E_nu ):
        if not self.U238_dict: 
            print('ERROR: no data loaded.')
            return
        total_spec = np.zeros(len(E_nu))
        dE_nu = E_nu[2]-E_nu[1]
        specs=[]
        for index,row in self.df_U238_branching.iterrows():
            for subindex,subrow in self.U238_dict[ row['Isotope '] ].iterrows():
                this_spec = np.zeros(len(E_nu))
                mask = E_nu<(subrow['End-point energy\n(keV)']/1.e3)
                #print('Isotope: {}\t Endpoint: {}\tA: {}\tZ: {}'.format(row['Isotope '],\
                #                                               subrow['End-point energy\n(keV)']/1.e3,\
                #                                               row['A'],\
                #                                               row['Z']))
                this_spec[mask] = self.nu_beta_spectrum( subrow['End-point energy\n(keV)']/1.e3,\
                                                   row['A'],\
                                                   row['Z'],\
                                                   E_nu[mask])
                if np.sum(this_spec>0.):
                    this_spec = this_spec / np.sum(this_spec*dE_nu) * \
                        row['Branching fraction (%)']/100. * subrow['Intensity\n(%)']/100.

                specs.append(this_spec)
                total_spec = total_spec + this_spec
        self.U238_spec = total_spec
        return specs, total_spec
    
    
    
    
    
    def Th232_spectrum( self, E_nu ):
        if not self.Th232_dict: 
            print('ERROR: no data loaded.')
            return
        total_spec = np.zeros(len(E_nu))
        dE_nu = E_nu[2]-E_nu[1]
        specs=[]
        for index,row in self.df_Th232_branching.iterrows():
            for subindex,subrow in self.Th232_dict[ row['Isotope'] ].iterrows():
                this_spec = np.zeros(len(E_nu))
                mask = E_nu<(subrow['End-point energy\n(keV)']/1.e3)
                #print('Isotope: {}\t Endpoint: {}\tA: {}\tZ: {}'.format(row['Isotope'],\
                #                                               subrow['End-point energy\n(keV)']/1.e3,\
                #                                               row['A'],\
                #                                               row['Z']))
                this_spec[mask] = self.nu_beta_spectrum( subrow['End-point energy\n(keV)']/1.e3,\
                                                   row['A'],\
                                                   row['Z'],\
                                                   E_nu[mask])
                if np.sum(this_spec>0.):
                    this_spec = this_spec / np.sum(this_spec*dE_nu) * \
                        row['Branching fraction (%)']/100. * subrow['Intensity\n(%)']/100.

                specs.append(this_spec)
                total_spec = total_spec + this_spec
        self.Th232_spec = total_spec
        return specs, total_spec
    
    
    
    def K40_spectrum( self, E_nu ):
        if not self.K40_dict: 
            print('ERROR: no data loaded.')
            return
        total_spec = np.zeros(len(E_nu))
        dE_nu = E_nu[2]-E_nu[1]
        specs=[]
        for index,row in self.df_K40_branching.iterrows():
            for subindex,subrow in self.K40_dict[ row['Isotope'] ].iterrows():
                this_spec = np.zeros(len(E_nu))
                mask = E_nu<(subrow['End-point energy\n(keV)']/1.e3)
                #print('Isotope: {}\t Endpoint: {}\tA: {}\tZ: {}'.format(row['Isotope'],\
                #                                               subrow['End-point energy\n(keV)']/1.e3,\
                #                                               row['A'],\
                #                                               row['Z']))
                this_spec[mask] = self.nu_beta_spectrum( subrow['End-point energy\n(keV)']/1.e3,\
                                                   row['A'],\
                                                   row['Z'],\
                                                   E_nu[mask])
                #print(this_spec)
                if np.sum(this_spec>0.):
                    this_spec = this_spec / np.sum(this_spec*dE_nu) * \
                        row['Branching fraction (%)']/100. * subrow['Intensity\n(%)']/100.

                specs.append(this_spec)
                total_spec = total_spec + this_spec
        self.K40_spec = total_spec
        return specs, total_spec


    