import numpy as np 
c_mps =  299792458 # Speed of Light in meters per second
amu_to_MeV = 931.49432 # 1 a.m.u in MeV
GeVpc_to_kgmps = 5.36E-22 # covertsion of momentum from GeV/c to kg m/s
amu_to_kg = 1.66054E-27

def eV_to_Joule(E_eV):
    """Converts Energy from Electron Volts to Joules"""
    return E_eV*1.6*10**(-19)


def gamma_from_energy_MeVpu(EMevpu):
    gamma = EMevpu/amu_to_MeV + 1
    return gamma 

def beta_from_gamma(gamma):
    beta = np.sqrt(1 - 1/(gamma*gamma))
    return beta


def momentum_from_energy_MeVpu_and_massnumber(E_MeVpu, massNumber):
    A = massNumber
    E = E_MeVpu
    P = A*np.sqrt(E*E + 2*E*amu_to_MeV)
    return P

# -------------------------------- To be Completed ------------------------------------------- #
# TODO: 
def mm_to_mrad(len_mm):
    pass 
# TODO: Conversion of Twiss Parameters