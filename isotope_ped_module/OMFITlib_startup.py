# -*-Python-*-


# define constants
from scipy import constants

speed_of_light = constants.speed_of_light
electron_charge = constants.e
proton_mass_joules = constants.physical_constants["proton mass energy equivalent"][0]  # in JOULES
proton_mass_kg = constants.proton_mass
electron_mass_joules = constants.physical_constants["electron mass energy equivalent"][0]  # in Joules
electron_mass_kg = constants.electron_mass
electron_mass_ev = constants.physical_constants["electron mass energy equivalent in MeV"][0] / 1e6  # in eV
boltzmann_const_JperK = constants.physical_constants["Boltzmann constant"][0]
# boltzmann_const_JperEV = constants.physical_constants[""]
electron_potential = 13.6
ion_potential = 6.8
eV2Kelvin = 1.16e4
Kelvin2eV = eV2Kelvin ** -1


# Gas constant to convert Tom's tools pressures to SI and atm units
R_pascals = 8.314  # m^3 * Pa * K^-1 * mol^-1

# ne_osborne 10^20 * te_osborne 10^3
# R_pascals * eV2kelvin / Avogadro * oom_correction
R_osborne2Pa = R_pascals / 6.022e23 * 11604 * 10e23


### THIS SECTION IS USED TO DEFINE COMPOUND SIGNALS ###
compound_func_dict = SortedDict()  # dictionary iterated through at the end


def calc_PTOT(dfIN):
    # Total Power in Ptot
    PTOT = np.array(dfIN['PNBI'] + dfIN['PECH'] + dfIN['POH'])
    PTOT_ERR = np.sqrt(np.array(dfIN['PNBI_ERR'] ** 2 + dfIN['PECH_ERR'] ** 2 + dfIN['POH_ERR'] ** 2))

    return PTOT, PTOT_ERR


compound_func_dict['PTOT'] = calc_PTOT


def calc_PCROSS(dfIN):
    # Power crossing separatrix PCROSS
    PCROSS = np.array(dfIN['PNBI'] + dfIN['PECH'] + dfIN['POH']) - dfIN['PRAD_CORE']
    # PCROSS_ERR = np.sqrt(dfIN['PTOT_ERR'] ** 2 + dfIN['PRAD_CORE_ERR'] ** 2)
    return PCROSS, np.zeros_like(PCROSS)


compound_func_dict['PCROSS'] = calc_PCROSS


def calc_GASTOT(dfIN):
    # Total gas
    GASTOT = np.array(dfIN['GASA_CAL'] + dfIN['GASB_CAL'])
    GASTOT_ERR = np.array(dfIN['GASA_CAL_ERR'] ** 2 + dfIN['GASB_CAL_ERR'] ** 2)
    return GASTOT, GASTOT_ERR


compound_func_dict['GASTOT'] = calc_GASTOT


def calc_NUSTAR(dfIN):
    # from Theresa Wilks collisionality area
    teped, neped = dfIN['PRMTAN_TEPED'].values, dfIN['PRMTAN_NEPED'].values
    q95, R0, aminor = dfIN['Q95'].values, dfIN['R0'].values, dfIN['AMINOR'].values
    # cerqnzt1 = df["CERQNZT1"].values
    # print(df['CERQNZT23'])
    cerqnzt23 = dfIN["CERQNZT23"].values
    # # alog is IDL for natural logarithm
    # df['NUSTAR'] = 6.92e-18 * q95 * R0 * neped * (1 + 30 * cerqnzt1 / neped) * 31.3 - np.log(neped * 0.5 / teped) / (
    #     (teped ** 2 * (aminor / R0) ** 1.5)
    # )
    # # From elm_control2.revplus (saskias)
    NUSTAR = (
        6.92e-18
        * q95
        * R0
        * neped
        * (1 + 30 * (cerqnzt23 / neped))
        * (31.3 - np.log(neped ** 0.5 / teped))
        / (teped ** 2 * (aminor / R0) ** 1.5)
    )

    return NUSTAR, np.zeros_like(NUSTAR)


compound_func_dict['NUSTAR'] = calc_NUSTAR


def calc_logNUSTAR(dfIN):
    # this is for color plotting
    return np.log(dfIN['NUSTAR']), np.zeros_like(dfIN.index)


compound_func_dict['logNUSTAR'] = calc_logNUSTAR


def calc_RHOSTARE(dfIN):
    # rho* electrons
    num = np.sqrt(electron_mass_ev * dfIN['PRMPED_TEPED'].values) * eV2Kelvin * boltzmann_const_JperK
    den = dfIN['A95'].values * dfIN['BTOT95'].values * electron_charge
    # print("rho*", num)
    # print("rho*", den)
    RHOSTARE = num / den
    return RHOSTARE, np.zeros_like(RHOSTARE)


compound_func_dict['RHOSTARE'] = calc_RHOSTARE


def calc_SYM_SHIFT(dfIN):
    # shift in symmetry points
    SYM_SHIFT = dfIN['PRMTAN_NESYM'].values - dfIN["PRMTAN_TESYM"].values
    return SYM_SHIFT, np.zeros_like(SYM_SHIFT)


compound_func_dict['SYM_SHIFT'] = calc_SYM_SHIFT


def calc_L_NE(dfIN):
    # Gradient Scale length "L_NE" = ne/grad(ne)
    grad_ne = (dfIN['PRMTAN_NEPED'].values - dfIN['PRMPED_NESEP'].values) / dfIN["PRMTAN_NEWID"].values
    ne_avg = (dfIN['PRMTAN_NEPED'].values + dfIN['PRMPED_NESEP'].values) / 2
    L_NE = ne_avg / grad_ne
    return L_NE, np.zeros_like(L_NE)


compound_func_dict['L_NE'] = calc_L_NE


def calc_L_TE(dfIN):
    # Temperature Gradient Scale length "L_TE" = Te/grad(Te)
    grad_Te = (dfIN['PRMTAN_TEPED'].values - dfIN['PRMPED_TESEP'].values) / dfIN["PRMTAN_TEWID"].values
    Te_avg = (dfIN['PRMTAN_TEPED'].values + dfIN['PRMPED_TESEP'].values) / 2
    L_TE = Te_avg / grad_Te
    return L_TE, np.zeros_like(L_TE)


compound_func_dict['L_TE'] = calc_L_TE


def calc_ETAE(dfIN):
    # Eta_e = Lne/Lte
    ETAE = dfIN['L_NE'].values / dfIN['L_TE'].values
    return ETAE, np.zeros_like(ETAE)


compound_func_dict['ETAE'] = calc_ETAE


def calc_GRADPE(dfIN):
    PESEP = dfIN['PRMPED_TESEP'].values * dfIN['PRMPED_NESEP'].values
    GRADPE = -(dfIN['PRMTAN_PEPED'].values - PESEP) / dfIN["PRMTAN_PEWID"].values
    return GRADPE, np.zeros_like(GRADPE)


compound_func_dict['GRADPE'] = calc_GRADPE


def calc_GRADPE_NORM(dfIN):
    GRADPE = -(dfIN['PRMTAN_PEPED'].values - (dfIN['PRMPED_TESEP'].values * dfIN['PRMPED_NESEP'].values)) / dfIN["PRMTAN_PEWID"].values
    NORM = -2 * 1e-7 * (dfIN['R0'] + dfIN['AMINOR'] - 0.06) * dfIN['Q95'] ** 2 / dfIN['BTOT95'] ** 2
    return GRADPE / NORM, np.zeros_like(GRADPE)


compound_func_dict['GRADPE_NORM'] = calc_GRADPE_NORM


def calc_L_PE(dfIN):
    GRADPE = -(dfIN['PRMTAN_PEPED'].values - (dfIN['PRMPED_TESEP'].values * dfIN['PRMPED_NESEP'].values)) / dfIN["PRMTAN_PEWID"].values
    PE_AVG = (dfIN['PRMTAN_PEPED'].values + (dfIN['PRMPED_TESEP'].values * dfIN['PRMPED_NESEP'].values)) / 2
    return PE_AVG / GRADPE, np.zeros_like(GRADPE)


compound_func_dict['L_PE'] = calc_L_PE


def calc_RATIO_NE_SEPPED(dfIN):
    return dfIN['PRMPED_NESEP'] / dfIN['PRMTAN_NEPED'], np.zeros_like(dfIN.index)


compound_func_dict['RATIO_NE_SEPPED'] = calc_RATIO_NE_SEPPED


def calc_RATIO_PE_SEPPED(dfIN):
    PESEP = dfIN['PRMPED_TESEP'] * dfIN['PRMPED_NESEP']
    PEPED = dfIN['PRMPED_TEPED'] * dfIN['PRMPED_NEPED']
    return PESEP / PEPED, np.zeros_like(dfIN.index)


compound_func_dict['RATIO_PE_SEPPED'] = calc_RATIO_PE_SEPPED


def calc_PRMPED_PEPED(dfIN):
    # Calculating this because the pedestal pressure from the TAN tree is in some weird units. possibly MPA?
    return dfIN['PRMPED_TEPED'] * dfIN['PRMPED_NEPED'], np.zeros_like(dfIN.index)


compound_func_dict['PRMPED_PEPED'] = calc_PRMPED_PEPED

### Special signals
special_signal_analysis = SortedDict()


def calc_ELMFREQUENCY(time, data):
    # factor of 1000 is for ms -> s
    return 'ELMFREQUENCY', 1000 * np.count_nonzero(data) / (time[-1] - time[0]), 0.0  # second value is error


special_signal_analysis['ELMSMARKER'] = calc_ELMFREQUENCY

print("Finished loading OMFITlib_startup")
