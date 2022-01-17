# -*-Python-*-
# Created by chabanr at 14 Dec 2021  13:39

"""
This script is used to populate the pandas dataframe in OUTPUTS['df'] with signals

defaultVars parameters
----------------------
:param kw1: kw1 can be passed to this script as <path to script>.run(kw1='hello')
"""

defaultVars(useless_rn=True)

skip_datagrab = False

### Database messing from saskia meeting 12/13/2021
import pandas as pd
from matplotlib import cm

from OMFITlib_startup import special_signal_analysis, compound_func_dict

tst_shots = shotInfo = {
    133826: (3700, 4400),
    133831: (2900, 3500),
    133832: (4000, 4600),
}

shotInfo = root['misc']['shotInfo']
Hrange = root['misc']['Hrange']

default_signals = [
    "IP",
    "BT",
    "DENSITY",
    "Q95",
    "BCENTR",
    "R0",
    "AMINOR",  # basics
    # "BETAN",  # for some reason not working regardless of case adding from EFIT
    "TRITOP",
    "PNBI",
    "PECH",
    "POH",
    "PRAD_TOT",
    "PRAD_CORE",
    # pedestal stuff
    "PRMPED_TEPED",
    "PRMPED_NEPED",
    "PRMPED_NESEP",
    "PRMPED_TESEP",
    "PRMTAN_TEPED",
    "PRMTAN_TEWID",
    "PRMTAN_TESYM",
    "PRMTAN_NEPED",
    "PRMTAN_NESYM",
    "PRMTAN_NEWID",
    # "PRMTAN_PESEP", There is no separatrix pressure anywhere
    "PRMTAN_PEWID",
    "PRMTAN_PEPED",
    "PRMTAN_PESYM",
    "PRMTAN_NELOWER",
    # "PRMPED_CHI", # heat flux at separatrix currently not working
    # "PRMPED_IPERPPED", # heat flux at pedestal top aso not working
    # "PRMPED_IPERP",  # pedestal  also not working
    # "RHOSTAR"  # dimensionless
    "GASA_CAL",
    "GASB_CAL",  # gasses
    "CERQNZT1",
    'CERQNZT23',
    'ASDEX_LOBAF',  # pressure gauges This doesnt work!
    'ELMSMARKER'  # where the elms happen
    #'ASDEX_MID', # These dont work
    #'ASDEX_PRI', # These dont work
]

# going to write a skip script
# to make it grab all the data delete the df in OUTPUTS
try:
    df = root['OUTPUTS']['df']
    shots_to_get = [x for x in shotInfo.keys() if x not in df.index]
except Exception:
    print("No dataframe exists in outputs. Collecting all data")
    shots_to_get = shotInfo.keys()

if skip_datagrab:  # if skipping the datagrab then just make the iterable empty
    shots_to_get = []
    # you already grabbed the df so no need to do it here

data = SortedDict()
empty_data = {sig: None for sig in default_signals}
empty_data.update({sig + '_ERR': None for sig in default_signals})

for shot in shots_to_get:
    tstart, tend = shotInfo[shot]
    tmid = int(np.mean([tstart, tend]))
    data[shot] = empty_data.copy()
    data[shot]['SHOT'] = shot
    # data_err[shot] = empty_data.copy()

    # gather the efit01
    efitout = from_mds_plus(device='DIII-D', shot=shot, times=[tmid], time_diff_warning_threshold=26, get_afile=False)
    gfile = efitout['gEQDSK'][tmid]
    ind95 = np.argmin(np.abs(gfile['AuxQuantities']['PSI_NORM'] - 0.95))

    # using flux 95 as the "top of the pedestal"
    data[shot]['BETAN'] = gfile['fluxSurfaces']['avg']['beta_n']
    data[shot]['BETAN95'] = gfile['fluxSurfaces']['avg']['beta_n'][ind95]
    data[shot]['BTOT95'] = gfile['fluxSurfaces']['avg']['Btot'][ind95]
    data[shot]['R0'] = gfile['fluxSurfaces']['R0']
    data[shot]['A95'] = gfile['fluxSurfaces']['avg']['a'][ind95]

    for signame in default_signals:
        try:
            # print(f"{shot}: {signame}")
            out = OMFITmdsValue(server='DIII-D', shot=shot, TDI=signame)
            time = out.dim_of(0)
            dat = out.data()
        except Exception:
            print(f"Failed data retrieval on {shot}: {signame}")
            dat = np.nan * np.ones(10)
            time = np.nan * np.ones(10)

        # Get the mean and std of the trace value
        try:
            my_ind = np.where((time >= tstart) & (time <= tend))[0]

            if signame in special_signal_analysis.keys():  # if it has a special analysis other than mean do it here
                special_name, special_out, special_out_err = special_signal_analysis[signame](time[my_ind], dat[my_ind])
                data[shot][special_name], data[shot][special_name + '_ERR'] = special_out, special_out_err
            else:
                data[shot][signame] = np.mean(dat[my_ind])
                data[shot][signame + '_ERR'] = np.nanstd(dat[my_ind] / 1e19) * 1e19

        except Exception:
            print(f"ERROR: {shot}: {signame}. Setting to None")
            data[shot][signame] = None
            data[shot][signame + '_ERR'] = None

    # ENDFOR SIGNAL

    # Get the Isotope
    shot_in_a_range = [shot in range(r[0], r[1] + 1) for r in Hrange]
    shot_is_hydrogen = any(shot_in_a_range)
    print(f"Shot in range: {shot_in_a_range}")
    if shot_is_hydrogen:
        print(f"{shot} in is Hydrogen")
        data[shot]['ISOTOPE'] = 'H'
        data[shot]['MASS'] = 1
    else:
        print(f"{shot} is Deuterium")
        data[shot]['ISOTOPE'] = 'D'
        data[shot]['MASS'] = 2

    # put the SS time in the dataframe
    data[shot]['SSTIME'] = shotInfo[shot]

# ENDFOR SHOT

if skip_datagrab or (len(shots_to_get) != 0):  # just put this here for ease
    # write the simple signals to a dataframe
    data_list = [data[shot] for shot in shotInfo.keys()]
    shot_list = shotInfo.keys()
    df = pd.DataFrame(data=data_list, index=shot_list)
    df.set_index('SHOT')

print(f"THere are a total of \n{df['ISOTOPE'].value_counts()}")

# Run through and calculate compound signals
# Run through the compound signals:
print("Calculating the Compound signals")
for compound_signal in compound_func_dict.keys():
    try:
        df[compound_signal], df[compound_signal + '_ERR'] = compound_func_dict[compound_signal](df)
    except Exception as e:
        print(f"Exception in calcualting compound signal: {compound_signal}")
        print(e)
        print(f"Moving on... \n")

# ENDFOR Compound Signal loop


root['OUTPUTS']['df'] = df
# root['OUTPUTS']['df_err'] = df_err

# dfH = df.where(df['ISOTOPE'] == 'H')
# dfD = df.where(df['ISOTOPE'] == 'D')

# if plot:
#     root['PLOTS']['plot_df'].run(abscissa=['PRMPED_NEPED'], ordinate=['PRMPED_NESEP'], colorby=['PTOT'])
