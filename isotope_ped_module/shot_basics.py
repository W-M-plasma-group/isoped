# -*-Python-*-
# Created by chabanr at 13 Dec 2021  17:20

"""
This script will put the shot basics or some isotope shots into a dataset

defaultVars parameters
----------------------
:param kw1: kw1 can be passed to this script as <path to script>.run(kw1='hello')
"""

defaultVars()

shotInfo = {
    # Hydrogen Shots
    ## Gohil H mode threshold should work for large rho* in H
    # 133773:  (2150, 2550), missing tanh and ped values
    133778: (1500, 2500),
    133779: (2250, 2750),
    133784: (1500, 2700),
    133785: (2000, 2700),
    # 133786:  (4840, 5100), missing tanh and ped values
    133787: (2300, 2700),
    # Osborne pedestal scaling experiment
    # 133821 - 133849 good shots at 133831
    133826: (3700, 4400),
    133831: (2900, 3500),
    133832: (4000, 4600),
    133834: (3100, 3600),
    133835: (3100, 3700),
    133836: (3000, 3800),
    133838: (2600, 3400),
    133841: (3150, 3750),
    # Kathreen's day (183485 - 1835XX)
    183487: (5150, 5500),
    183489: (5100, 5600),
    183490: (5000, 5500),
    183491: (5050, 5600),  # elm Frequency very high
    183496: (3050, 3450),
    183497: (3050, 3450),
    # 183498:  (5000, 5600),  EFIT Problem?
    183499: (3200, 3500),
    # 183500:  (3150, 3450),  Also EFIT time distance problem
    # 183501:  (3350, 3850),
    # 183502:  (4100, 4400),
    183503: (3400, 4000),
    183504: (3200, 3800),
    # 183505:  (3400, 3950),  also efit
    183508: (4000, 4500),
    # 183509:  (3500, 4000), also efit
    ## RMP ELM Experiments you should never have included. High Plasma current
    # 183975: (2700, 3700),
    # 183976: (2700, 3700),
    # 183977: (3900, 4450),
    # 183979: (2200, 3700),
    # 183980: (2000, 3700),
    # 183981: (2300, 3600),
    # 183982: (2200, 3700),
    # 183984: (2200, 3700),
    # 183985: (2900, 3700),
    # 183986: (2000, 2800),
    # Livia Casali's Day
    184005: (1300, 1700),
    184006: (1400, 1700),
    184007: (2000, 3100),
    184013: (1300, 1700),
    184015: (2000, 3500),
    184016: (2500, 3800),
    184017: (2500, 3800),
    184018: (2500, 3800),
    # 184023: (4100, 5200),
    # 184024: (4400, 5600),
    # 184025: (4200, 5100),
    # 184026: (4200, 5600),
    # 184027: (4200, 5600),
    # 184028: (4200, 5600),
    # 184029: (4200, 5600),
    # 184030: (2100, 2700),
    # Deuterium shots:
    # Low trinagularity
    131996: (3600, 4000),
    131997: (3600, 4000),
    131999: (3800, 4000),
    132000: (3600, 4000),
    132002: (3500, 4200),
    # Ip & Bt scan at const q
    132003: (3500, 4800),
    132006: (3500, 4400),
    # 132007: theres a problem
    132008: (1200, 1600),  # locked mode at 2.2 sec
    # Upper triangularity scan ip=1.17 Bt = 2.1T
    132009: (2500, 3500),  # (3300, 4300) <- my guestimate of ss time
    132010: (3700, 4400),
    132011: (4150, 4800),
    # High traingularity
    132014: (2600, 3500),
    132016: (2650, 3200),
    132017: (2900, 3300),
    132018: (2900, 3300),
    # Ip scan at constant q
    132019: (1800, 2000),  # there are several dropouts weird shot.
    132020: (1500, 2100),
    # 132021: (2600, 3000),  # something weird here
    # Different day
    160245: (2500, 3000),
    160521: (2000, 3400),
    160522: (3200, 4200),
    160523: (3200, 4300),
    160524: (3400, 4400),
    160525: (3700, 4400),
    160526: (1800, 2700),
    160527: (3000, 3500),
    160528: (3900, 4500),
    160529: (1600, 2200),
    160530: (1600, 2200),
    160531: (1700, 2400),
    160532: (3300, 4000),
    # From Saskias notes all 2008 experiments in H
}

Hrange = [
    (133750, 133849),
    (157700, 160526),
    # (171477, 175901),
    (183487, 184045),
]

root['misc']['shotInfo'] = shotInfo
root['misc']['Hrange'] = Hrange