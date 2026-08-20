### missing parameters

# in resistances.R
rbs =  soil[['rbs']];
L =  meteo[['L']];

# in heatfluxes.R resolve
#es_fun saturated pressure function es(hPa)=f(T(C))
#s_fun slope of the saturated pressure function (s(hPa/C) = f(T(C), es(hPa))

# iin ebal.R
# input preparation
nl          = canopy.nlayers;
GAM         = soil.GAM;
Ps          = gap.Ps;
kV          = canopy.kV;
xl          = canopy.xl;
rss = soil.rss;
# in bichemical MD12.R
Temp = meteo[['Tleaf']];
eb      = meteo[['eb']]
Rdparam =  leafbio[['RdPerVcmax25']] ##Check if its correct, original indicate RdPerVcmax25;

# in bichemical.R
Q = Cs = meteo[['Q']]                                                  # [umol m-2 s-1] absorbed PAR flux
Cs = meteo[['Cs']]
