import Comove

targname="TYC 6110-1007-1"
radvel=-15.12

vlim=5.0
srad=25.0


output_location = Comove.findfriends(targname,radvel,velocity_limit=vlim,search_radius=srad,radec=[None,None],output_directory=None)
#print('output is in ' + output_location)
#print('Done')



##Example Target Names and radial velocities
    # Search targets
    #targname="HR 8799"
    #radvel = -13.0 * u.kilometer / u.second
    #targname="HD 110082"                        # Faux-Octans planet
    #radvel = 3.63 * u.kilometer / u.second
    #targname="V1298 Tau"
    #radvel = 16.15 * u.kilometer / u.second
    #targname="Kappa And"
    #radvel = -12.7 * u.kilometer / u.second
    #targname="HD 63433"                        # UMa planet
    #radvel = -15.81 * u.kilometer / u.second
    #targname='HD 107700'                        # Middle of Coma Ber cluster
    #radvel= 0.5 * u.kilometer / u.second
    #targname="HD 984"
    #radvel = 0.55 * u.kilometer / u.second
    #targname="TYC 7520-00369-1"
    #radvel = -2.33 * u.kilometer / u.second
    #targname="TYC 5909-00319-1"
    #radvel = 23.68 * u.kilometer / u.second
    #targname="HIP 73765"
    #radvel = -12.10 * u.kilometer / u.second
    #targname="K2-284"
    #radvel = 16.96 * u.kilometer / u.second
    #targname="Kepler-51"
    #radvel = -4.3 * u.kilometer / u.second
    #targname="TYC 4337-00128-1"
    #radvel = -9.21 * u.kilometer / u.second
    #targname="HD 92487"
    #radvel = 0.0 * u.kilometer / u.second # Guess for initial run, no RV available
    #targname="TYC 3496-1082-1"
    #radvel = -7.36 * u.kilometer / u.second
    #targname="HD 269132"                       # NGC 1901 central star
    #radvel = 1.43 * u.kilometer / u.second
    #targname="HIP 65469"
    #radvel = -7.33 * u.kilometer / u.second
    #targname="TYC 6110-1007-1"
    #radvel = -15.12 * u.kilometer / u.second
