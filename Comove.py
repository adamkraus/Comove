import pkg_resources
##figure out where the big fits files are in this installation
datapath = pkg_resources.resource_filename('Comove','resources')
import math as math
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
Simbad.reset_votable_fields()
Simbad.TIMEOUT = 15
Simbad.server = "simbad.harvard.edu"
Simbad.add_votable_fields('typed_id')
customSimbad = Simbad()
customSimbad.add_votable_fields('rvz_radvel','rvz_error','rvz_bibcode')
from astropy.coordinates import SkyCoord
from astropy import coordinates
from astropy.coordinates import ICRS
from astroquery.gaia import Gaia
from astroquery.exceptions import NoResultsWarning
import galpy.util.bovy_coords as bc
import matplotlib.pyplot as plt
from matplotlib import cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from astroquery.mast import Catalogs
from astroquery.irsa import Irsa
Irsa.TIMEOUT = 600
from astropy.coordinates import SkyCoord
from scipy.interpolate import interp1d
from scipy.io.idl import readsav
from astroquery.vizier import Vizier
from astropy.utils.data import conf
conf.remote_timeout = 60.0
Vizier.TIMEOUT = 600
import os,warnings,sys
import urllib.request
import csv
import pickle
import matplotlib as mpl

if 'dustmaps.bayestar' in sys.modules: print('Bayestar already imported, skipping 30-second load time.')

if 'dustmaps.bayestar' not in sys.modules:
    print('Bayestar not imported, doing so now. Will require 30 seconds or so.')
    from dustmaps.config import config
    datadir = '~/Dropbox/Malmquist/'
    bayestarver = 'bayestar2019'
    testname = datadir + 'bayestar/' + bayestarver + '.h5'
    config['data_dir'] = datadir
    from dustmaps.bayestar import BayestarQuery
    bayestar = BayestarQuery(version=bayestarver)
    if ((os.path.isfile(os.path.expanduser(testname))) == True) : print('Already downloaded Bayestar files.')
    if ((os.path.isfile(os.path.expanduser(testname))) == False):
        import dustmaps.bayestar # Only uncomment if running in a new place to download dust maps again
        dustmaps.bayestar.fetch()

mpl.rcParams['lines.linewidth']   = 2
mpl.rcParams['axes.linewidth']    = 2
mpl.rcParams['xtick.major.width'] =2
mpl.rcParams['ytick.major.width'] =2
mpl.rcParams['ytick.labelsize'] = 10
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['legend.numpoints'] = 1
mpl.rcParams['axes.labelweight']='semibold'
mpl.rcParams['axes.titlesize']=9
mpl.rcParams['axes.titleweight']='semibold'
mpl.rcParams['font.weight'] = 'semibold'
plt.rcParams['figure.facecolor'] = 'white'

def findfriends(targname,radial_velocity,velocity_limit=5.0,search_radius=25.0,rvcut=5.0,convergcut=5.0,radec=[None,None],output_directory = None,showplots=False,verbose=False,DoGALEX=True,DoWISE=True,DoROSAT=True):


    radvel= radial_velocity * u.kilometer / u.second

    if (convergcut == None): convergcut = 0.0

    if output_directory == None:
        outdir = './' + targname.replace(" ", "") + '_friends/'
    else: 
        outdir = output_directory
    if os.path.isdir(outdir) == True:
        print('Output directory ' + outdir +' Already Exists!!')
        print('Either Move it, Delete it, or input a different [output_directory] Please!')
        return
    os.mkdir(outdir)

    if velocity_limit < 0.00001 :
        print('input velocity_limit is too small, try something else')
        print('velocity_limit: ' + str(velocity_limit))
    if search_radius < 0.0000001:
        print('input search_radius is too small, try something else')
        print('search_radius: ' + str(search_radius))

    # Search parameters
    vlim=velocity_limit * u.kilometer / u.second
    searchradpc=search_radius * u.parsec

    if (radec[0] != None) & (radec[1] != None):
        usera,usedec = radec[0],radec[1]
    else:  ##use the target name to get simbad ra and dec.
        print('Asking Simbad for RA and DEC')
        result_table = Simbad.query_object(targname)
        usera,usedec = result_table['RA'][0],result_table['DEC'][0]
    
    if verbose == True:
        print('Target name: ',targname)
        print('Coordinates: ' + str(usera) +' '+str(usedec))
        print()

    c = SkyCoord( ra=usera , dec=usedec , unit=(u.hourangle, u.deg) , frame='icrs')
    if verbose == True: print(c)

    # Find precise coordinates and distance from Gaia, define search radius and parallax cutoff
    print('Asking Gaia for precise coordinates')
    sqltext = "SELECT * FROM gaiadr3.gaia_source WHERE CONTAINS( \
               POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec), \
               CIRCLE('ICRS'," + str(c.ra.value) +","+ str(c.dec.value) +","+ str(6.0/3600.0) +"))=1;"
    job = Gaia.launch_job_async(sqltext , dump_to_file=False)
    Pgaia = job.get_results()
    if verbose == True:
        print(sqltext)
        print()
        print(Pgaia['source_id','ra','dec','phot_g_mean_mag','parallax','ruwe'].pprint_all())
        print()
        print(Pgaia['phot_g_mean_mag'].mask)

    minpos = Pgaia['phot_g_mean_mag'].tolist().index(min( Pgaia['phot_g_mean_mag'][~Pgaia['phot_g_mean_mag'].mask] ))

    Pcoord = SkyCoord( ra=Pgaia['ra'][minpos]*u.deg , dec=Pgaia['dec'][minpos]*u.deg , \
                      distance=(1000.0/Pgaia['parallax'][minpos])*u.parsec , frame='icrs' , \
                      radial_velocity=radvel , \
                      pm_ra_cosdec=Pgaia['pmra'][minpos]*u.mas/u.year , pm_dec=Pgaia['pmdec'][minpos]*u.mas/u.year )

    searchraddeg = np.arcsin(searchradpc/Pcoord.distance).to(u.deg)
    minpar = (1000.0 * u.parsec) / (Pcoord.distance + searchradpc) * u.mas
    if verbose == True:
        print(Pcoord)
        print()
        print('Search radius in deg: ',searchraddeg)
        print('Minimum parallax: ',minpar)


    # Query Gaia with search radius and parallax cut
    # Note, a cut on parallax_error was added because searches at low galactic latitude 
    # return an overwhelming number of noisy sources that scatter into the search volume - ALK 20210325
    print('Querying Gaia for neighbors')

    Pllbb     = bc.radec_to_lb(Pcoord.ra.value , Pcoord.dec.value , degree=True)
    if ( np.abs(Pllbb[1]) > 10.0): plxcut = max( 0.5 , (1000.0/Pcoord.distance.value/10.0) )
    else: plxcut = 0.5
    print('Parallax cut: ',plxcut)

    if (searchradpc < Pcoord.distance):
        sqltext = "SELECT * FROM gaiadr3.gaia_source WHERE CONTAINS( \
            POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec), \
            CIRCLE('ICRS'," + str(Pcoord.ra.value) +","+ str(Pcoord.dec.value) +","+ str(searchraddeg.value) +"))\
            =1 AND parallax>" + str(minpar.value) + " AND parallax_error<" + str(plxcut) + ";"
    if (searchradpc >= Pcoord.distance):
        sqltext = "SELECT * FROM gaiadr3.gaia_source WHERE parallax>" + str(minpar.value) + " AND parallax_error<" + str(plxcut) + ";"
        print('Note, using all-sky search')
    if verbose == True:
        print(sqltext)
        print()

    job = Gaia.launch_job_async(sqltext , dump_to_file=False)
    r = job.get_results()
   
    if verbose == True: print('Number of records: ',len(r['ra']))


    # Construct coordinates array for all stars returned in cone search

    gaiacoord = SkyCoord( ra=r['ra'] , dec=r['dec'] , distance=(1000.0/r['parallax'])*u.parsec , \
                         frame='icrs' , \
                         pm_ra_cosdec=r['pmra'] , pm_dec=r['pmdec'] )

    sep = gaiacoord.separation(Pcoord)
    sep3d = gaiacoord.separation_3d(Pcoord)

    if verbose == True:
        print('Printing angular separations in degrees as sanity check')
        print(sep.degree)



    Pllbb     = bc.radec_to_lb(Pcoord.ra.value , Pcoord.dec.value , degree=True)
    Ppmllpmbb = bc.pmrapmdec_to_pmllpmbb( Pcoord.pm_ra_cosdec.value , Pcoord.pm_dec.value , \
                                         Pcoord.ra.value , Pcoord.dec.value , degree=True )
    Pvxvyvz   = bc.vrpmllpmbb_to_vxvyvz(Pcoord.radial_velocity.value , Ppmllpmbb[0] , Ppmllpmbb[1] , \
                                   Pllbb[0] , Pllbb[1] , Pcoord.distance.value/1000.0 , XYZ=False , degree=True)

    Cll = (math.atan2(Pvxvyvz[1],Pvxvyvz[0]) * 180.0/np.pi) % 360
    Cbb = math.atan2(Pvxvyvz[2],np.sqrt(Pvxvyvz[0]**2+Pvxvyvz[1]**2)) * 180.0/np.pi

    Cradec = bc.lb_to_radec(Cll,Cbb,degree=True,epoch=2000.0)
    Ccoord = SkyCoord( ra=Cradec[0]*u.deg , dec=Cradec[1]*u.deg , distance=999999.9 , frame='icrs' )
    print('Convergent point: ',Ccoord)

    Cangle = gaiacoord.separation(Ccoord)
    zz = np.where( (Cangle.degree > 90.0) )
    if (np.array(zz).size > 0): Cangle[zz] = (180.0-Cangle[zz].degree)*u.deg

    if verbose == True:
        print('Science Target Name: ',targname)
        print('Science Target RA/DEC: ',Pcoord.ra.value,Pcoord.dec.value)
        print('Science Target Galactic Coordinates: ',Pllbb)
        print('Science Target UVW: ',Pvxvyvz)
        print('Convergent point RA/DEC: ',Ccoord.ra.value,Ccoord.dec.value)
        print()

    Gllbb = bc.radec_to_lb(gaiacoord.ra.value , gaiacoord.dec.value , degree=True)
    Gxyz = bc.lbd_to_XYZ( Gllbb[:,0] , Gllbb[:,1] , gaiacoord.distance/1000.0 , degree=True)
    Gvrpmllpmbb = bc.vxvyvz_to_vrpmllpmbb( \
                    Pvxvyvz[0]*np.ones(len(Gxyz[:,0])) , Pvxvyvz[1]*np.ones(len(Gxyz[:,1])) , Pvxvyvz[2]*np.ones(len(Gxyz[:,2])) , \
                    Gxyz[:,0] , Gxyz[:,1] , Gxyz[:,2] , XYZ=True)
    Gpmrapmdec = bc.pmllpmbb_to_pmrapmdec( Gvrpmllpmbb[:,1] , Gvrpmllpmbb[:,2] , Gllbb[:,0] , Gllbb[:,1] , degree=True)

    # Code in case I want to do chi^2 cuts someday
    Gvtanerr = 1.0 * np.ones(len(Gxyz[:,0]))
    Gpmerr = Gvtanerr * 206265000.0 * 3.154e7 / (gaiacoord.distance.value * 3.086e13)


    Gchi2 = ( (Gpmrapmdec[:,0]-gaiacoord.pm_ra_cosdec.value)**2 + (Gpmrapmdec[:,1]-gaiacoord.pm_dec.value)**2 )**0.5
    Gchi2 = Gchi2 / Gpmerr
    if verbose == True:
        print('Predicted PMs if comoving:')
        print(Gpmrapmdec , "\n")
        print('Actual PMRAs from Gaia:')
        print(gaiacoord.pm_ra_cosdec.value , "\n")
        print('Actual PMDECs from Gaia:')
        print(gaiacoord.pm_dec.value , "\n")
        print('Predicted PM errors:')
        print(Gpmerr , "\n")
        print('Chi^2 values:')
        print(Gchi2)


    # Query external list(s) of RVs

    zz = np.where( (sep3d.value < searchradpc.value) & (Gchi2 < vlim.value) )
    yy = zz[0][np.argsort(sep3d[zz])]
    
    RV    = np.empty(np.array(r['ra']).size)
    RVerr = np.empty(np.array(r['ra']).size)
    RVsrc = np.array([ '                             None' for x in range(np.array(r['ra']).size) ])
    RV[:]    = np.nan
    RVerr[:] = np.nan


    print('Populating RV table')
    for x in range(0 , np.array(yy).size):
        if np.isnan(r['radial_velocity'][yy[x]]) == False:        # First copy over DR3 RVs
            RV[yy[x]] = r['radial_velocity'][yy[x]]
            if (np.ma.is_masked(r['teff_gspphot'][yy[x]]) == False):
                if (r['teff_gspphot'][yy[x]] >= 8500.0) & (np.ma.is_masked(r['grvs_mag'][yy[x]]) == False): 
                    RV[yy[x]] = r['radial_velocity'][yy[x]] - (7.98 - 1.135 * r['grvs_mag'][yy[x]])
                    if verbose == True:
                        print('Applying hot-star RV correction with GRVS: ',r['ra','teff_gspphot','grvs_mag','phot_rp_mean_mag','radial_velocity'][yy[x]]) 
                elif (r['teff_gspphot'][yy[x]] >= 8500.0) & (np.ma.is_masked(r['phot_rp_mean_mag'][yy[x]]) == False): 
                    RV[yy[x]] = r['radial_velocity'][yy[x]] - (7.98 - 1.135 * r['phot_rp_mean_mag'][yy[x]])
                    if verbose == True:
                        print('Applying hot-star RV correction with GRP: ',r['ra','teff_gspphot','grvs_mag','phot_rp_mean_mag','radial_velocity'][yy[x]]) 
            RVerr[yy[x]] = r['radial_velocity_error'][yy[x]]
            RVsrc[yy[x]] = 'Gaia DR3'
    if os.path.isfile('LocalRV.csv'):
        with open('LocalRV.csv') as csvfile:                          # Now check for a local RV that would supercede
            readCSV = csv.reader(csvfile, delimiter=',')
            next(readCSV)
            for row in readCSV:
                ww = np.where(r['DESIGNATION'] == row[0])[0]
                if ( (np.array(ww).size == 1) & (RVerr[ww] > float(row[3]))  ):
                    RV[ww]    = row[2]
                    RVerr[ww] = row[3]
                    RVsrc[ww] = row[4]
                    if verbose == True: 
                        print('Using stored RV: ',row)
                        print(r['ra','dec','phot_g_mean_mag'][ww])
                        print(RV[ww])
                        print(RVerr[ww])
                        print(RVsrc[ww])



    # Create Gaia CMD plot

    mamajek  = np.loadtxt(datapath+'/sptGBpRp.txt')
    pleiades = np.loadtxt(datapath+'/PleGBpRp.txt')
    tuchor   = np.loadtxt(datapath+'/TucGBpRp.txt')
    usco     = np.loadtxt(datapath+'/UScGBpRp.txt')
    chai     = np.loadtxt(datapath+'/ChaGBpRp.txt')

    zz = np.where( (sep3d.value < searchradpc.value) & (Gchi2 < vlim.value) & (np.isnan(r['bp_rp']) == False) ) # Note, this causes an error because NaNs
    yy = zz[0][np.argsort(sep3d[zz])]
    zz2= np.where( (sep3d.value < searchradpc.value) & (Gchi2 < vlim.value) & (sep.degree > 0.00001) & \
                 (r['phot_bp_rp_excess_factor'] < (1.3 + 0.06*r['bp_rp']**2)) & \
                 (Cangle.degree > convergcut) & \
                 (np.isnan(r['bp_rp']) == False) )                                                              # Note, this causes an error because NaNs
    yy2= zz2[0][np.argsort((-Gchi2)[zz2])]


    figname=outdir + targname.replace(" ", "") + "cmd.png"
    if verbose == True: print(figname)

    fig,ax1 = plt.subplots(figsize=(12,8))

    ax1.axis([ math.floor(min(r['bp_rp'][zz])) , \
               math.ceil(max(r['bp_rp'][zz])), \
               math.ceil(max((r['phot_g_mean_mag'][zz] - (5.0*np.log10(gaiacoord.distance[zz].value)-5.0))))+1, \
               math.floor(min((r['phot_g_mean_mag'][zz] - (5.0*np.log10(gaiacoord.distance[zz].value)-5.0))))-1 ] )
    ax1.set_xlabel(r'$B_p-R_p$ (mag)' , fontsize=16)
    ax1.set_ylabel(r'$M_G$ (mag)' , fontsize=16)
    ax1.tick_params(axis='both',which='major',labelsize=12)

    ax2 = ax1.twiny()
    ax2.set_xlim(ax1.get_xlim())
    spttickvals = np.array([ -0.037 , 0.377 , 0.782 , 0.980 , 1.84 , 2.50 , 3.36 , 4.75 ])
    sptticklabs = np.array([ 'A0' , 'F0' , 'G0' , 'K0' , 'M0' , 'M3' , 'M5' , 'M7' ])
    xx = np.where( (spttickvals >= math.floor(min(r['bp_rp'][zz]))) & (spttickvals <= math.ceil(max(r['bp_rp'][zz]))) )[0]
    ax2.set_xticks(spttickvals[xx])
    ax2.set_xticklabels( sptticklabs[xx] )
    ax2.set_xlabel('SpT' , fontsize=16, labelpad=15)
    ax2.tick_params(axis='both',which='major',labelsize=12)

    ax1.plot(    chai[:,1] ,     chai[:,0]  , zorder=1 , label='Cha-I (0-5 Myr)')
    ax1.plot(    usco[:,1] ,     usco[:,0]  , zorder=2 , label='USco (11 Myr)')
    ax1.plot(  tuchor[:,1] ,   tuchor[:,0]  , zorder=3 , label='Tuc-Hor (40 Myr)')
    ax1.plot(pleiades[:,1] , pleiades[:,0]  , zorder=4 , label='Pleiades (125 Myr)')
    ax1.plot( mamajek[:,2] ,  mamajek[:,1]  , zorder=5 , label='Mamajek MS')

    for x in range(0 , np.array(yy2).size):
        msize  = (17-12.0*(sep3d[yy2[x]].value/searchradpc.value))**2
        mcolor = Gchi2[yy2[x]]
        medge  = 'black'
        mzorder= 7
        if (r['ruwe'][yy2[x]] < 1.2):
            mshape='o'
        if (r['ruwe'][yy2[x]] >= 1.2):
            mshape='s'
        if (rvcut != None): 
            if (np.isnan(RV[yy2[x]])==False) & (np.abs(RV[yy2[x]]-Gvrpmllpmbb[yy2[x],0]) > rvcut) & (np.abs(RV[yy2[x]]-Gvrpmllpmbb[yy2[x],0])/RVerr[yy2[x]] > 2.0):
                mshape='+'
                mcolor='black'
                mzorder=6
            if (np.isnan(RV[yy2[x]])==False) & (np.abs(RV[yy2[x]]-Gvrpmllpmbb[yy2[x],0]) <= rvcut):
                medge='blue'
        if (mcolor == 'black'):
            ddd = ax1.scatter([ r['bp_rp'][yy2[x]] ] , [ (r['phot_g_mean_mag'][yy2[x]] - (5.0*np.log10(gaiacoord.distance[yy2[x]].value)-5.0)) ] , \
                s=msize , c=mcolor , marker=mshape , edgecolors=medge , zorder=mzorder , \
                vmin=0.0 , vmax=vlim.value , cmap='cubehelix' , label='_nolabel' )
        else:
            ccc = ax1.scatter([ r['bp_rp'][yy2[x]] ] , [ (r['phot_g_mean_mag'][yy2[x]] - (5.0*np.log10(gaiacoord.distance[yy2[x]].value)-5.0)) ] , \
                s=msize , c=mcolor , marker=mshape , edgecolors=medge , zorder=mzorder , \
                vmin=0.0 , vmax=vlim.value , cmap='cubehelix' , label='_nolabel' )

    temp1 = ax1.scatter([] , [] , c='white' , edgecolors='black', marker='o' , s=12**2 , label = 'RUWE < 1.2')
    temp2 = ax1.scatter([] , [] , c='white' , edgecolors='black', marker='s' , s=12**2 , label = 'RUWE >= 1.2')
    temp3 = ax1.scatter([] , [] , c='white' , edgecolors='blue' , marker='o' , s=12**2 , label = 'RV Comoving')
    temp4 = ax1.scatter([] , [] , c='black' , marker='+' , s=12**2 , label = 'RV Outlier')
    

    ax1.plot(r['bp_rp'][yy[0]] , (r['phot_g_mean_mag'][yy[0]] - (5.0*np.log10(gaiacoord.distance[yy[0]].value)-5.0)) , \
             'rx' , markersize=18 , mew=3 , markeredgecolor='red' , zorder=10 , label=targname)

    ax1.arrow( 1.3 , 2.5 , 0.374, 0.743 , length_includes_head=True , head_width=0.07 , head_length = 0.10 )
    ax1.text(  1.4 , 2.3, r'$A_V=1$' , fontsize=12)



    ax1.legend(fontsize=11)
    cb = plt.colorbar(ccc , ax=ax1)
    cb.set_label(label='Velocity Difference (km/s)',fontsize=14)
    plt.savefig(figname , bbox_inches='tight', pad_inches=0.2 , dpi=200)
    if showplots == True: plt.show()
    plt.close('all')


    # Create PM plot


    zz2= np.where( (sep3d.value < searchradpc.value) & (Gchi2 < vlim.value) & (sep.degree > 0.00001) & (Cangle.degree > convergcut) )
    yy2= zz2[0][np.argsort((-Gchi2)[zz2])]
    zz3= np.where( (sep3d.value < searchradpc.value) & (sep.degree > 0.00001) )

    figname=outdir + targname.replace(" ", "") + "pmd.png"

    fig,ax1 = plt.subplots(figsize=(12,8))

    ax1.axis([ (max(r['pmra'][zz2]) + 0.05*np.ptp(r['pmra'][zz2]) ) , \
           (min(r['pmra'][zz2]) - 0.05*np.ptp(r['pmra'][zz2]) ) , \
           (min(r['pmdec'][zz2])- 0.05*np.ptp(r['pmra'][zz2]) ) , \
           (max(r['pmdec'][zz2])+ 0.05*np.ptp(r['pmra'][zz2]) ) ] )
    ax1.tick_params(axis='both',which='major',labelsize=16)

    if  ((max(r['pmra'][zz2]) + 0.05*np.ptp(r['pmra'][zz2])) > 0.0) & \
            ((min(r['pmra'][zz2]) - 0.05*np.ptp(r['pmra'][zz2])) < 0.0) & \
            ((min(r['pmdec'][zz2])- 0.05*np.ptp(r['pmra'][zz2])) < 0.0) & \
            ((max(r['pmdec'][zz2])+ 0.05*np.ptp(r['pmra'][zz2])) > 0.0):
        ax1.plot( [0.0,0.0] , [-1000.0,1000.0] , 'k--' , linewidth=1 )
        ax1.plot( [-1000.0,1000.0] , [0.0,0.0] , 'k--' , linewidth=1 )

    ax1.errorbar( (r['pmra'][yy2]) , (r['pmdec'][yy2]) , \
            yerr=(r['pmdec_error'][yy2]) , xerr=(r['pmra_error'][yy2]) , fmt='none' , ecolor='k' )

    ax1.scatter( [ (r['pmra'][zz3]) ] , [ (r['pmdec'][zz3]) ] , \
              s=(0.5)**2 , marker='o' , c='black' , zorder=2 , label='Field' )

    for x in range(0 , np.array(yy2).size):
        msize  = (17-12.0*(sep3d[yy2[x]].value/searchradpc.value))**2
        mcolor = Gchi2[yy2[x]]
        medge  = 'black'
        mzorder= 7
        if (r['ruwe'][yy2[x]] < 1.2):
            mshape='o'
        if (r['ruwe'][yy2[x]] >= 1.2):
            mshape='s'
        if (rvcut != None): 
            if (np.isnan(RV[yy2[x]])==False) & (np.abs(RV[yy2[x]]-Gvrpmllpmbb[yy2[x],0]) > rvcut) & (np.abs(RV[yy2[x]]-Gvrpmllpmbb[yy2[x],0])/RVerr[yy2[x]] > 2.0):
                mshape='+'
                mcolor='black'
                mzorder=6
            if (np.isnan(RV[yy2[x]])==False) & (np.abs(RV[yy2[x]]-Gvrpmllpmbb[yy2[x],0]) <= rvcut):
                medge='blue'
        if (mcolor == 'black'):
            ddd = ax1.scatter([ r['pmra'][yy2[x]] ] , [ r['pmdec'][yy2[x]] ] , \
                s=msize , c=mcolor , marker=mshape , edgecolors=medge , zorder=mzorder , \
                vmin=0.0 , vmax=vlim.value , cmap='cubehelix' , label='_nolabel' )
        else:
            ccc = ax1.scatter([ r['pmra'][yy2[x]] ] , [ r['pmdec'][yy2[x]] ], \
                s=msize , c=mcolor , marker=mshape , edgecolors=medge , zorder=mzorder , \
                vmin=0.0 , vmax=vlim.value , cmap='cubehelix' , label='_nolabel' )


    temp1 = ax1.scatter([] , [] , c='white' , edgecolors='black', marker='o' , s=12**2 , label = 'RUWE < 1.2')
    temp2 = ax1.scatter([] , [] , c='white' , edgecolors='black', marker='s' , s=12**2 , label = 'RUWE >= 1.2')
    temp3 = ax1.scatter([] , [] , c='white' , edgecolors='blue' , marker='o' , s=12**2 , label = 'RV Comoving')
    temp4 = ax1.scatter([] , [] , c='black' , marker='+' , s=12**2 , label = 'RV Outlier')

    ax1.plot( Pgaia['pmra'][minpos] , Pgaia['pmdec'][minpos] , \
         'rx' , markersize=18 , mew=3 , markeredgecolor='red' , zorder=3 , label=targname)

    ax1.set_xlabel(r'$\mu_{RA}$ (mas/yr)' , fontsize=22 , labelpad=10)
    ax1.set_ylabel(r'$\mu_{DEC}$ (mas/yr)' , fontsize=22 , labelpad=10)
    ax1.legend(fontsize=12)

    cb = plt.colorbar(ccc , ax=ax1)
    cb.set_label(label='Tangential Velocity Difference (km/s)',fontsize=18 , labelpad=10)
    plt.savefig(figname , bbox_inches='tight', pad_inches=0.2 , dpi=200)
    if showplots == True: plt.show()
    plt.close('all')


    # Create RV plot

    zz2= np.where( (sep3d.value < searchradpc.value) & (Gchi2 < vlim.value) & (sep.degree > 0.00001) & \
             (Cangle.degree > convergcut) & \
             (np.isnan(RV) == False) )
    yy2= zz2[0][np.argsort((-Gchi2)[zz2])]

    zz3= np.where( (sep3d.value < searchradpc.value) & (Gchi2 < vlim.value) & (sep.degree > 0.00001) & \
             (np.isnan(RV) == False) & (np.isnan(r['phot_g_mean_mag']) == False) & \
             (Cangle.degree > convergcut) & \
             (np.abs(RV-Gvrpmllpmbb[:,0]) < 20.0) ) # Just to set Y axis

    fig,ax1 = plt.subplots(figsize=(12,8))
    ax1.axis([ -20.0 , +20.0, \
           max( np.append( np.array(r['phot_g_mean_mag'][zz3] - (5.0*np.log10(gaiacoord.distance[zz3].value)-5.0)) ,  0.0 )) + 0.3 , \
           min( np.append( np.array(r['phot_g_mean_mag'][zz3] - (5.0*np.log10(gaiacoord.distance[zz3].value)-5.0)) , 15.0 )) - 0.3   ])
    ax1.tick_params(axis='both',which='major',labelsize=16)

    ax1.plot( [0.0,0.0] , [-20.0,25.0] , 'k--' , linewidth=1 )

    ax1.errorbar( (RV[yy2]-Gvrpmllpmbb[yy2,0]) , \
           (r['phot_g_mean_mag'][yy2] - (5.0*np.log10(gaiacoord.distance[yy2].value)-5.0)) , \
            yerr=None,xerr=(RVerr[yy2]) , fmt='none' , ecolor='k' )

    nrvcut  = 0
    nrvpass = 0
    for x in range(0 , np.array(yy2).size):
        msize  = (17-12.0*(sep3d[yy2[x]].value/searchradpc.value))**2
        mcolor = Gchi2[yy2[x]]
        medge  = 'black'
        mzorder= 2
        if (r['ruwe'][yy2[x]] < 1.2):
            mshape='o'
        if (r['ruwe'][yy2[x]] >= 1.2):
            mshape='s'
        ccc = ax1.scatter( [ (RV[yy2[x]]-Gvrpmllpmbb[yy2[x],0]) ] , \
                [ (r['phot_g_mean_mag'][yy2[x]] - (5.0*np.log10(gaiacoord.distance[yy2[x]].value)-5.0)) ] , \
                s=msize , c=mcolor , marker=mshape , edgecolors=medge , zorder=mzorder , \
                vmin=0.0 , vmax=vlim.value , cmap='cubehelix' , label='_nolabel' )
        if (rvcut != None): 
            if (np.isnan(RV[yy2[x]])==False) & (np.abs(RV[yy2[x]]-Gvrpmllpmbb[yy2[x],0]) >  rvcut): nrvcut  = nrvcut  + 1
            if (np.isnan(RV[yy2[x]])==False) & (np.abs(RV[yy2[x]]-Gvrpmllpmbb[yy2[x],0]) <= rvcut): nrvpass = nrvpass + 1

    if (rvcut != None): print('Number with RV outside/inside selection range: ',nrvcut,' ',nrvpass)

    temp1 = ax1.scatter([] , [] , c='white' , edgecolors='black', marker='o' , s=12**2 , label = 'RUWE < 1.2')
    temp2 = ax1.scatter([] , [] , c='white' , edgecolors='black', marker='s' , s=12**2 , label = 'RUWE >= 1.2')
    temp3 = ax1.scatter([] , [] , c='white' , edgecolors='blue' , marker='o' , s=12**2 , label = 'RV Comoving')

    if ( (Pgaia['phot_g_mean_mag'][minpos] - (5.0*np.log10(Pcoord.distance.value)-5.0)) < \
                                     (max( np.append( np.array(r['phot_g_mean_mag'][zz3] - (5.0*np.log10(gaiacoord.distance[zz3].value)-5.0)) , 0.0 )) + 0.3) ):
        ax1.plot( [0.0] , (Pgaia['phot_g_mean_mag'][minpos] - (5.0*np.log10(Pcoord.distance.value)-5.0)) , \
                  'rx' , markersize=18 , mew=3 , markeredgecolor='red' , zorder=3 , label=targname)


    ax1.set_ylabel(r'$M_G$ (mag)' , fontsize=22 , labelpad=10)
    ax1.set_xlabel(r'$v_{r,obs}-v_{r,pred}$ (km/s)' , fontsize=22 , labelpad=10)
    ax1.legend(fontsize=12)

    cb = plt.colorbar(ccc , ax=ax1)
    cb.set_label(label='Tangential Velocity Difference (km/s)',fontsize=18 , labelpad=10)

    figname=outdir + targname.replace(" ", "") + "drv.png"
    plt.savefig(figname , bbox_inches='tight', pad_inches=0.2 , dpi=200)
    if showplots == True: plt.show()
    plt.close('all')



    
    # Create XYZ plot

    Pxyz = bc.lbd_to_XYZ( Pllbb[0] , Pllbb[1] , Pcoord.distance.value/1000.0 , degree=True)

    fig,axs = plt.subplots(2,2)
    fig.set_figheight(16)
    fig.set_figwidth(16)
    fig.subplots_adjust(hspace=0.03,wspace=0.03)

    zz2= np.where( (sep3d.value < searchradpc.value) & (Gchi2 < vlim.value) & (sep.degree > 0.00001) & (Cangle.degree > convergcut) )
    yy2= zz2[0][np.argsort((-Gchi2)[zz2])]

    for x in range(0 , np.array(yy2).size):
        msize  = (17-12.0*(sep3d[yy2[x]].value/searchradpc.value))**2
        mcolor = Gchi2[yy2[x]]
        medge  = 'black'
        mzorder= 3
        if (r['ruwe'][yy2[x]] < 1.2):
            mshape='o'
        if (r['ruwe'][yy2[x]] >= 1.2):
            mshape='s'
        if (rvcut != None): 
            if (np.isnan(RV[yy2[x]])==False) & (np.abs(RV[yy2[x]]-Gvrpmllpmbb[yy2[x],0]) > rvcut) & (np.abs(RV[yy2[x]]-Gvrpmllpmbb[yy2[x],0])/RVerr[yy2[x]] > 2.0):
                mshape='+'
                mcolor='black'
                mzorder=2
            if (np.isnan(RV[yy2[x]])==False) & (np.abs(RV[yy2[x]]-Gvrpmllpmbb[yy2[x],0]) <= rvcut):
                medge='blue'
        ccc = axs[0,0].scatter( [ 1000.0*Gxyz[yy2[x],0] ] , [ 1000.0*Gxyz[yy2[x],1] ] , \
                s=msize , c=mcolor , marker=mshape , edgecolors=medge , zorder=mzorder , \
                vmin=0.0 , vmax=vlim.value , cmap='cubehelix' , label='_nolabel' )
        ccc = axs[0,1].scatter( [ 1000.0*Gxyz[yy2[x],2] ] , [ 1000.0*Gxyz[yy2[x],1] ] , \
                s=msize , c=mcolor , marker=mshape , edgecolors=medge , zorder=mzorder , \
                vmin=0.0 , vmax=vlim.value , cmap='cubehelix' , label='_nolabel' )
        if (mcolor == 'black'):
            ccc = axs[1,0].scatter( [ 1000.0*Gxyz[yy2[x],0] ] , [ 1000.0*Gxyz[yy2[x],2] ] , \
                s=msize , c=mcolor , marker=mshape , edgecolors=medge , zorder=mzorder , \
                vmin=0.0 , vmax=vlim.value , cmap='cubehelix' , label='_nolabel' )
        else:
            ddd = axs[1,0].scatter( [ 1000.0*Gxyz[yy2[x],0] ] , [ 1000.0*Gxyz[yy2[x],2] ] , \
                s=msize , c=mcolor , marker=mshape , edgecolors=medge , zorder=mzorder , \
                vmin=0.0 , vmax=vlim.value , cmap='cubehelix' , label='_nolabel' )

    temp1 = axs[0,0].scatter([] , [] , c='white' , edgecolors='black', marker='o' , s=12**2 , label = 'RUWE < 1.2')
    temp2 = axs[0,0].scatter([] , [] , c='white' , edgecolors='black', marker='s' , s=12**2 , label = 'RUWE >= 1.2')
    temp3 = axs[0,0].scatter([] , [] , c='white' , edgecolors='blue' , marker='o' , s=12**2 , label = 'RV Comoving')
    temp4 = axs[0,0].scatter([] , [] , c='black' , marker='+' , s=12**2 , label = 'RV Outlier')

    axs[0,0].plot( 1000.0*Pxyz[0] , 1000.0*Pxyz[1] , 'rx' , markersize=18 , mew=3 , markeredgecolor='red')
    axs[0,1].plot( 1000.0*Pxyz[2] , 1000.0*Pxyz[1] , 'rx' , markersize=18 , mew=3 , markeredgecolor='red')
    axs[1,0].plot( 1000.0*Pxyz[0] , 1000.0*Pxyz[2] , 'rx' , markersize=18 , mew=3 , markeredgecolor='red' , zorder=1 , label = targname)

    axs[0,0].set_xlim( [1000.0*Pxyz[0]-(search_radius+1.0) , 1000.0*Pxyz[0]+(search_radius+1.0)] )
    axs[0,0].set_ylim( [1000.0*Pxyz[1]-(search_radius+1.0) , 1000.0*Pxyz[1]+(search_radius+1.0)] )
    axs[0,1].set_xlim( [1000.0*Pxyz[2]-(search_radius+1.0) , 1000.0*Pxyz[2]+(search_radius+1.0)] )
    axs[0,1].set_ylim( [1000.0*Pxyz[1]-(search_radius+1.0) , 1000.0*Pxyz[1]+(search_radius+1.0)] )
    axs[1,0].set_xlim( [1000.0*Pxyz[0]-(search_radius+1.0) , 1000.0*Pxyz[0]+(search_radius+1.0)] )
    axs[1,0].set_ylim( [1000.0*Pxyz[2]-(search_radius+1.0) , 1000.0*Pxyz[2]+(search_radius+1.0)] )
    
    axs[0,0].set_xlabel(r'$X$ (pc)',fontsize=20,labelpad=10)
    axs[0,0].set_ylabel(r'$Y$ (pc)',fontsize=20,labelpad=10)

    axs[1,0].set_xlabel(r'$X$ (pc)',fontsize=20,labelpad=10)
    axs[1,0].set_ylabel(r'$Z$ (pc)',fontsize=20,labelpad=10)

    axs[0,1].set_xlabel(r'$Z$ (pc)',fontsize=20,labelpad=10)
    axs[0,1].set_ylabel(r'$Y$ (pc)',fontsize=20,labelpad=10)

    axs[0,0].xaxis.set_ticks_position('top')
    axs[0,1].xaxis.set_ticks_position('top')
    axs[0,1].yaxis.set_ticks_position('right')

    axs[0,0].xaxis.set_label_position('top')
    axs[0,1].xaxis.set_label_position('top')
    axs[0,1].yaxis.set_label_position('right')

    for aa in [0,1]:
        for bb in [0,1]:
            axs[aa,bb].tick_params(top=True,bottom=True,left=True,right=True,direction='in',labelsize=18)

    fig.delaxes(axs[1][1])
    strsize = 26
    if (len(targname) > 12.0): strsize = np.floor(24 / (len(targname)/14.5))
    fig.legend( bbox_to_anchor=(0.92,0.37) , prop={'size':strsize})

    cbaxes = fig.add_axes([0.55,0.14,0.02,0.34])
    cb = plt.colorbar( ddd , cax=cbaxes )
    cb.set_label( label='Velocity Difference (km/s)' , fontsize=24 , labelpad=20 )
    cb.ax.tick_params(labelsize=18)

    figname=outdir + targname.replace(" ", "") + "xyz.png"
    plt.savefig(figname , bbox_inches='tight', pad_inches=0.2 , dpi=200)

    if showplots == True: plt.show()
    plt.close('all')


    # Create XYZ plot for only RV confirmed/disproven

    Pxyz = bc.lbd_to_XYZ( Pllbb[0] , Pllbb[1] , Pcoord.distance.value/1000.0 , degree=True)

    fig,axs = plt.subplots(2,2)
    fig.set_figheight(16)
    fig.set_figwidth(16)
    fig.subplots_adjust(hspace=0.03,wspace=0.03)

    Gchi3 = np.sqrt( Gchi2**2 + (RV-Gvrpmllpmbb[:,0])**2 )
    vlim3 = vlim.value
    if (rvcut != None): vlim3 = np.sqrt((vlim.value)**2 + rvcut**2)

    zz2= np.where( (sep3d.value < searchradpc.value) & (Gchi2 < vlim.value) & (sep.degree > 0.00001) & (Cangle.degree > convergcut) & \
        (np.isnan(RV)==False) )
    yy2= zz2[0][np.argsort((-Gchi3)[zz2])]

    for x in range(0 , np.array(yy2).size):
        msize  = (17-12.0*(sep3d[yy2[x]].value/searchradpc.value))**2
        mcolor = Gchi3[yy2[x]]
        medge  = 'black'
        mzorder= 3
        if (r['ruwe'][yy2[x]] < 1.2):
            mshape='o'
        if (r['ruwe'][yy2[x]] >= 1.2):
            mshape='s'
        if (rvcut != None):
            if (np.isnan(RV[yy2[x]])==False) & (np.abs(RV[yy2[x]]-Gvrpmllpmbb[yy2[x],0]) > rvcut) & (np.abs(RV[yy2[x]]-Gvrpmllpmbb[yy2[x],0])/RVerr[yy2[x]] > 2.0):
                mshape='+'
                mcolor='black'
                mzorder=2
            if (np.isnan(RV[yy2[x]])==False) & (np.abs(RV[yy2[x]]-Gvrpmllpmbb[yy2[x],0]) <= rvcut):
                medge='blue'
        ccc = axs[0,0].scatter( [ 1000.0*Gxyz[yy2[x],0] ] , [ 1000.0*Gxyz[yy2[x],1] ] , \
                s=msize , c=mcolor , marker=mshape , edgecolors=medge , zorder=mzorder , \
                vmin=0.0 , vmax=vlim3 , cmap='cubehelix' , label='_nolabel' )
        ccc = axs[0,1].scatter( [ 1000.0*Gxyz[yy2[x],2] ] , [ 1000.0*Gxyz[yy2[x],1] ] , \
                s=msize , c=mcolor , marker=mshape , edgecolors=medge , zorder=mzorder , \
                vmin=0.0 , vmax=vlim3 , cmap='cubehelix' , label='_nolabel' )
        if (mcolor == 'black'):
            ccc = axs[1,0].scatter( [ 1000.0*Gxyz[yy2[x],0] ] , [ 1000.0*Gxyz[yy2[x],2] ] , \
                s=msize , c=mcolor , marker=mshape , edgecolors=medge , zorder=mzorder , \
                vmin=0.0 , vmax=vlim3 , cmap='cubehelix' , label='_nolabel' )
        else:
            ddd = axs[1,0].scatter( [ 1000.0*Gxyz[yy2[x],0] ] , [ 1000.0*Gxyz[yy2[x],2] ] , \
                s=msize , c=mcolor , marker=mshape , edgecolors=medge , zorder=mzorder , \
                vmin=0.0 , vmax=vlim3 , cmap='cubehelix' , label='_nolabel' )

    temp1 = axs[0,0].scatter([] , [] , c='white' , edgecolors='black', marker='o' , s=12**2 , label = 'RUWE < 1.2')
    temp2 = axs[0,0].scatter([] , [] , c='white' , edgecolors='black', marker='s' , s=12**2 , label = 'RUWE >= 1.2')
    temp3 = axs[0,0].scatter([] , [] , c='white' , edgecolors='blue' , marker='o' , s=12**2 , label = 'RV Comoving')
    temp4 = axs[0,0].scatter([] , [] , c='black' , marker='+' , s=12**2 , label = 'RV Outlier')

    axs[0,0].plot( 1000.0*Pxyz[0] , 1000.0*Pxyz[1] , 'rx' , markersize=18 , mew=3 , markeredgecolor='red')
    axs[0,1].plot( 1000.0*Pxyz[2] , 1000.0*Pxyz[1] , 'rx' , markersize=18 , mew=3 , markeredgecolor='red')
    axs[1,0].plot( 1000.0*Pxyz[0] , 1000.0*Pxyz[2] , 'rx' , markersize=18 , mew=3 , markeredgecolor='red' , zorder=1 , label = targname)

    axs[0,0].set_xlim( [1000.0*Pxyz[0]-(search_radius+1.0) , 1000.0*Pxyz[0]+(search_radius+1.0)] )
    axs[0,0].set_ylim( [1000.0*Pxyz[1]-(search_radius+1.0) , 1000.0*Pxyz[1]+(search_radius+1.0)] )
    axs[0,1].set_xlim( [1000.0*Pxyz[2]-(search_radius+1.0) , 1000.0*Pxyz[2]+(search_radius+1.0)] )
    axs[0,1].set_ylim( [1000.0*Pxyz[1]-(search_radius+1.0) , 1000.0*Pxyz[1]+(search_radius+1.0)] )
    axs[1,0].set_xlim( [1000.0*Pxyz[0]-(search_radius+1.0) , 1000.0*Pxyz[0]+(search_radius+1.0)] )
    axs[1,0].set_ylim( [1000.0*Pxyz[2]-(search_radius+1.0) , 1000.0*Pxyz[2]+(search_radius+1.0)] )
    
    axs[0,0].set_xlabel(r'$X$ (pc)',fontsize=20,labelpad=10)
    axs[0,0].set_ylabel(r'$Y$ (pc)',fontsize=20,labelpad=10)

    axs[1,0].set_xlabel(r'$X$ (pc)',fontsize=20,labelpad=10)
    axs[1,0].set_ylabel(r'$Z$ (pc)',fontsize=20,labelpad=10)

    axs[0,1].set_xlabel(r'$Z$ (pc)',fontsize=20,labelpad=10)
    axs[0,1].set_ylabel(r'$Y$ (pc)',fontsize=20,labelpad=10)

    axs[0,0].xaxis.set_ticks_position('top')
    axs[0,1].xaxis.set_ticks_position('top')
    axs[0,1].yaxis.set_ticks_position('right')

    axs[0,0].xaxis.set_label_position('top')
    axs[0,1].xaxis.set_label_position('top')
    axs[0,1].yaxis.set_label_position('right')

    for aa in [0,1]:
        for bb in [0,1]:
            axs[aa,bb].tick_params(top=True,bottom=True,left=True,right=True,direction='in',labelsize=18)

    fig.delaxes(axs[1][1])
    strsize = 26
    if (len(targname) > 12.0): strsize = np.floor(24 / (len(targname)/14.5))
    fig.legend( bbox_to_anchor=(0.92,0.37) , prop={'size':strsize})

    cbaxes = fig.add_axes([0.55,0.14,0.02,0.34])
    cb = plt.colorbar( ddd , cax=cbaxes )
    cb.set_label( label='Velocity Difference (km/s)' , fontsize=24 , labelpad=20 )
    cb.ax.tick_params(labelsize=18)

    figname=outdir + targname.replace(" ", "") + "xyzRV.png"
    plt.savefig(figname , bbox_inches='tight', pad_inches=0.2 , dpi=200)

    if showplots == True: plt.show()
    plt.close('all')





    # Create sky map
    # Hacked from cartopy.mpl.gridliner
    _DEGREE_SYMBOL = u'\u00B0'
    def _east_west_formatted(longitude, num_format='g'):
        fmt_string = u'{longitude:{num_format}}{degree}'
        return fmt_string.format(longitude=(longitude if (longitude >= 0) else (longitude + 360)) , \
                                            num_format=num_format,degree=_DEGREE_SYMBOL)
    def _north_south_formatted(latitude, num_format='g'):
        fmt_string = u'{latitude:{num_format}}{degree}'
        return fmt_string.format(latitude=latitude, num_format=num_format,degree=_DEGREE_SYMBOL)
    LONGITUDE_FORMATTER = mticker.FuncFormatter(lambda v, pos:
                                                _east_west_formatted(v))
    LATITUDE_FORMATTER = mticker.FuncFormatter(lambda v, pos:
                                               _north_south_formatted(v))

    zz = np.where( (sep3d.value < searchradpc.value) & (Gchi2 < vlim.value) & (sep.degree > 0.00001) & (Cangle.degree > convergcut) )
    yy = zz[0][np.argsort((-Gchi2)[zz])]

    searchcircle = Pcoord.directional_offset_by( (np.arange(0,360)*u.degree) , searchraddeg*np.ones(360))
    circleRA = searchcircle.ra.value
    circleDE = searchcircle.dec.value
    ww = np.where(circleRA > 180.0)
    circleRA[ww] = circleRA[ww] - 360.0

    RAlist = gaiacoord.ra[yy].value
    DElist = gaiacoord.dec[yy].value
    ww = np.where( RAlist > 180.0 )
    RAlist[ww] = RAlist[ww] - 360.0

    polelat = ((Pcoord.dec.value+90) if (Pcoord.dec.value<0) else (90-Pcoord.dec.value))
    polelong= (Pcoord.ra.value if (Pcoord.dec.value<0.0) else (Pcoord.ra.value+180.0))
    polelong= (polelong if polelong < 180 else polelong - 360.0)

    if verbose == True:
        print('Alignment variables: ',polelat,polelong,Pcoord.ra.value)
        print(Pcoord.dec.value+searchraddeg.value)
    rotated_pole = ccrs.RotatedPole( \
        pole_latitude=polelat , \
        pole_longitude=polelong , \
        central_rotated_longitude=90.0 )#\
    #    (Pcoord.ra.value if (Pcoord.dec.value > 0.0) else (Pcoord.ra.value+180.0)) )

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(1, 1, 1, projection=rotated_pole)

    ax.gridlines(draw_labels=True,x_inline=True,y_inline=True, \
                 xformatter=LONGITUDE_FORMATTER,yformatter=LATITUDE_FORMATTER)
    ax.plot( circleRA , circleDE , c="gray" , ls="--" , transform=ccrs.Geodetic())
    
    figname=outdir + targname.replace(" ", "") + "sky.png"

    base=plt.cm.get_cmap('cubehelix')

    for x in range(0 , np.array(yy).size):
        msize  = (17-12.0*(sep3d[yy[x]].value/searchradpc.value))
        mcolor = base(Gchi2[yy[x]]/vlim.value)
        medge  = 'black'
        mzorder= 3
        if (r['ruwe'][yy[x]] < 1.2):
            mshape='o'
        if (r['ruwe'][yy[x]] >= 1.2):
            mshape='s'
        if (rvcut != None):
            if (np.isnan(RV[yy[x]])==False) & (np.abs(RV[yy[x]]-Gvrpmllpmbb[yy[x],0]) > rvcut) & (np.abs(RV[yy[x]]-Gvrpmllpmbb[yy[x],0])/RVerr[yy[x]] > 2.0):
                mshape='+'
                mcolor='black'
                mzorder=2
            if (np.isnan(RV[yy[x]])==False) & (np.abs(RV[yy[x]]-Gvrpmllpmbb[yy[x],0]) <= rvcut):
                medge='blue'
        if (mcolor == 'black'):
            ddd = ax.plot( RAlist[x] , DElist[x] , marker=mshape ,  \
                markeredgecolor=medge , ms = msize , mfc = mcolor , transform=ccrs.Geodetic() )
        else:
            ccc = ax.plot( RAlist[x] , DElist[x] , marker=mshape ,  \
                markeredgecolor=medge , ms = msize , mfc = mcolor , transform=ccrs.Geodetic() )
        
    ax.plot( (Pcoord.ra.value-360.0) , Pcoord.dec.value , \
            'rx' , markersize=18 , mew=3 , transform=ccrs.Geodetic())

    plt.savefig(figname , bbox_inches='tight', pad_inches=0.2 , dpi=200)
    
    if showplots == True: plt.show()
    plt.close('all')

    ## Query GALEX and 2MASS data

    zz = np.where( (sep3d.value < searchradpc.value) & (Gchi2 < vlim.value) )
    yy = zz[0][np.argsort((-Gchi2)[zz])]
    
    NUVmag = np.empty(np.array(r['ra']).size)
    NUVerr = np.empty(np.array(r['ra']).size)
    NUVmag[:] = np.nan
    NUVerr[:] = np.nan

    print('Searching on neighbors in GALEX')
    ##suppress the stupid noresultswarning from the catalogs package
    warnings.filterwarnings("ignore",category=NoResultsWarning)

    for x in range(0 , np.array(yy).size):
        querystring=((str(gaiacoord.ra[yy[x]].value) if (gaiacoord.ra[yy[x]].value > 0) \
                      else str(gaiacoord.ra[yy[x]].value+360.0)) + " " + str(gaiacoord.dec[yy[x]].value))
        print('GALEX query ',x,' of ',np.array(yy).size, end='\r')
        if verbose == True: print('GALEX query ',x,' of ',np.array(yy).size)
        if verbose == True: print(querystring)
        if (DoGALEX == True): 
            galex = Catalogs.query_object(querystring , catalog="Galex" , radius=0.0028 , TIMEOUT=600)
            if ((np.where(galex['nuv_magerr'] > 0.0)[0]).size > 0):
                ww = np.where( (galex['nuv_magerr'] == min(galex['nuv_magerr'][np.where(galex['nuv_magerr'] > 0.0)])))
                NUVmag[yy[x]] = galex['nuv_mag'][ww][0]
                NUVerr[yy[x]] = galex['nuv_magerr'][ww][0]
                if verbose == True: print(galex['distance_arcmin','ra','nuv_mag','nuv_magerr'][ww])

        
    Jmag = np.empty(np.array(r['ra']).size)
    Jerr = np.empty(np.array(r['ra']).size)
    Jmag[:] = np.nan
    Jerr[:] = np.nan

    print('Searching on neighbors in 2MASS')

    for x in range(0 , np.array(yy).size):
        if ( np.isnan(NUVmag[yy[x]]) == False ):
            querycoord = SkyCoord((str(gaiacoord.ra[yy[x]].value) if (gaiacoord.ra[yy[x]].value > 0) else \
                     str(gaiacoord.ra[yy[x]].value+360.0)) , str(gaiacoord.dec[yy[x]].value) , \
                     unit=(u.deg,u.deg) , frame='icrs')
            print('2MASS query ',x,' of ',np.array(yy).size, end='\r')
            if verbose == True: print('2MASS query ',x,' of ',np.array(yy).size)
            if verbose == True: print(querycoord)
            tmass = []
            if (DoGALEX == True): 
                tmass = Irsa.query_region(querycoord , catalog='fp_psc' , radius='0d0m10s' )
                if ((np.where(tmass['j_m'] > -10.0)[0]).size > 0):
                    ww = np.where( (tmass['j_m'] == min(tmass['j_m'][np.where(tmass['j_m'] > 0.0)])))
                    Jmag[yy[x]] = tmass['j_m'][ww][0]
                    Jerr[yy[x]] = tmass['j_cmsig'][ww][0]
                    if verbose == True: print(tmass['j_m','j_cmsig'][ww])
        


    # Create GALEX plots
    mamajek = np.loadtxt(datapath+'/sptGBpRp.txt')
    f = interp1d( mamajek[:,2] , mamajek[:,0] , kind='cubic')

    zz2 = np.where( (sep3d.value < searchradpc.value) & (Gchi2 < vlim.value) )
    yy2 = zz[0][np.argsort(sep3d[zz])]
    zz = np.where( (sep3d.value < searchradpc.value) & (Gchi2 < vlim.value) & (sep.degree > 0.00001) & (Cangle.degree > convergcut) )
    yy = zz[0][np.argsort((-Gchi2)[zz])]

    fnuvj = (3631.0 * 10**6 * 10**(-0.4 * NUVmag)) / (1594.0 * 10**6 * 10**(-0.4 * Jmag))
    spt = f(r['bp_rp'].filled(np.nan))
    sptstring = ["nan" for x in range(np.array(r['bp_rp']).size)]
    for x in range(0 , np.array(zz2).size):
        if (round(spt[yy2[x]],1) >= 17.0) and (round(spt[yy2[x]],1) < 27.0):
            sptstring[yy2[x]] = 'M' + ('% 3.1f' % (round(spt[yy2[x]],1)-17.0)).strip()
        if (round(spt[yy2[x]],1) >= 16.0) and (round(spt[yy2[x]],1) < 17.0):
            sptstring[yy2[x]] = 'K' + ('% 3.1f' % (round(spt[yy2[x]],1)-9.0)).strip()
        if (round(spt[yy2[x]],1) >= 10.0) and (round(spt[yy2[x]],1) < 16.0):
            sptstring[yy2[x]] = 'K' + ('% 3.1f' % (round(spt[yy2[x]],1)-10.0)).strip()
        if (round(spt[yy2[x]],1) >= 0.0) and (round(spt[yy2[x]],1) < 10.0):
            sptstring[yy2[x]] = 'G' + ('% 3.1f' % (round(spt[yy2[x]],1)-0.0)).strip()
        if (round(spt[yy2[x]],1) >= -10.0) and (round(spt[yy2[x]],1) < 0.0):
            sptstring[yy2[x]] = 'F' + ('% 3.1f' % (round(spt[yy2[x]],1)+10.0)).strip()
        if (round(spt[yy2[x]],1) >= -20.0) and (round(spt[yy2[x]],1) < -10.0):
            sptstring[yy2[x]] = 'A' + ('% 3.1f' % (round(spt[yy2[x]],1)+20.0)).strip()       
        if (round(spt[yy2[x]],1) >= -30.0) and (round(spt[yy2[x]],1) < -20.0):
            sptstring[yy2[x]] = 'B' + ('% 3.1f' % (round(spt[yy2[x]],1)+30.0)).strip()  
    


    figname=outdir + targname.replace(" ", "") + "galex.png"
    if verbose == True: print(figname)
    ##Muck with the axis to get two x axes

    fig,ax1 = plt.subplots(figsize=(12,8))
    ax1.set_yscale('log')
    ax1.axis([5.0 , 24.0 , 0.000004 , 0.02])
    ax2 = ax1.twiny()
    ax2.set_xlim(ax1.get_xlim())
    ax1.set_xticks(np.array([5.0 , 10.0 , 15.0 , 17.0 , 22.0 , 24.0]))
    ax1.set_xticklabels(['G5','K0','K5','M0','M5','M7'])
    ax1.set_xlabel('SpT' , fontsize=20, labelpad=15)
    ax1.tick_params(axis='both',which='major',labelsize=16)
    ax2.set_xticks(np.array([5.0 , 10.0 , 15.0 , 17.0 , 22.0 , 24.0]))
    ax2.set_xticklabels(['0.85','0.98','1.45','1.84','3.36','4.75'])
    ax2.set_xlabel(r'$B_p-R_p$ (mag)' , fontsize=20, labelpad=15)
    ax2.tick_params(axis='both',which='major',labelsize=16)
    ax1.set_ylabel(r'$F_{NUV}/F_{J}$' , fontsize=22, labelpad=0)

    ##Hyades
    hyades = readsav(datapath +'/HYsaved.sav')
    hyadesfnuvj = (3631.0 * 10**6 * 10**(-0.4 * hyades['clnuv'])) / (1594.0 * 10**6 * 10**(-0.4 * hyades['clJ']))
    ax1.plot(hyades['clspt'] , hyadesfnuvj , 'x' , markersize=4 , mew=1 , markeredgecolor='black' , zorder=1 , label='Hyades' )

    for x in range(0 , np.array(yy).size):
        msize  = (17-12.0*(sep3d[yy[x]].value/searchradpc.value))**2
        mcolor = Gchi2[yy[x]]
        medge  = 'black'
        mzorder= 3
        if (r['ruwe'][yy[x]] < 1.2):
            mshape='o'
        if (r['ruwe'][yy[x]] >= 1.2):
            mshape='s'
        if (rvcut != None):
            if (np.isnan(RV[yy[x]])==False) & (np.abs(RV[yy[x]]-Gvrpmllpmbb[yy[x],0]) > rvcut) & (np.abs(RV[yy[x]]-Gvrpmllpmbb[yy[x],0])/RVerr[yy[x]] > 2.0):
                mshape='+'
                mcolor='black'
                mzorder=2
            if (np.isnan(RV[yy[x]])==False) & (np.abs(RV[yy[x]]-Gvrpmllpmbb[yy[x],0]) <= rvcut):
                medge='blue'
        if (mcolor == 'black'):
            ddd = ax1.scatter( [ spt[yy[x]] ] , [ fnuvj[yy[x]] ] , \
                s=msize , c=mcolor , marker=mshape , edgecolors=medge , zorder=mzorder , \
                vmin=0.0 , vmax=vlim.value , cmap='cubehelix' , label='_nolabel' )
        else:
            ccc = ax1.scatter( [ spt[yy[x]] ] , [ fnuvj[yy[x]] ] , \
                s=msize , c=mcolor , marker=mshape , edgecolors=medge , zorder=mzorder , \
                vmin=0.0 , vmax=vlim.value , cmap='cubehelix' , label='_nolabel' )

    temp1 = ax1.scatter([] , [] , c='white' , edgecolors='black', marker='o' , s=12**2 , label = 'RUWE < 1.2')
    temp2 = ax1.scatter([] , [] , c='white' , edgecolors='black', marker='s' , s=12**2 , label = 'RUWE >= 1.2')
    temp3 = ax1.scatter([] , [] , c='white' , edgecolors='blue' , marker='o' , s=12**2 , label = 'RV Comoving')
    temp4 = ax1.scatter([] , [] , c='black' , marker='+' , s=12**2 , label = 'RV Outlier')



    # Plot science target
    if (spt[yy2[0]] > 5): ax1.plot(spt[yy2[0]] , fnuvj[yy2[0]] , 'rx' , markersize=18 , mew=3 , markeredgecolor='red' , zorder=3 , label=targname )

    ax1.legend(fontsize=16 , loc='lower left')
    cb = fig.colorbar(ccc , ax=ax1)
    cb.set_label(label='Tangential Velocity Offset (km/s)',fontsize=13)
    if (DoGALEX == True): plt.savefig(figname , bbox_inches='tight', pad_inches=0.2 , dpi=200)
    if showplots == True: plt.show()
    plt.close('all')
    
    
    # Query CatWISE for W1+W2 and AllWISE for W3+W4

    zz = np.where( (sep3d.value < searchradpc.value) & (Gchi2 < vlim.value) )
    yy = zz[0][np.argsort((-Gchi2)[zz])]

    WISEmag = np.empty([np.array(r['ra']).size,4])
    WISEerr = np.empty([np.array(r['ra']).size,4])
    WISEmag[:] = np.nan
    WISEerr[:] = np.nan

    print('Searching on neighbors in WISE')
    ##there's an annoying nan warning here, hide it for now as it's not a problem
    warnings.filterwarnings("ignore",category=UserWarning)

    for x in range(0 , np.array(yy).size):
        querycoord = SkyCoord((str(gaiacoord.ra[yy[x]].value) if (gaiacoord.ra[yy[x]].value > 0) else \
                     str(gaiacoord.ra[yy[x]].value+360.0)) , str(gaiacoord.dec[yy[x]].value) , \
                     unit=(u.deg,u.deg) , frame='icrs')
        print('WISE query ',x,' of ',np.array(yy).size, end='\r')
        if verbose == True: print('WISE query ',x,' of ',np.array(yy).size)
        if verbose == True: print(querycoord)
    
        wisecat = []
        if (DoWISE == True): 
            wisecat = Irsa.query_region(querycoord,catalog='catwise_2020' , radius='0d0m10s')
            if ((np.where(wisecat['w1mpro'] > -10.0)[0]).size > 0):
                ww = np.where( (wisecat['w1mpro'] == min( wisecat['w1mpro'][np.where(wisecat['w1mpro'] > -10.0)]) ))
                WISEmag[yy[x],0] = wisecat['w1mpro'][ww][0]
                WISEerr[yy[x],0] = wisecat['w1sigmpro'][ww][0]
            if ((np.where(wisecat['w2mpro'] > -10.0)[0]).size > 0):
                ww = np.where( (wisecat['w2mpro'] == min( wisecat['w2mpro'][np.where(wisecat['w2mpro'] > -10.0)]) ))
                WISEmag[yy[x],1] = wisecat['w2mpro'][ww][0]
                WISEerr[yy[x],1] = wisecat['w2sigmpro'][ww][0]
 
        if (DoWISE == True): 
            wisecat = Irsa.query_region(querycoord,catalog='allwise_p3as_psd' , radius='0d0m10s')
            if ((np.where(wisecat['w1mpro'] > -10.0)[0]).size > 0):
                ww = np.where( (wisecat['w1mpro'] == min( wisecat['w1mpro'][np.where(wisecat['w1mpro'] > -10.0)]) ))
                if (np.isnan(WISEmag[yy[x],0]) == True) | (wisecat['w1mpro'][ww][0] < 11.0):				# Note, only if CatWISE absent/saturated
                    WISEmag[yy[x],0] = wisecat['w1mpro'][ww][0]
                    WISEerr[yy[x],0] = wisecat['w1sigmpro'][ww][0]
            if ((np.where(wisecat['w2mpro'] > -10.0)[0]).size > 0):
                ww = np.where( (wisecat['w2mpro'] == min( wisecat['w2mpro'][np.where(wisecat['w2mpro'] > -10.0)]) ))
                if (np.isnan(WISEmag[yy[x],1]) == True) | (wisecat['w2mpro'][ww][0] < 11.0):				# Note, only if CatWISE absent/saturated
                    WISEmag[yy[x],1] = wisecat['w2mpro'][ww][0]
                    WISEerr[yy[x],1] = wisecat['w2sigmpro'][ww][0]
            if ((np.where(wisecat['w3mpro'] > -10.0)[0]).size > 0):
                ww = np.where( (wisecat['w3mpro'] == min( wisecat['w3mpro'][np.where(wisecat['w3mpro'] > -10.0)]) ))
                WISEmag[yy[x],2] = wisecat['w3mpro'][ww][0]
                WISEerr[yy[x],2] = wisecat['w3sigmpro'][ww][0]
            if ((np.where(wisecat['w4mpro'] > -10.0)[0]).size > 0):
                ww = np.where( (wisecat['w4mpro'] == min( wisecat['w4mpro'][np.where(wisecat['w4mpro'] > -10.0)]) ))
                WISEmag[yy[x],3] = wisecat['w4mpro'][ww][0]
                WISEerr[yy[x],3] = wisecat['w4sigmpro'][ww][0]
        
        if verbose == True: print(yy[x],WISEmag[yy[x],:],WISEerr[yy[x],:])

    # Create WISE plots

    W13 = WISEmag[:,0]-WISEmag[:,2]
    W13err = ( WISEerr[:,0]**2 + WISEerr[:,2]**2 )**0.5

    zz = np.argwhere( np.isnan(W13err) )
    W13[zz] = np.nan
    W13err[zz] = np.nan

    zz = np.where( (W13err > 0.15) )
    W13[zz] = np.nan
    W13err[zz] = np.nan
    warnings.filterwarnings("default",category=UserWarning)




    zz2 = np.where( (sep3d.value < searchradpc.value) & (Gchi2 < vlim.value))
    yy2 = zz2[0][np.argsort(sep3d[zz2])]
    zz = np.where( (sep3d.value < searchradpc.value) & (Gchi2 < vlim.value) & (sep.degree > 0.00001) & (Cangle.degree > convergcut) )
    yy = zz[0][np.argsort((-Gchi2)[zz])]

    figname=outdir + targname.replace(" ", "") + "wise.png"
    if verbose == True: print(figname)
    plt.figure(figsize=(12,8))

    if (verbose == True) & ((np.where(np.isfinite(W13+W13err))[0]).size > 0): print('Max y value: ' , (max((W13+W13err)[np.isfinite(W13+W13err)])+0.1) )
    plt.axis([ 5.0 , 24.0 , \
              max( [(min(np.append((W13-W13err)[ np.isfinite(W13-W13err) ],-0.1))-0.1) , -0.3]) , \
              max( [(max(np.append((W13+W13err)[ np.isfinite(W13+W13err) ],+0.0))+0.2) , +0.6]) ])

    ax1 = plt.gca()
    ax2 = ax1.twiny()
    ax2.set_xlim(5.0,24.0)

    ax1.set_xticks(np.array([5.0 , 10.0 , 15.0 , 17.0 , 22.0 , 24.0]))
    ax1.set_xticklabels(['G5','K0','K5','M0','M5','M7'])
    ax1.set_xlabel('SpT' , fontsize=20, labelpad=15)
    ax1.tick_params(axis='both',which='major',labelsize=16)

    ax2.set_xticks(np.array([5.0 , 10.0 , 15.0 , 17.0 , 22.0 , 24.0]))
    ax2.set_xticklabels(['0.85','0.98','1.45','1.84','3.36','4.75'])
    ax2.set_xlabel(r'$B_p-R_p$ (mag)' , fontsize=20, labelpad=15)
    ax2.tick_params(axis='both',which='major',labelsize=16)

    ax1.set_ylabel(r'$W1-W3$ (mag)' , fontsize=22, labelpad=0)

    # Plot field sequence from Tuc-Hor (Kraus et al. 2014)
    fldspt = [ 5 , 7 , 10 , 12 , 15 , 17 , 20 , 22 , 24 ]
    fldW13 = [ 0 , 0 ,  0 , .02, .06, .12, .27, .40, .60]
    plt.plot(fldspt , fldW13  , zorder=0 , label='Photosphere')

    # Plot neighbors
    ax1.errorbar( spt[yy] , W13[yy] , yerr=W13err[yy] , fmt='none' , ecolor='k')


    for x in range(0 , np.array(yy).size):
        msize  = (17-12.0*(sep3d[yy[x]].value/searchradpc.value))**2
        mcolor = Gchi2[yy[x]]
        medge  = 'black'
        mzorder= 3
        if (r['ruwe'][yy[x]] < 1.2):
            mshape='o'
        if (r['ruwe'][yy[x]] >= 1.2):
            mshape='s'
        if (rvcut != None):
            if (np.isnan(RV[yy[x]])==False) & (np.abs(RV[yy[x]]-Gvrpmllpmbb[yy[x],0]) > rvcut) & (np.abs(RV[yy[x]]-Gvrpmllpmbb[yy[x],0])/RVerr[yy[x]] > 2.0):
                mshape='+'
                mcolor='black'
                mzorder=2
            if (np.isnan(RV[yy[x]])==False) & (np.abs(RV[yy[x]]-Gvrpmllpmbb[yy[x],0]) <= rvcut):
                medge='blue'
        if (mcolor == 'black'):
            ddd = ax1.scatter( [ spt[yy[x]] ] , [ W13[yy[x]] ] , \
                s=msize , c=mcolor , marker=mshape , edgecolors=medge , zorder=mzorder , \
                vmin=0.0 , vmax=vlim.value , cmap='cubehelix' , label='_nolabel' )
        else:
            ccc = ax1.scatter( [ spt[yy[x]] ] , [ W13[yy[x]] ] , \
                s=msize , c=mcolor , marker=mshape , edgecolors=medge , zorder=mzorder , \
                vmin=0.0 , vmax=vlim.value , cmap='cubehelix' , label='_nolabel' )

    temp1 = ax1.scatter([] , [] , c='white' , edgecolors='black', marker='o' , s=12**2 , label = 'RUWE < 1.2')
    temp2 = ax1.scatter([] , [] , c='white' , edgecolors='black', marker='s' , s=12**2 , label = 'RUWE >= 1.2')
    temp3 = ax1.scatter([] , [] , c='white' , edgecolors='blue' , marker='o' , s=12**2 , label = 'RV Comoving')
    temp4 = ax1.scatter([] , [] , c='black' , marker='+' , s=12**2 , label = 'RV Outlier')


    # Plot science target
    if (spt[yy2[0]] > 5):
        plt.plot(spt[yy2[0]] , W13[yy2[0]] , 'rx' , markersize=18 , mew=3 , markeredgecolor='red' , zorder=3 , label=targname )

    plt.legend(fontsize=16 , loc='upper left')
    cb = plt.colorbar(ccc , ax=ax1)
    cb.set_label(label='Velocity Offset (km/s)',fontsize=14)
    if (DoWISE == True): plt.savefig(figname , bbox_inches='tight', pad_inches=0.2 , dpi=200)
    if showplots == True: plt.show()
    plt.close('all')

    # Cross-reference with ROSAT

    v = Vizier(columns=["**", "+_R"] , catalog='J/A+A/588/A103/cat2rxs' )

    zz = np.where( (sep3d.value < searchradpc.value) & (Gchi2 < vlim.value) )
    yy = zz[0][np.argsort(sep3d[zz])]

    ROSATflux = np.empty([np.array(r['ra']).size])
    ROSATflux[:] = np.nan

    print('Searching on neighbors in ROSAT')
    for x in range(0 , np.array(yy).size):
        querycoord = SkyCoord((str(gaiacoord.ra[yy[x]].value) if (gaiacoord.ra[yy[x]].value > 0) else \
                     str(gaiacoord.ra[yy[x]].value+360.0)) , str(gaiacoord.dec[yy[x]].value) , \
                     unit=(u.deg,u.deg) , frame='icrs')
        print('ROSAT query ',x,' of ',np.array(yy).size, end='\r')
        if verbose == True: print('ROSAT query ',x,' of ',np.array(yy).size)
        if verbose == True: print(querycoord)
        if (DoROSAT == True): 
            rosatcat = v.query_region(querycoord , radius='0d1m0s' )
            if (len(rosatcat) > 0):
                rosatcat = rosatcat['J/A+A/588/A103/cat2rxs']
                if verbose == True: print(rosatcat)
                if ((np.where(rosatcat['CRate'] > -999)[0]).size > 0):
                    ww = np.where( (rosatcat['CRate'] == max(rosatcat['CRate'][np.where(rosatcat['CRate'] > -999)])))
                    ROSATflux[yy[x]] = rosatcat['CRate'][ww][0]
                if verbose == True: print(x,yy[x],ROSATflux[yy[x]])


    # Create output table with results
    print('Creating Output Tables with Results')
    if verbose == True: 
        print('Reminder, there were this many input entries: ',len(Gxyz[:,0]))
        print('The search radius in velocity space is: ',vlim)
        print()

    zz = np.where( (sep3d.value < searchradpc.value) & (Gchi2 < vlim.value) )
    sortlist = np.argsort(sep3d[zz])
    yy = zz[0][sortlist]

    fmt1 = "%11.7f %11.7f %6.3f %6.3f %11.3f %8.4f %8.4f %8.2f %8.2f %8.2f %8.3f %4s %8.6f %6.2f %7.3f %7.3f %35s"
    fmt2 = "%11.7f %11.7f %6.3f %6.3f %11.3f %8.4f %8.4f %8.2f %8.2f %8.2f %8.3f %4s %8.6f %6.2f %7.3f %7.3f %35s"
    filename=outdir + targname.replace(" ", "") + ".txt"
    
    warnings.filterwarnings("ignore",category=UserWarning)
    if verbose == True: 
        print('Also creating SIMBAD query table')
        print(filename)
        print('RA            DEC        Gmag   Bp-Rp  Voff(km/s) Sep(deg)   3D(pc) Vr(pred)  Vr(obs)    Vrerr Plx(mas)  SpT    FnuvJ  W1-W3    RUWE  XCrate RVsrc')
    with open(filename,'w') as file1:
        file1.write('RA            DEC        Gmag   Bp-Rp  Voff(km/s) Sep(deg)   3D(pc) Vr(pred)  Vr(obs)    Vrerr Plx(mas)  SpT    FnuvJ  W1-W3    RUWE  XCrate RVsrc \n')
    for x in range(0 , np.array(zz).size):
            if verbose == True:
                print(fmt1 % (gaiacoord.ra[yy[x]].value,gaiacoord.dec[yy[x]].value, \
                  r['phot_g_mean_mag'][yy[x]], r['bp_rp'][yy[x]] , \
                  Gchi2[yy[x]] , sep[yy[x]].value , sep3d[yy[x]].value , \
                  Gvrpmllpmbb[yy[x],0] , RV[yy[x]] , RVerr[yy[x]] , \
                  r['parallax'][yy[x]], \
                  sptstring[yy[x]] , fnuvj[yy[x]] , W13[yy[x]] , r['ruwe'][yy[x]] , ROSATflux[yy[x]] , RVsrc[yy[x]]) )
            with open(filename,'a') as file1:
                  file1.write(fmt2 % (gaiacoord.ra[yy[x]].value,gaiacoord.dec[yy[x]].value, \
                      r['phot_g_mean_mag'][yy[x]], r['bp_rp'][yy[x]] , \
                      Gchi2[yy[x]],sep[yy[x]].value,sep3d[yy[x]].value , \
                      Gvrpmllpmbb[yy[x],0] , RV[yy[x]] , RVerr[yy[x]] , \
                      r['parallax'][yy[x]], \
                      sptstring[yy[x]] , fnuvj[yy[x]] , W13[yy[x]] , r['ruwe'][yy[x]] , ROSATflux[yy[x]] , RVsrc[yy[x]]) )
                  file1.write("\n")

    filename=outdir + targname.replace(" ", "") + ".csv"
    with open(filename,mode='w') as result_file:
        wr = csv.writer(result_file)
        wr.writerow(['RA','DEC','Gmag','Bp-Rp','Voff(km/s)','Sep(deg)','3D(pc)','Vr(pred)','Vr(obs)','Vrerr','Plx(mas)','SpT','FnuvJ','W1-W3','RUWE','XCrate','RVsrc'])
        for x in range(0 , np.array(zz).size):
            wr.writerow(( "{0:.7f}".format(gaiacoord.ra[yy[x]].value) , "{0:.7f}".format(gaiacoord.dec[yy[x]].value) , \
                      "{0:.3f}".format(r['phot_g_mean_mag'][yy[x]]), "{0:.3f}".format(r['bp_rp'][yy[x]]) , \
                      "{0:.3f}".format(Gchi2[yy[x]]) , "{0:.4f}".format(sep[yy[x]].value) , "{0:.4f}".format(sep3d[yy[x]].value) , \
                      "{0:.2f}".format(Gvrpmllpmbb[yy[x],0]) , "{0:.2f}".format(RV[yy[x]]) , "{0:.2f}".format(RVerr[yy[x]]) , \
                      "{0:.3f}".format(r['parallax'][yy[x]]), \
                      sptstring[yy[x]] , "{0:.6f}".format(fnuvj[yy[x]]) , "{0:.2f}".format(W13[yy[x]]) , \
                      "{0:.3f}".format(r['ruwe'][yy[x]]) , "{0:.3f}".format(ROSATflux[yy[x]]) , RVsrc[yy[x]].strip()) )

    if verbose == True: print('All output can be found in ' + outdir)

    return outdir





def binprob(targname,targfilt,targDmag,targDmagerr,targsep,targDPM=None,targDPMerr=None,targDPI=None,Pradvel=None,Pdist=None,Pdisterr=None,PdistU=None,PdistL=None,Pmass=None,PT=None):

### Defining standard color-magnitude and extinction relations

    SpT_Mama = []
    T_Mama   = []
    M_Mama   = []
    R_Mama   = []
    MG_Mama  = []
    MJ_Mama  = []
    MH_Mama  = []
    MK_Mama  = []
    Phot_Mama= []
    Filt_Mama= []

    mamaurl = 'http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt'
    mamafile= urllib.request.urlopen(mamaurl)
    print()
    print('This is querying a table maintained by Eric Mamajek, based on an initial compilation by Pecaut & Mamajek (2013).')
    print('Cite them in your paper, and footnote the URL:')
    print(mamaurl)
    print('Go there and read the warnings. Note, they come after the spells.')
    print()
    print('You should also footnote the date that the compilation was accessed:')

    for i in range(0,3):
        txt = mamafile.readline()
    print(mamafile.readline())
    for i in range(4,22): txt = mamafile.readline()
    hdr = str(mamafile.readline())
    print(hdr)
    for i in range(24,142):
        line = str(mamafile.readline())
        strloc = hdr.find('SpT')
        SpT_Mama.append(           line[(strloc-1):(strloc+4)].strip() )
        strloc = hdr.find('Teff')
        T_Mama.append(      float( line[(strloc+0):(strloc+5)].strip()) )
        strloc = hdr.find('R_Rsun')
        if (line[(strloc+1):(strloc+3)] != '..'):
            R_Mama.append(      float( line[(strloc+0):(strloc+6)].strip()) )
        else:
            R_Mama.append(np.nan)
        strloc = hdr.find('Msun')
        if (line[(strloc+0):(strloc+2)] != '..'):
            M_Mama.append(  float( line[(strloc+0):(strloc+5)].strip()) )
        else:
            M_Mama.append(np.nan)
        strloc = hdr.find('M_G')
        if (line[(strloc+1):(strloc+3)] != '..'):
            MG_Mama.append( float( line[(strloc-1):(strloc+5)].replace(':',' ').strip()) )
        else:
            MG_Mama.append(np.nan)
        strloc = hdr.find('M_Ks')
        if (line[(strloc+0):(strloc+2)] != '..'):
            MK_Mama.append( float( line[(strloc-1):(strloc+5)].strip()) )
        else:
            MK_Mama.append(np.nan)
        strloc = hdr.find('H-Ks')
        if (line[(strloc+0):(strloc+2)] != '..') and (np.isnan(MK_Mama[-1]) == False):
            MH_Mama.append( float( line[(strloc-1):(strloc+5)].strip()) + MK_Mama[-1] )
        else:
            MH_Mama.append(np.nan)
        strloc = hdr.find('J-H')
        if (line[(strloc+0):(strloc+2)] != '..') and (np.isnan(MH_Mama[-1]) == False):
            MJ_Mama.append( float( line[(strloc-1):(strloc+5)].strip()) + MH_Mama[-1] )
        else:
            MJ_Mama.append(np.nan)

    SpT_Mama = np.array(SpT_Mama)
    T_Mama   = np.array(T_Mama)
    R_Mama   = np.array(R_Mama)
    M_Mama   = np.array(M_Mama)
    MG_Mama  = np.array(MG_Mama)
    MJ_Mama  = np.array(MJ_Mama)
    MH_Mama  = np.array(MH_Mama)
    MK_Mama  = np.array(MK_Mama)
    logg_Mama = np.log10( 27400.0 * M_Mama / R_Mama**2)
    print('Done parsing Mamajek table.')

    print('Now parsing color tables from Kraus+2021')
    krausurl = '../../SynthPhot/TableSynColors.txt'
    f = open(krausurl,"r")
    Krausphot = []
    Krausfilt = np.array([ 'G'     , 'Ks'    , 'r'     ,     'i' , 'z'     , 'Bp'    , 'Rp' , 'Kp' , 'LP600' , \
              '[467]' , '[562]' , '[692]' , '[716]' , '[832]' , '[880]' ])
    for s in f:
        Krausphot.append( [ float(s[15:22].strip()) , \
                        (float(s[15:22].strip())-float(s[23:31].strip())) , \
                        float(s[15:22].strip())-float(s[32:40].strip()) , \
                        float(s[15:22].strip())-float(s[41:49].strip()) , \
                        float(s[15:22].strip())-float(s[50:58].strip()) , \
                        float(s[15:22].strip())-float(s[59:67].strip()) , \
                        float(s[15:22].strip())-float(s[68:76].strip()) , \
                        float(s[15:22].strip())-float(s[77:85].strip()) , \
                        float(s[15:22].strip())-float(s[86:94].strip()) , \
                        float(s[15:22].strip())-float(s[95:103].strip()) , \
                        float(s[15:22].strip())-float(s[104:112].strip()) , \
                        float(s[15:22].strip())-float(s[113:121].strip()) , \
                        float(s[15:22].strip())-float(s[122:130].strip()) , \
                        float(s[15:22].strip())-float(s[131:139].strip()) , \
                        float(s[15:22].strip())-float(s[140:148].strip()) ] )
    f.close
    Krausphot = np.array(Krausphot)

    print('Now parsing extinction tables from Kraus+2021')
    krausurl = '../../SynthPhot/TableAXAV.txt'
    f = open(krausurl,"r")
    KrausAXAV = []
    s1=[]
    for s in f:
        s1.append(s[20:25])
    f.close
    s1 = np.array([ 0.0,np.float(s1[1]),np.float(s1[2]),np.float(s1[3]),np.float(s1[5]),np.float(s1[6]),\
                np.float(s1[7]),np.float(s1[8]),np.float(s1[9]),np.float(s1[10]),np.float(s1[11]),\
                np.float(s1[12]),np.float(s1[13]),np.float(s1[14]),0.276,0.176,np.float(s1[4]) ])
    KrausAXAV = np.stack([s1 for n in range(0,72)],axis=0)

    krausurl = '../../SynthPhot/TableATeff.txt'
    f = open(krausurl,"r")
    nline=0
    for s in f:
        KrausAXAV[nline,0] = np.float( s[6:13].strip() )
        KrausAXAV[nline,1] = np.float( s[23:31].strip() )
        KrausAXAV[nline,2] = np.float( s[32:40].strip() )
        KrausAXAV[nline,3] = np.float( s[41:49].strip() )
        KrausAXAV[nline,7] = np.float( s[50:58].strip() )
        nline=nline+1
    f.close

    print('Parsing all to phot table')
    filtarr = np.array([ 'G' , 'Bp' , 'Rp' , 'r' , 'i' , 'z' , 'LP600' , \
              '[467]' , '[562]' , '[692]' , '[716]' , '[832]' , '[880]' , 'J' , 'H' , 'Ks' ])
    photarr = np.zeros( (SpT_Mama.size,filtarr.size) )

    GlocP = np.where( filtarr   == 'G')[0][0]
    KlocP = np.where( filtarr   == 'Ks')[0][0]
    GlocK = np.where( Krausfilt == 'G')[0][0]
    KlocK = np.where( Krausfilt == 'Ks')[0][0]
    for i in np.arange(0,SpT_Mama.size):
        photarr[i,0]  = MG_Mama[i]
        photarr[i,13] = MJ_Mama[i]
        photarr[i,14] = MH_Mama[i]
        photarr[i,15] = MK_Mama[i]
    for i in np.arange(0,SpT_Mama.size):
        for j in np.arange(1,13):
            filtloc = np.where( Krausfilt == filtarr[j])[0][0]
            photarr[i,j]  = np.interp( photarr[i,GlocP] , \
                              Krausphot[:,GlocK] , Krausphot[:,filtloc] , \
                              left = np.nan , right = np.nan )

    print('Parsing all to AXAV table.')
    AXAVarr = np.zeros( (SpT_Mama.size,filtarr.size) )
    for i in np.arange(0,SpT_Mama.size):
        for j in np.arange(0,16):
            AXAVarr[i,j] = np.interp( T_Mama[i] , np.flip(KrausAXAV[:,0]) , np.flip(KrausAXAV[:,(j+1)]))



### Defining Primary Star Properties

    print('Looking up primary star RA/DEC in SIMBAD')
    result_table = customSimbad.query_object(targname)
    ('Target name: ',targname)
    print(result_table['RA','DEC','RVZ_RADVEL','RVZ_ERROR','RVZ_BIBCODE'])
    c = SkyCoord( ra=result_table['RA'][0] , dec=result_table['DEC'][0] , unit=(u.hourangle, u.deg) , frame='icrs')
    print(c)
    print()
    
    print('Look up primary star astrometry/photometry in Gaia FULL DR3')
    sqltext = "SELECT * FROM gaiadr3.gaia_source WHERE CONTAINS( \
        POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec), \
        CIRCLE('ICRS'," + str(c.ra.value) +","+ str(c.dec.value) +","+ str(6.0/3600.0) +"))=1;"
    job = Gaia.launch_job_async(sqltext , dump_to_file=False)
    Pgaia = job.get_results()
    print(Pgaia['ra','dec','phot_g_mean_mag','pmra','pmdec','parallax','ruwe'].pprint_all())
    minpos = Pgaia['phot_g_mean_mag'].tolist().index(min(Pgaia['phot_g_mean_mag']))
    Pcoord = SkyCoord( ra=Pgaia['ra'][minpos]*u.deg , dec=Pgaia['dec'][minpos]*u.deg , \
        distance=1.0*u.parsec , frame='icrs' , \
        radial_velocity=0.0*u.kilometer/u.second , \
        pm_ra_cosdec=Pgaia['pmra'][minpos]*u.mas/u.year , pm_dec=Pgaia['pmdec'][minpos]*u.mas/u.year )
    Pgaia = Pgaia[minpos]
    print(Pcoord)
    print()

    print('Looking up primary star photometry in 2MASS')
    PcoordTM = Pcoord.apply_space_motion(dt=(-15.0*u.year))  
    print(PcoordTM)
    tmass = Irsa.query_region(PcoordTM,catalog='fp_psc' , radius='0d0m10s')
    if ((np.where(tmass['j_m'] > -10.0)[0]).size > 0):
        ww = np.where( (tmass['j_m'] == min(tmass['j_m'][np.where(tmass['j_m'] > 0.0)])))
        Ptmass = tmass[ww][0]
    print(Ptmass['ra','dec','j_m','j_cmsig','h_m','h_cmsig','k_m','k_cmsig','dist','angle'])
    print()

    print('Parsing primary star RV')
    print(result_table['RVZ_RADVEL'].filled(np.nan)[0])
    if (Pradvel != None): 
        print('Using user-provided RV.')
    if (Pradvel == None):
        if (np.isnan(result_table['RVZ_RADVEL'].filled(np.nan)[0]) == False): 
            Pradvel = result_table['RVZ_RADVEL'][0]
            print('Using RV from SIMBAD: ',Pradvel)
        if (np.isnan(result_table['RVZ_RADVEL'].filled(np.nan)[0]) == True) : 
            Pradvel = 0.0
            print('No RV provided, setting to 0.0. Be cautious of projection effects.')
    print('Adopted radvel: ',Pradvel)
    print()

    print('Parsing primary star distance')
    if (Pdist != None): 
        print('Using user-provided distance.')
    if (Pdist == None):
        if ( np.isnan(Pgaia['parallax']) == False):
            v = Vizier(columns=["id","r_med_geo","r_lo_geo","r_hi_geo"] ,\
                            ).query_constraints(id=Pgaia['source_id'],catalog='I/352/gedr3dis')
            vtable = v['I/352/gedr3dis']
            distfrac = (vtable['B_rgeo'][0]-vtable['b_rgeo'][0])/(2.0*vtable['rgeo'][0])
            if (distfrac < 0.1):
                Pdist = vtable['rgeo'][0]
                PdistU= vtable['B_rgeo'][0]- vtable['rgeo'][0]
                PdistL= vtable['rgeo'][0]  - vtable['b_rgeo'][0]
                Pdisterr = 0.5 * (PdistU + PdistL)
                print('Distances of Bailer-Jones21: ',Pdist,Pdisterr,PdistU,PdistL)
    if (Pdist == None):
        print('Can infer primary distance from G-K color. Not implemented yet. Will probably crash imminently.')

    Pcoord = SkyCoord( ra=Pgaia['ra']*u.deg , dec=Pgaia['dec']*u.deg , \
                  distance=Pdist*u.parsec , frame='icrs' , \
                  radial_velocity=Pradvel*u.kilometer/u.second , \
                  pm_ra_cosdec=Pgaia['pmra']*u.mas/u.year , pm_dec=Pgaia['pmdec']*u.mas/u.year )
    print(Pcoord)

    print('Parsing primary star extinction')
    PAV = bayestar(Pcoord , mode='median') * 2.742
    if (np.isnan(PAV) == True): 
        print('NaN was returned, resetting to 0.0')
        PAV = 0.0
    print('Extinction to primary star: ',PAV)

    print('Parsing primary star mass')
    print(Pmass)
    if (Pmass != None):
        print('Using user-provided mass.')
    zz = np.where( (np.isnan(MG_Mama) == False) & (np.isnan(M_Mama) == False) )
    if (Pmass == None):
        PMG = Pgaia['phot_g_mean_mag'] - (5.0*np.log10(Pdist)-5.0) - PAV*0.822
        Pmass = np.interp( PMG , MG_Mama[zz] , M_Mama[zz])
        print('Mass from M_G: ',Pmass)
    zz = np.where( (np.isnan(MK_Mama) == False) & (np.isnan(M_Mama) == False) )
    if (Pmass == None):
        PMK = Ptmass['k_m'] - (5.0*np.log10(Pdist)-5.0) - PAV*0.120
        Pmass = np.interp( PMK , MK_Mama[zz] , M_Mama[zz])
        print('Mass from M_K: ',Pmass)
    if (Pmass == None):
        print('Can infer primary mass from G-K color. Not implemented yet. Will probably crash imminently.')
    print(Pmass)

    print('Parsing primary star temperature')
    if (PT != None):
        print('Using user-privided Teff')
    if (PT == None):
        zz = np.where( (np.isnan(T_Mama) == False) & (np.isnan(M_Mama) == False) )
        PT = np.interp( Pmass , np.flip(M_Mama[zz]) , np.flip(T_Mama[zz]) )
    print(PT)

    print('Parsing primary AX/AV extinction ratios in dmag filters')
    PAXAV = np.zeros(targDmag.size)
    zz = np.where( (np.isnan(T_Mama) == False) )
    for i in range(0,targDmag.size):
        filtloc = np.where( filtarr == targfilt[i] )[0][0]
        PAXAV[i] = np.interp( PT , np.flip(T_Mama[zz]) , np.flip(np.ravel(AXAVarr[zz,filtloc])) )
    print(PAXAV)

    print('Computing primary star magnitudes in dmag filters')
    Pmag = np.zeros(targDmag.size)
    zz = np.where( (np.isnan(T_Mama) == False) )
    filtloc = np.where( filtarr == 'G' )[0][0]
    MamaG = np.interp( PT , np.flip(T_Mama[zz]) , np.flip(photarr[:,filtloc]))
    for i in range(0,targDmag.size):
        filtloc = np.where( filtarr == targfilt[i] )[0][0]
        MamaX = np.interp( PT , np.flip(T_Mama[zz]) , np.flip(photarr[:,filtloc]))
        Pmag[i] = Pgaia['phot_g_mean_mag'] - PAV*0.822 - (MamaG-MamaX) + PAV*PAXAV[i]
    print(Pmag)
    print(targfilt)



#### Simulate binary population

    print('Simulating Binary Population') # Note, someday can make this mass dependent

    Pmass2=1.0
    nsim = 1000000
    print('Primary mass used: ' , Pmass2)
    BF0   = np.interp( Pmass2 , [0.1 , 0.2 , 0.6 , 1.0 , 2.0] , [0.25 , 0.35 , 0.45 , 0.45 , 0.70] )
    gamma = np.interp( Pmass2 , [0.1 , 0.2 , 0.6 , 1.0 , 2.0] , [4.0  , 1.02 , 0.18 , 0.00 ,-2.30] )
    mu    = 1.48 * np.log10(Pmass2) + 2.11
    sigma = -0.93*(np.log10(Pmass2)+0.18)**2 + 1.01

    print('Total binary fraction adopted:   ',BF0)
    print('Mass ratio power law adopted:    ',gamma)
    print('Semimajor axis mu/sigma adopted: ',mu,sigma)

    # Populate random orbital elements

    binfrac = BF0			# Note, assuming BD secondaries aren't included in binary fraction.
    ag , bg = (0.08/Pmass)**(gamma+1.0) , 1.0**(gamma+1.0)
    qq      = ( ag + (bg-ag)*np.random.random(                  nsim ))**(1.0/(gamma+1.0))
    aa      = 10**(np.random.normal(        mu , sigma      , nsim ))

#    binfrac = 0.46 * (1.0 - 0.080/Pmass) / (1.0 - 0.1) # Account for BD secondaries that can't be simulated
#    qq     = np.random.power(      0.080/Pmass , 1.0       , nsim )
#    aa     = 10**(np.random.normal(       1.70 , 1.52      , nsim ))
    ee      = np.random.uniform(           0.0  , 0.95      , nsim )
    littleo = np.random.uniform(           0.0  , np.pi     , nsim )
    bigO    = np.random.uniform(           0.0  , 2.0*np.pi , nsim )
    ii      = np.arccos(np.random.uniform( 0.0  , 1.0       , nsim ))
    MM      = np.random.uniform(           0.0  , 2.0*np.pi , nsim )

    zz = np.where( (qq >= 0.8) & (qq <= 1.0))[0]
    print('Number with q of 0.8 to 1.0: ',zz.size)
    zz = np.where( (qq >= 0.6) & (qq <= 0.8))[0]
    print('Number with q of 0.6 to 0.8: ',zz.size)
    zz = np.where( (qq >= 0.4) & (qq <= 0.6))[0]
    print('Number with q of 0.4 to 0.6: ',zz.size)
    zz = np.where( (qq >= 0.2) & (qq <= 0.4))[0]
    print('Number with q of 0.2 to 0.4: ',zz.size)
    zz = np.where( (qq >= 0.1) & (qq <= 0.2))[0]
    print('Number with q of 0.1 to 0.2: ',zz.size)
    zz = np.where( (qq >= 0.08) & (qq <= 0.1))[0]
    print('Number with q of 0.08 to 0.1: ',zz.size)



    i = np.where( (MM<0.000001) )[0]
    for j in i: MM[j] = MM[j] + 0.000001
    i = np.where( (np.abs(MM-np.pi)<0.000001) )[0]
    for j in i: MM[j] = MM[j] + 0.000002
    i = np.where( (np.abs(MM-2.0*np.pi)<0.000001) )[0]
    for j in i: MM[j] = MM[j] - 0.000001

    print('Computing projected separations from orbital elements')
    EEL  = 0.0
    EE   = MM
    niter= 0
    while np.max(np.abs(EE-EEL)) > 0.00001:
        niter=niter+1
        EEL  = EE
        gEE  = EE - ee*np.sin(EE) - MM
        ggEE = 1.0 - ee*np.cos(EE)
        EE   = EE - gEE/ggEE
    print('Iterations required for Keplers Eqn:',niter)
    rr = aa*(1 - ee*np.cos(EE))
    ff = np.arccos( ((aa*(1.0-ee**2))/rr - 1.0)/ee )
    XX = rr * (np.cos(bigO)*np.cos(littleo+ff) - np.sin(bigO)*np.sin(littleo+ff)*np.cos(ii) )
    YY = rr * (np.sin(bigO)*np.cos(littleo+ff) - np.cos(bigO)*np.sin(littleo+ff)*np.cos(ii) )
    rhoAU = (XX**2 + YY**2)**0.5
    rhoAS = rhoAU/Pdist

    print('Diagnostics')
    print(qq[0:10])
    # Convert mass ratios to mass and Teff, compute AG and AK
    sM   = Pmass*qq
    print(sM[0:10])

    zz = np.where( (np.isnan(M_Mama) == False) & (np.isnan(T_Mama) == False) )
    sT   = np.interp( sM , np.flip(M_Mama[zz]) , np.flip(T_Mama[zz]))
    print(sT[0:10])

    sAXAV = np.zeros([nsim,targfilt.size])
    zz = np.where( (np.isnan(T_Mama) == False) )
    for i in range(0,targfilt.size):
        filtloc = np.where( filtarr == targfilt[i] )[0][0]
        sAXAV[:,i] = np.interp( sT , np.flip(T_Mama[zz]) , np.flip(np.ravel(AXAVarr[zz,filtloc])) )

    filtloc = np.where( filtarr == 'G' )[0][0]
    sAGAV = np.interp( sT , np.flip(T_Mama[zz]) , np.flip(np.ravel(AXAVarr[zz,filtloc])) )
    filtloc = np.where( filtarr == 'Ks' )[0][0]
    sAKAV = np.interp( sT , np.flip(T_Mama[zz]) , np.flip(np.ravel(AXAVarr[zz,filtloc])) )

    SMX = np.zeros([nsim,targfilt.size])
    SX = np.zeros([nsim,targfilt.size])
    dX = np.zeros([nsim,targfilt.size])
    # Compute photometry
    zz = np.where( (np.isnan(T_Mama) == False) )
    for i in range(0,targfilt.size):
        filtloc = np.where( filtarr == targfilt[i] )[0][0]
        SMX[:,i] = np.interp( sT , np.flip(T_Mama[zz]) , np.flip(np.ravel(photarr[zz,filtloc])) )
        SX[:,i]  = np.interp( sT , np.flip(T_Mama[zz]) , np.flip(np.ravel(photarr[zz,filtloc])) ) + \
                                (5.0*np.log10(Pdist)-5.0) - PAV*sAXAV[:,i]
        dX[:,i]  = SX[:,i] - Pmag[i]
    zz = np.where( (np.isnan(T_Mama) == False) & (np.isnan(MG_Mama) == False) )
    sMG   = np.interp( sT , np.flip(T_Mama[zz]) , np.flip(MG_Mama[zz]))
    sG    = sMG + (5.0*np.log10(Pdist)-5.0) + PAV*sAGAV
    with np.printoptions(threshold=np.inf):
        print(SMX[0:10,:])
        print(SX[0:10,:])
        print(dX[0:10,:])

    vorb   = 29.78*np.sqrt((Pmass + sM)/rhoAU)
    if (targDPMerr == None): sPMerr = np.interp( sG , [ 15.0 , 17.0 , 20.0 , 21.0] , [0.03 , 0.07 , 0.5 , 1.4])
    if (targDPMerr != None): sPMerr = targDPMerr * np.ones(vorb.size)
    sPMerr = np.sqrt( sPMerr**2 + (vorb*210.0/Pdist)**2 )
    sPIerr = np.interp( sG , [ 15.0 , 17.0 , 20.0 , 21.0] , [0.03 , 0.07 , 0.5 , 1.3])


### Create binary figures

    if (os.path.isdir('binprob/' + targname) == False): os.mkdir('binprob/' + targname)

    if (targfilt.size > 1):
        for i in range(1,targfilt.size):
            fig,ax1 = plt.subplots(figsize=(9,8))
            ax1.axis([ 0.0 , (max(dX[:,0])+0.2) , 0.0 , (max(dX[:,i])+0.2) ])
            ax1.set_xlabel(r'$\Delta ' + targfilt[0] + '$ (mag)' , fontsize=16)
            ax1.set_ylabel(r'$\Delta ' + targfilt[i] + '$ (mag)' , fontsize=16)
            ax1.tick_params(axis='both',which='major',labelsize=12)

            ccc = ax1.scatter( dX[:,0] + np.random.normal( 0.0 , (targDmagerr[0]+0.1) , nsim ) , \
                           dX[:,i] + np.random.normal( 0.0 , (targDmagerr[i]+0.1) , nsim )  , \
                           s=2 , c='black' , marker='+' , label='Simulated bins')

            plt.plot( targDmag[0] , targDmag[i] , 'rx' , \
                 markersize=12 , mew=3 , markeredgecolor='red' , zorder=3 , label=targname )
            plt.savefig( 'binprob/' + targname + '/' + (targname + '_Bin_'+targfilt[0]+'_'+targfilt[i]+'_D.png') , \
                            bbox_inches='tight', pad_inches=0.2 , dpi=200)
            plt.close('all')

    if (targDPM is not None):
        fig,ax1 = plt.subplots(figsize=(9,8))
        pmsize = np.amax(targDPM) + 5.0
        ax1.axis([ Pgaia['pmra']-pmsize , Pgaia['pmra']+pmsize , Pgaia['pmdec']-pmsize , Pgaia['pmdec']+pmsize ])
        ax1.set_xlabel(r'PMRA (mas)' , fontsize=16)
        ax1.set_ylabel(r'PMDE (mas)' , fontsize=16)
        ax1.tick_params(axis='both',which='major',labelsize=12)

        ccc = ax1.scatter( Pgaia['pmra'] + np.random.normal( 0.0 , sPMerr )  , \
                  Pgaia['pmdec'] + np.random.normal( 0.0 , sPMerr ) , \
                  s=2 , c='black' , marker='+' , label='Simulated bins')
        plt.plot( Pgaia['pmra'] + targDPM[0] , Pgaia['pmdec'] + targDPM[1] , 'rx' , \
                 markersize=12 , mew=3 , markeredgecolor='red' , zorder=3 , label=targname )

        plt.savefig('binprob/' + targname + '/' + targname + '_BinPMD.png',bbox_inches='tight', pad_inches=0.2 , dpi=200)
        plt.close('all')

    if (targDPI is not None):
        fig,ax1 = plt.subplots(figsize=(9,8))
        ax1.axis([ Pgaia['parallax']-4.0 , Pgaia['parallax']+4.0 , 0.0 , (max(dX[:,0])+0.2) ])
        ax1.set_xlabel(r'Parallax (mas)' , fontsize=16)
        ax1.set_ylabel(r'$\Delta ' + targfilt[0] + '$ (mag)' , fontsize=16)
        ax1.tick_params(axis='both',which='major',labelsize=12)

        ccc = ax1.scatter((1000.0/Pdist) + np.random.normal( 0.0 , sPIerr )  , \
                  dX[:,0] + np.random.normal( 0.0 , (targDmagerr[0]+0.1) , nsim ) , \
                  s=2 , c='black' , marker='+' , label='Simulated bins')

        plt.plot( Pgaia['parallax'] + targDPI , targDmag[0] , 'rx' , \
                 markersize=12 , mew=3 , markeredgecolor='red' , zorder=3 , label=targname )


        plt.savefig('binprob/' + targname + '/' + targname + '_BinGPID.png',bbox_inches='tight', pad_inches=0.2 , dpi=200)
        plt.close('all')

    if (targfilt.size > 0): # Note, I think one contrast and a separation are always needed.
    
        fig,ax1 = plt.subplots(figsize=(9,8))
        ax1.axis([ 10**-5.0, (100000.0/Pdist) , 0.0 , (max(dX[:,0])+0.2) ])
        ax1.set_xlabel(r'Proj. Sep. (arcsec)' , fontsize=16)
        ax1.set_ylabel(r'$\Delta ' + targfilt[0] + '$ (mag)' , fontsize=16)
        ax1.tick_params(axis='both',which='major',labelsize=12)
        ax1.set_xscale('log')

        ccc = ax1.scatter(rhoAS  , \
                  dX[:,0] + np.random.normal( 0.0 , (targDmagerr[0]+0.1) , nsim ) , \
                  s=2 , c='black' , marker='+' , label='Simulated bins')
        plt.plot( targsep , targDmag[0] , 'rx' , \
                 markersize=12 , mew=3 , markeredgecolor='red' , zorder=3 , label=targname )

        plt.savefig('binprob/' + targname + '/' + targname + '_BinGrhoD.png',bbox_inches='tight', pad_inches=0.2 , dpi=200)
        plt.close('all')


### Query Gaia for field population

    smin = (40000.0/Pdist)
    smax = max([3600.0 , 206265.0/Pdist])
    print('Search inner/outer radii in arcsec: ',smin,smax)

    # Query Gaia FULL DR3 unless already downloaded
    testname = 'binprob/GaiaDL/' + targname + '_DR3.pickle'
    if ((os.path.isfile(os.path.expanduser(testname))) == True) : 
        print('Already downloaded Gaia DR3 query.')
        with open(testname , 'rb') as gaiadata:
            r = pickle.load(gaiadata)
    if ((os.path.isfile(os.path.expanduser(testname))) == False):
        sqltext = "SELECT * FROM gaiadr3.gaia_source WHERE CONTAINS( \
            POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec), \
            CIRCLE('ICRS'," + str(Pcoord.ra.value) +","+ str(Pcoord.dec.value) +","+ str(smax/3600.0) +"))=1;"
        print(sqltext)
        job = Gaia.launch_job_async(sqltext , dump_to_file=False)
        r = job.get_results()
        with open(testname , 'wb') as gaiadata:
            pickle.dump(r['ra','dec','parallax','parallax_error','parallax_over_error','pmra','pmdec',\
                          'phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag', \
                          'phot_bp_mean_flux_over_error','phot_rp_mean_flux_over_error'] , gaiadata)


    r2 = r
    #print(r[0:10].pprint_all())
    print('Number of entries: ',np.array(r['ra']).size)



### Compute field population properties

    r = r2

    zz = np.where(np.isnan(r['phot_g_mean_mag']) == False)[0]
    print('Number of entries without Gmag: ', (np.array(r['ra']).size - zz.size))
    r = r[zz]

    # Predict mags in other filters
    fT = np.zeros([len(r['ra'])])
    fT[:] = np.nan
    fdist = np.zeros([len(r['ra'])])
    fdist[:] = np.nan
    fAV = np.zeros([len(r['ra'])])
    fAV[:] = np.nan
    fMG = np.zeros([len(r['ra'])])
    fMG[:] = np.nan
    fAGAV = np.zeros([len(r['ra'])])
    fAGAV[:] = 0.822
    fAXAV = np.zeros([len(r['ra']),targfilt.size])
    fAXAV[:,:] = np.nan
    fMX = np.zeros([len(r['ra']),targfilt.size])
    fMX[:,:] = np.nan
    fX = np.zeros([len(r['ra']),targfilt.size])
    fX[:,:] = np.nan
    fdX = np.zeros([len(r['ra']),targfilt.size])
    fdX[:,:] = np.nan
    fdXerr = np.zeros([len(r['ra']),targfilt.size])
    fdXerr[:,:] = np.nan
    fcalctype = np.array(["    " for i in range(len(r['ra']))])
    frhoAS = np.empty(np.array(r['ra']).size)
    frhoAS[:] = np.nan 

    yy = np.where( (np.isnan(fdX[:,0]) == True) & (r['parallax_over_error'] > 10.0))[0]
    print('Number with distance-based properties: ',yy.size)
    if (yy.size > 0):
        fdist[yy] = (1000.0/r['parallax'][yy])
        fcoord = SkyCoord( ra=np.array(r['ra'][yy])*u.deg , dec=np.array(r['dec'][yy])*u.deg , \
                               distance=np.array(fdist[yy])*u.parsec , frame='icrs' )
        frhoAS[yy] = Pcoord.separation(fcoord).to(u.arcsecond).value
        # Rough props
        fAV[yy] = bayestar(fcoord , mode='median') * 2.742
        zz = np.where( np.isnan(fAV[yy]) == True )[0]
        if (zz.size > 0):
            print('Nan is reset to zero for: ',zz.size)
            fAV[yy[zz]] = 0.0

        fMG[yy] = r['phot_g_mean_mag'][yy] - (5.0*np.log10(fdist[yy])-5.0) - fAV[yy]*fAGAV[yy]
        zz = np.where( (np.isnan(MG_Mama) == False) & (np.isnan(T_Mama) == False) )[0]
        fT[yy]  = np.interp( fMG[yy] , MG_Mama[zz] , T_Mama[zz])
        # Compute extinctions
        zz = np.where( (np.isnan(T_Mama) == False) )[0]
        for j in range(0,targfilt.size):
            filtloc = np.where( filtarr == targfilt[j] )[0][0]
            fAXAV[yy,j] = np.interp( fT[yy] , np.flip(T_Mama[zz]) , np.flip(np.ravel(AXAVarr[zz,filtloc])) )
        filtloc = np.where( filtarr == 'G' )[0][0]
        fAGAV[yy] = np.interp( fT[yy] , np.flip(T_Mama[zz]) , np.flip(np.ravel(AXAVarr[zz,filtloc])) )
        # Compute final MG and final temperature
        fMG[yy] = r['phot_g_mean_mag'][yy] - (5.0*np.log10(fdist[yy])-5.0) - fAV[yy]*fAGAV[yy]
        zz = np.where( (np.isnan(MG_Mama) == False) & (np.isnan(T_Mama) == False) )[0]
        fT[yy]   = np.interp( fMG[yy] , MG_Mama[zz] , T_Mama[zz])
        zz = np.where( (np.isnan(T_Mama) == False) )[0]
        for j in range(0,targfilt.size):
            filtloc = np.where( filtarr == targfilt[j] )[0][0]
            fAXAV[yy,j] = np.interp( fT[yy] , np.flip(T_Mama[zz]) , np.flip(np.ravel(AXAVarr[zz,filtloc])) )
        for j in range(0,targfilt.size):
            filtloc = np.where( filtarr == targfilt[j] )[0][0]
            zz = np.where( (np.isnan(T_Mama) == False) & (np.isnan(np.ravel(photarr[:,filtloc])) == False))[0]
            fMX[yy,j] = np.interp( fT[yy] , np.flip(T_Mama[zz]) , \
                            np.flip(np.ravel(photarr[zz,filtloc])) , left=np.nan , right=np.nan )
            fX[yy,j]  = fMX[yy,j] + (5.0*np.log10(fdist[yy])-5.0) - fAV[yy]*fAXAV[yy,j]
            fdX[yy,j] = fX[yy,j] - Pmag[j]
            fdXerr[yy,j] = 0.1
            fcalctype[yy] = 'Dist'

    yy = np.where( (np.isnan(fdX[:,0]) == True) & (r['phot_bp_mean_flux_over_error'][:] > 10.0) & \
                                              (r['phot_rp_mean_flux_over_error'][:] > 10.0))[0]
    print('Number with BpRp-based properties: ',yy.size)
    if (yy.size > 0):     
        # Estimate preliminary fT from color
        filtloc1 = np.where( filtarr == 'Bp' )[0][0]
        filtloc2 = np.where( filtarr == 'Rp' )[0][0]
        fBpRp = (r['phot_bp_mean_mag'][yy] - r['phot_rp_mean_mag'][yy])
        zz = np.where( (np.isnan(np.ravel(photarr[:,filtloc1])) == False) & \
                       (np.isnan(np.ravel(photarr[:,filtloc2])) == False) & \
                       (np.isnan(T_Mama) == False) )
        fT[yy] = np.interp( fBpRp , (np.ravel(photarr[zz,filtloc1]) - np.ravel(photarr[zz,filtloc2])) , T_Mama[zz])
        # Estimate preliminary extinctions coefficients and preliminary distance from fT
        zz = np.where( (np.isnan(MG_Mama) == False) & (np.isnan(T_Mama) == False) )
        filtloc = np.where( filtarr == 'G' )[0][0]
        fAGAV[yy] = np.interp( fT[yy] , np.flip(T_Mama[zz]) , np.flip(np.ravel(AXAVarr[zz,filtloc] )) )
        fABAV     = np.interp( fT[yy] , np.flip(T_Mama[zz]) , np.flip(np.ravel(AXAVarr[zz,filtloc1])) )
        fARAV     = np.interp( fT[yy] , np.flip(T_Mama[zz]) , np.flip(np.ravel(AXAVarr[zz,filtloc2])) )
        for j in range(0,targfilt.size):
            filtloc = np.where( filtarr == targfilt[j] )[0][0]
            fAXAV[yy,j] = np.interp( fT[yy] , np.flip(T_Mama[zz]) , np.flip(np.ravel(AXAVarr[zz,filtloc])) )
        fMG[yy] = np.interp( fT[yy] , np.flip(T_Mama[zz]) , np.flip(MG_Mama[zz]) )
        fdist[yy] = 10**(( r['phot_g_mean_mag'][yy]-fMG[yy]+5.0 )/5.0)
        # Estimate AV and compute projected separation
        fcoord = SkyCoord( ra=np.array(r['ra'][yy])*u.deg , dec=np.array(r['dec'][yy])*u.deg , \
                               distance=np.array(fdist[yy])*u.parsec , frame='icrs' )
        fAV[yy] = bayestar(fcoord , mode='median') * 2.742
        zz = np.where( np.isnan(fAV[yy]) == True )[0]
        if (zz.size > 0):
            print('Nan is reset to zero for: ',zz.size)
            fAV[yy[zz]] = 0.0

        frhoAS[yy] = Pcoord.separation(fcoord).to(u.arcsecond).value
        # Estimate final fT
        zz = np.where( (np.isnan(np.ravel(photarr[:,filtloc1])) == False) & \
                       (np.isnan(np.ravel(photarr[:,filtloc2])) == False) & \
                       (np.isnan(T_Mama) == False) )
        fBpRp = (r['phot_bp_mean_mag'][yy] - r['phot_rp_mean_mag'][yy]) - fAV[yy]*(fABAV-fARAV)
        fT[yy] = np.interp( fBpRp , (np.ravel(photarr[zz,filtloc1]) - np.ravel(photarr[zz,filtloc2])) , T_Mama[zz])
        # Estimate final extinction coefficients and distance
        zz = np.where( (np.isnan(T_Mama) == False) )
        for j in range(0,targfilt.size):
            filtloc = np.where( filtarr == targfilt[j] )[0][0]
            fAXAV[yy,j] = np.interp( fT[yy] , np.flip(T_Mama[zz]) , np.flip(np.ravel(AXAVarr[zz,filtloc])) )
        fMG[yy] = np.interp( fT[yy] , np.flip(T_Mama[zz]) , np.flip(MG_Mama[zz]) )
        fdist[yy] = 10**(( r['phot_g_mean_mag'][yy]-fMG[yy]+fAV[yy]*fAGAV[yy]+5.0 )/5.0)
        fcoord = SkyCoord( ra=np.array(r['ra'][yy])*u.deg , dec=np.array(r['dec'][yy])*u.deg , \
                               distance=np.array(fdist[yy])*u.parsec , frame='icrs' )
        fAV[yy] = bayestar(fcoord , mode='median') * 2.742
        zz = np.where( np.isnan(fAV[yy]) == True )[0]
        if (zz.size > 0):
            print('Nan is reset to zero for: ',zz.size)
            fAV[yy[zz]] = 0.0

        # Compute final abs mag, apparent mag, and contrast
        for j in range(0,targfilt.size):
            filtloc = np.where( filtarr == targfilt[j] )[0][0]
            zz = np.where( (np.isnan(T_Mama) == False) & (np.isnan(np.ravel(photarr[:,filtloc])) == False))[0]
            fMX[yy,j] = np.interp( fT[yy] , np.flip(T_Mama[zz]) , \
                            np.flip(np.ravel(photarr[zz,filtloc])) , left=np.nan , right=np.nan )
            fX[yy,j]  = fMX[yy,j] + (5.0*np.log10(fdist[yy])-5.0) - fAV[yy]*fAXAV[yy,j]
            fdX[yy,j] = fX[yy,j] - Pmag[j]
            fdXerr[yy,j] = 0.1
        fcalctype[yy] = 'BpRp'

    yy = np.where( (np.isnan(fdX[:,0]) == True) & (r['phot_rp_mean_flux_over_error'][:] > 10.0))[0]
    print('Number with mGRp-based properties: ',yy.size)
    if (yy.size > 0):     
        # Estimate preliminary fT from color
        filtloc1 = np.where( filtarr == 'G' )[0][0]
        filtloc2 = np.where( filtarr == 'Rp' )[0][0]
        fmGRp = (r['phot_g_mean_mag'][yy] - r['phot_rp_mean_mag'][yy])
        zz = np.where( (np.isnan(np.ravel(photarr[:,filtloc1])) == False) & \
                       (np.isnan(np.ravel(photarr[:,filtloc2])) == False) & \
                       (np.isnan(T_Mama) == False) )
        fT[yy] = np.interp( fmGRp , (np.ravel(photarr[zz,filtloc1]) - np.ravel(photarr[zz,filtloc2])) , T_Mama[zz])
        # Estimate preliminary extinctions coefficients and preliminary distance from fT
        zz = np.where( (np.isnan(MG_Mama) == False) & (np.isnan(T_Mama) == False) )
        filtloc = np.where( filtarr == 'G' )[0][0]
        fAGAV[yy] = np.interp( fT[yy] , np.flip(T_Mama[zz]) , np.flip(np.ravel(AXAVarr[zz,filtloc] )) )
        fARAV     = np.interp( fT[yy] , np.flip(T_Mama[zz]) , np.flip(np.ravel(AXAVarr[zz,filtloc2])) )
        for j in range(0,targfilt.size):
            filtloc = np.where( filtarr == targfilt[j] )[0][0]
            fAXAV[yy,j] = np.interp( fT[yy] , np.flip(T_Mama[zz]) , np.flip(np.ravel(AXAVarr[zz,filtloc])) )
        fMG[yy] = np.interp( fT[yy] , np.flip(T_Mama[zz]) , np.flip(MG_Mama[zz]) )
        fdist[yy] = 10**(( r['phot_g_mean_mag'][yy]-fMG[yy]+5.0 )/5.0)
        # Estimate AV and compute projected separation
        fcoord = SkyCoord( ra=np.array(r['ra'][yy])*u.deg , dec=np.array(r['dec'][yy])*u.deg , \
                               distance=np.array(fdist[yy])*u.parsec , frame='icrs' )
        fAV[yy] = bayestar(fcoord , mode='median') * 2.742
        zz = np.where( np.isnan(fAV[yy]) == True )[0]
        if (zz.size > 0):
            print('Nan is reset to zero for: ',zz.size)
            fAV[yy[zz]] = 0.0

        frhoAS[yy] = Pcoord.separation(fcoord).to(u.arcsecond).value
        # Estimate final fT
        zz = np.where( (np.isnan(np.ravel(photarr[:,filtloc1])) == False) & \
                       (np.isnan(np.ravel(photarr[:,filtloc2])) == False) & \
                       (np.isnan(T_Mama) == False) )
        fmGRp = (r['phot_g_mean_mag'][yy] - r['phot_rp_mean_mag'][yy]) - fAV[yy]*(fAGAV[yy]-fARAV)
        fT[yy] = np.interp( fmGRp , (np.ravel(photarr[zz,filtloc1]) - np.ravel(photarr[zz,filtloc2])) , T_Mama[zz])
        # Estimate final extinction coefficients and distance
        zz = np.where( (np.isnan(T_Mama) == False) )
        for j in range(0,targfilt.size):
            filtloc = np.where( filtarr == targfilt[j] )[0][0]
            fAXAV[yy,j] = np.interp( fT[yy] , np.flip(T_Mama[zz]) , np.flip(np.ravel(AXAVarr[zz,filtloc])) )
        fMG[yy] = np.interp( fT[yy] , np.flip(T_Mama[zz]) , np.flip(MG_Mama[zz]) )
        fdist[yy] = 10**(( r['phot_g_mean_mag'][yy]-fMG[yy]+fAV[yy]*fAGAV[yy]+5.0 )/5.0)
        fcoord = SkyCoord( ra=np.array(r['ra'][yy])*u.deg , dec=np.array(r['dec'][yy])*u.deg , \
                               distance=np.array(fdist[yy])*u.parsec , frame='icrs' )
        fAV[yy] = bayestar(fcoord , mode='median') * 2.742
        zz = np.where( np.isnan(fAV[yy]) == True )[0]
        if (zz.size > 0):
            print('Nan is reset to zero for: ',zz.size)
            fAV[yy[zz]] = 0.0

        # Compute final abs mag, apparent mag, and contrast
        for j in range(0,targfilt.size):
            filtloc = np.where( filtarr == targfilt[j] )[0][0]
            zz = np.where( (np.isnan(T_Mama) == False) & (np.isnan(np.ravel(photarr[:,filtloc])) == False))[0]
            fMX[yy,j] = np.interp( fT[yy] , np.flip(T_Mama[zz]) , \
                            np.flip(np.ravel(photarr[zz,filtloc])) , left=np.nan , right=np.nan )
            fX[yy,j]  = fMX[yy,j] + (5.0*np.log10(fdist[yy])-5.0) - fAV[yy]*fAXAV[yy,j]
            fdX[yy,j] = fX[yy,j] - Pmag[j]
            fdXerr[yy,j] = 0.1
        fcalctype[yy] = 'mGRp' 

    yy = np.where( (np.isnan(fdX[:,0]) == True) )[0]
    print('Number with Teff-based properties: ',yy.size)
    if (yy.size > 0):     
        fT[yy] = 4500.0
        # Compute extinction coefficients
        zz = np.where( (np.isnan(T_Mama) == False) )
        filtloc = np.where( filtarr == 'G' )[0][0]
        fAGAV[yy] = np.interp( fT[yy] , np.flip(T_Mama[zz]) , np.flip(np.ravel(AXAVarr[zz,filtloc])) )
        for j in range(0,targfilt.size):
            filtloc = np.where( filtarr == targfilt[j] )[0][0]
            fAXAV[yy,j] = np.interp( fT[yy] , np.flip(T_Mama[zz]) , np.flip(np.ravel(AXAVarr[zz,filtloc])) )
        # Compute MG and then fdist
        zz = np.where( (np.isnan(MG_Mama) == False) & (np.isnan(T_Mama) == False) )
        fMG[yy] = np.interp( fT[yy] , np.flip(T_Mama[zz]) , np.flip(MG_Mama[zz]) )
        fdist[yy] = 10**(( r['phot_g_mean_mag'][yy]-fMG[yy]+5.0 )/5.0)
        # Compute actual AV, then update fdist and compute projected sep
        fcoord = SkyCoord( ra=np.array(r['ra'][yy])*u.deg , dec=np.array(r['dec'][yy])*u.deg , \
                               distance=np.array(fdist[yy])*u.parsec , frame='icrs' )
        print(fcoord[0])
        fAV[yy] = bayestar(fcoord , mode='median') * 2.742
        zz = np.where( np.isnan(fAV[yy]) == True )[0]
        if (zz.size > 0):
            print('Nan is reset to zero for: ',zz.size)
            fAV[yy[zz]] = 0.0

        frhoAS[yy] = Pcoord.separation(fcoord).to(u.arcsecond).value
        fdist[yy] = 10**(( r['phot_g_mean_mag'][yy]-fMG[yy]+fAV[yy]*fAGAV[yy]+5.0 )/5.0)
        fcoord = SkyCoord( ra=np.array(r['ra'][yy])*u.deg , dec=np.array(r['dec'][yy])*u.deg , \
                               distance=np.array(fdist[yy])*u.parsec , frame='icrs' )
        fAV[yy] = bayestar(fcoord , mode='median') * 2.742
        zz = np.where( np.isnan(fAV[yy]) == True )[0]
        if (zz.size > 0):
            print('Nan is reset to zero for: ',zz.size)
            fAV[yy[zz]] = 0.0

        # Compute final abs mag, apparent mag, and contrast
        for j in range(0,targfilt.size):
            filtloc = np.where( filtarr == targfilt[j] )[0][0]
            zz = np.where( (np.isnan(T_Mama) == False) & (np.isnan(np.ravel(photarr[:,filtloc])) == False))[0]
            fMX[yy,j] = np.interp( fT[yy] , np.flip(T_Mama[zz]) , \
                            np.flip(np.ravel(photarr[zz,filtloc])) , left=np.nan , right=np.nan )
            fX[yy,j]  = fMX[yy,j] + (5.0*np.log10(fdist[yy])-5.0) - fAV[yy]*fAXAV[yy,j]
            fdX[yy,j] = fX[yy,j] - Pmag[j]
            fdXerr[yy,j] = 0.5
        fcalctype[yy] = 'Teff'


### Field interloper plots

    if (targfilt.size > 1):
        for i in range(1,targfilt.size):
            fig,ax1 = plt.subplots(figsize=(9,8))
            ax1.axis([ 0.0 , (max(fdX[:,0])+0.2) , 0.0 , (max(fdX[:,i])+0.2) ])
            ax1.set_xlabel(r'$\Delta ' + targfilt[0] + '$ (mag)' , fontsize=16)
            ax1.set_ylabel(r'$\Delta ' + targfilt[i] + '$ (mag)' , fontsize=16)
            ax1.tick_params(axis='both',which='major',labelsize=12)

            zz = np.where( fcalctype == 'Teff')[0]
            ccc = ax1.scatter( fdX[zz,0] + np.random.normal( 0.0 , (targDmagerr[0]+fdXerr[zz,0]) , fdX[zz,0].size ) , \
                           fdX[zz,i] + np.random.normal( 0.0 , (targDmagerr[i]+fdXerr[zz,i]) , fdX[zz,0].size)  , \
                           s=2 , c='black' , marker='o' , label='Field (Assumed Teff)')
            zz = np.where( fcalctype == 'mGRp')[0]
            ccc = ax1.scatter( fdX[zz,0] + np.random.normal( 0.0 , (targDmagerr[0]+fdXerr[zz,0]) , fdX[zz,0].size ) , \
                           fdX[zz,i] + np.random.normal( 0.0 , (targDmagerr[i]+fdXerr[zz,i]) , fdX[zz,0].size)  , \
                           s=2 , c='orange' , marker='o' , label='Field (Using G-Rp)')
            zz = np.where( fcalctype == 'BpRp')[0]
            ccc = ax1.scatter( fdX[zz,0] + np.random.normal( 0.0 , (targDmagerr[0]+fdXerr[zz,0]) , fdX[zz,0].size ) , \
                           fdX[zz,i] + np.random.normal( 0.0 , (targDmagerr[i]+fdXerr[zz,i]) , fdX[zz,0].size)  , \
                           s=2 , c='green' , marker='o' , label='Field (Using Bp-Rp)')
            zz = np.where( fcalctype == 'Dist')[0]
            ccc = ax1.scatter( fdX[zz,0] + np.random.normal( 0.0 , (targDmagerr[0]+fdXerr[zz,0]) , fdX[zz,0].size ) , \
                           fdX[zz,i] + np.random.normal( 0.0 , (targDmagerr[i]+fdXerr[zz,i]) , fdX[zz,0].size)  , \
                           s=2 , c='blue' , marker='o' , label='Field (Using Dist)')
        
        
            plt.plot( targDmag[0] , targDmag[i] , 'rx' , \
                 markersize=15 , mew=3 , markeredgecolor='red' , zorder=3 , label=(targname+' cand'))

            lgnd = plt.legend()
            for lh in lgnd.legendHandles:
                lh._sizes = [25.0]
            plt.savefig( 'binprob/' + targname + '/' + (targname + '_Fld_'+targfilt[0]+'_'+targfilt[i]+'_D.png') , \
                            bbox_inches='tight', pad_inches=0.2 , dpi=200)
            plt.close('all')

    if (targDPM is not None):
        fig,ax1 = plt.subplots(figsize=(9,8))
        pmsize = np.amax(targDPM) + 5.0
        ax1.axis([ Pgaia['pmra']-pmsize , Pgaia['pmra']+pmsize , Pgaia['pmdec']-pmsize , Pgaia['pmdec']+pmsize ])
        ax1.set_xlabel(r'PMRA (mas)' , fontsize=16)
        ax1.set_ylabel(r'PMDE (mas)' , fontsize=16)
        ax1.tick_params(axis='both',which='major',labelsize=12)

        ccc = ax1.scatter( r['pmra']  , r['pmdec'] , \
                  s=2 , c='black' , marker='+' , label='Simulated bins')
        plt.plot( Pgaia['pmra'] + targDPM[0] , Pgaia['pmdec'] + targDPM[1] , 'rx' , \
                 markersize=12 , mew=3 , markeredgecolor='red' , zorder=3 , label=targname )

        plt.savefig('binprob/' + targname + '/' + targname + '_FldPMD.png',bbox_inches='tight', pad_inches=0.2 , dpi=200)
        plt.close('all')



    if (targDPI is not None):
        fig,ax1 = plt.subplots(figsize=(9,8))
        ax1.axis([ Pgaia['parallax']-4.0 , Pgaia['parallax']+4.0 , 0.0 , (max(fdX[:,0])+0.2) ])
        ax1.set_xlabel(r'Parallax (mas)' , fontsize=16)
        ax1.set_ylabel(r'$\Delta ' + targfilt[0] + '$ (mag)' , fontsize=16)
        ax1.tick_params(axis='both',which='major',labelsize=12)

        ccc = ax1.scatter(r['parallax']  , \
                  r['phot_g_mean_mag'] - Pgaia['phot_g_mean_mag'] , \
                  s=2 , c='black' , marker='+' , label='Simulated bins')
        plt.plot( Pgaia['parallax'] + targDPI , targDmag[0] , 'rx' , \
                 markersize=12 , mew=3 , markeredgecolor='red' , zorder=3 , label=targname )

        plt.savefig('binprob/' + targname + '/' + targname + '_FldGPID.png',bbox_inches='tight', pad_inches=0.2 , dpi=200)
        plt.close('all')

    print(fdX)
    print(smax)
    if (targfilt.size > 0):
        fig,ax1 = plt.subplots(figsize=(9,8))
        ax1.axis([ 1.0, smax , 0.0 , (max(fdX[:,0])+0.2) ])
        ax1.set_xlabel(r'Proj. Sep. (arcsec)' , fontsize=16)
        ax1.set_ylabel(r'$\Delta ' + targfilt[0] + '$ (mag)' , fontsize=16)
        ax1.tick_params(axis='both',which='major',labelsize=12)
        ax1.set_xscale('log')

        ccc = ax1.scatter(frhoAS  , fdX[:,0] + np.random.normal( 0.0 , (targDmagerr[0]+fdXerr[:,0]) , fdX[:,0].size ) , \
                  s=2 , c='black' , marker='+' , label='Simulated bins')
        plt.plot( targsep , targDmag[0] , 'rx' , \
                 markersize=12 , mew=3 , markeredgecolor='red' , zorder=3 , label=targname )

        plt.savefig('binprob/' + targname + '/' + targname + '_FldGrhoD.png',bbox_inches='tight', pad_inches=0.2 , dpi=200)
        plt.close('all')


### Create vectorized KDE

    logBF = np.zeros(nsim)
    logFF = np.zeros(frhoAS.size)

    # Projected separation term
    DeltalogBF = np.log10(np.e) * ( (-0.5) * ((np.log10(targsep)-np.log10(rhoAS))/0.2)**2 ) / (np.sqrt(2.0*np.pi)*0.2)
    DeltalogFF = np.ones(frhoAS.size) * np.log10((frhoAS.size / (np.pi * (smax**2 - smin**2))) * 2.0 * np.pi * targsep**2 / np.log10(np.e))
    logBF = logBF + DeltalogBF
    logFF = logFF + DeltalogFF
    logBFA = logBF
    logBFP = logBF
    logBFS = logBF
    logFFA = logFF
    logFFP = logFF
    logFFS = logFF

    # Contrast term
    for i in range(0,targfilt.size):
        DeltalogBF = np.log10(np.e) * ( (-0.5) * ((targDmag[i] - dX[:,i])/(targDmagerr[i]+0.2))**2) \
                    / (np.sqrt(2.0*np.pi)*(targDmagerr[i]+0.2))
        DeltalogFF = np.log10(np.e) * ( (-0.5) * ((targDmag[i] - fdX[:,i])/(targDmagerr[i]+fdXerr[:,i]))**2) \
                    / (np.sqrt(2.0*np.pi)*(targDmagerr[i]+fdXerr[:,i]))
        logBF  = logBF  + DeltalogBF
        logFF  = logFF  + DeltalogFF
        logBFP = logBFP + DeltalogBF
        logFFP = logFFP + DeltalogFF
        if (i == 0):
            logBFA = logBFA + DeltalogBF
            logFFA = logFFA + DeltalogFF
            logBFS = logBFS + DeltalogBF
            logFFS = logFFS + DeltalogFF
    
    # Proper motion term
    if (targDPM is not None):
        DeltalogBF = np.log10(np.e) * ( (-0.5) * ((targDPM[0]/sPMerr)**2 + (targDPM[1]/sPMerr)**2)) \
                - np.log10(2.0*np.pi*sPMerr**2)
#        DeltalogFF = np.log10(np.e) * ( (-0.5) * (((Pgaia['pmra']  + targDPM[0] - r['pmra'] )/1.0)**2 + \
#                                              ((Pgaia['pmdec'] + targDPM[1] - r['pmdec'])/1.0)**2)) \
        print('zzzzzz')
        print(np.median(sPMerr))
        print(np.array([1.0 , np.median(sPMerr)]))
        print('Note, KDE bandwidth for field PMD plot is: ',np.amax(np.array([1.0 , np.median(sPMerr)])) )

        DeltalogFF = np.log10(np.e) * ( (-0.5) * (((Pgaia['pmra']  + targDPM[0] - r['pmra'] )/np.amax(np.array([1.0 , np.median(sPMerr)])) )**2 + \
                                              ((Pgaia['pmdec'] + targDPM[1] - r['pmdec'])/np.amax(np.array([1.0 , np.median(sPMerr)])) )**2)) \
                - np.log10(2.0*np.pi*1.0**2)
        logBF  = logBF  + DeltalogBF
        logFF  = logFF  + DeltalogFF
        logBFA = logBFA + DeltalogBF
        logFFA = logFFA + DeltalogFF

    # Parallax term
    if (targDPI is not None):
        DeltalogBF = np.log10(np.e) * ((-0.5)*(targDPI/sPIerr)**2) \
                                        - np.log10(2.0*np.pi*sPIerr**2)
        DeltalogFF = np.log10(np.e) * ((-0.5)*((targDPI+Pgaia['parallax']-r['parallax'])/(0.1+r['parallax_error']))**2) \
                                        - np.log10(2.0*np.pi*1.0**2)
        logBF  = logBF  + DeltalogBF
        logFF  = logFF  + DeltalogFF
        logBFA = logBFA + DeltalogBF
        logFFA = logFFA + DeltalogFF

    BF  = 10**logBF
    BFP = 10**logBFP
    BFA = 10**logBFA
    BFS = 10**logBFS
    FF  = 10**logFF
    FFP = 10**logFFP
    FFA = 10**logFFA
    FFS = 10**logFFS

    yy = np.where(np.isnan(BF) == False)[0]
    BFtot = binfrac*np.sum(BF[yy])   / (yy.size)
    yy = np.where(np.isnan(BFP) == False)[0]
    BFtotP = binfrac*np.sum(BFP[yy]) / (yy.size)
    yy = np.where(np.isnan(BFA) == False)[0]
    BFtotA = binfrac*np.sum(BFA[yy]) / (yy.size)
    yy = np.where(np.isnan(BFS) == False)[0]
    BFtotS = binfrac*np.sum(BFS[yy]) / (yy.size)

    yy = np.where(np.isnan(FF) == False)[0]
    FFtot = np.sum(FF[yy])   / (yy.size)
    yy = np.where(np.isnan(FFP) == False)[0]
    FFtotP = np.sum(FFP[yy]) / (yy.size)
    yy = np.where(np.isnan(FFA) == False)[0]
    FFtotA = np.sum(FFA[yy]) / (yy.size)
    yy = np.where(np.isnan(FFS) == False)[0]
    FFtotS = np.sum(FFS[yy]) / (yy.size)
    
    print('PDF for binary companions: ',BFtot) # units of companions per mag^2 per (mas/yr)^2 per dex of sep
    print('PDF for field interlopers: ',FFtot) # units of companions per mag^2 per (mas/yr)^2 per dex of sep

    print('All: ')
    print('Probability of field:  ',FFtot/(FFtot+BFtot))
    print('Probability of binary: ',BFtot/(FFtot+BFtot))

    print('Survey detection: ')
    print('Probability of field:  ',FFtotS/(FFtotS+BFtotS))
    print('Probability of binary: ',BFtotS/(FFtotS+BFtotS))

    print('Photometry: ')
    print('Probability of field:  ',FFtotP/(FFtotP+BFtotP))
    print('Probability of binary: ',BFtotP/(FFtotP+BFtotP))

    print('Astrometry: ')
    print('Probability of field:  ',FFtotA/(FFtotA+BFtotA))
    print('Probability of binary: ',BFtotA/(FFtotA+BFtotA))

    print('Returning an array of Pbin (all, survey, phot, astro) and then Pfield (same order).')

    return np.array([ [BFtot/(FFtot+BFtot) , BFtotS/(FFtotS+BFtotS) , BFtotP/(FFtotP+BFtotP) , BFtotA/(FFtotA+BFtotA)] , \
                      [FFtot/(FFtot+BFtot) , FFtotS/(FFtotS+BFtotS) , FFtotP/(FFtotP+BFtotP) , FFtotA/(FFtotA+BFtotA)] ])
