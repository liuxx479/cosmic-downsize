from scipy import *
import numpy as np
from scipy import interpolate,stats
from scipy.integrate import quad
import scipy.optimize as op
import sys, os
from scipy.spatial import cKDTree
## rm down.out ; for i in {1..22}; do echo $i; python cosmicdown.py $i >> down.out; done

machine = 'perseus'#'local'#
operation = 1
######### local
if machine == 'local':
    folder = '/Users/jia/cosmicdown/'
    files_arr = sort(os.listdir(folder+'npy/'))

######## perseus
elif machine == 'perseus':
    ##N=int(sys.argv[1])
    from emcee.utils import MPIPool    
    folder=('/home/jialiu/tigress/cosmicdown/')
    files_arr = sort(os.listdir(folder+'npy/'))
    prefix ='Subhalo'#'Behroozi'#''#'Peak'#

h = 0.68
Msun = 1.9891e33
kpc = 3.08567758e18*1e3
kB = 1.38064852e-16 # erg/K
mP = 1.6726219e-24 # grams
G = 6.67408e-8 # 6.67408e-11 m^3 kg^-1 s^-2
yr = 31556926.0 # seconds
Cnfw = 5.0

OmegaM = 0.307115#1#Planck15 TT,TE,EE+lowP+lensing+ext
OmegaV = 1.0-OmegaM
OmegaB = 0.05
fbary = 0.05/OmegaM
muH = 1.18#1/0.76
mu_tot = 0.62#1/(0.76*2+0.24/4*3)

Hcgs = lambda z: 68.0*sqrt(OmegaM*(1+z)**3+OmegaV)*3.24e-20
rho_cz = lambda z: 0.375*Hcgs(z)**2/pi/G
rho_cz0 = rho_cz(0)
n_bary = lambda z: rho_cz(z)*OmegaB/mP*0.76#2.2e-7 #rho_crit * fbary * Omega_m

mdot_fcn = lambda M, z: 46.1*(M/1e12)**1.1*(1+1.11*z)*sqrt(OmegaM*(1+z)**3+(1-OmegaM)) ## unit Msun/yr
Tvir_fcn = lambda Mvir, Rvir: 0.58*G*Mvir*mP/Rvir/kB
V_fcn = lambda Mvir, Rvir: sqrt(G*Mvir/2.0/Rvir)

################### interpolate Mstar from Mhalo include subhalo mass function ###########
logM1C, logepsC, alphaC, delC, gamC = [11.493, -1.600, -2.138, 3.572, 0.487]
logM1S, logepsS, alphaS, delS, gamS = [10.775, -0.957, -2.474, 3.586, 0.423]
logM1B, logepsB, alphaB, delB, gamB = [11.514, -1.777, -1.412, 3.508, 0.316]## Behroozi2013 params 

fC_fcn = lambda x: delC*(log10(1.0+exp(x)))**gamC/(1.0+exp(10.0**(-x))) - log10(10**(alphaC*x)+1.0)
mstarC_fcn = lambda logM: (logepsC+logM1C) + fC_fcn(logM - logM1C) - fC_fcn(0.0)

fS_fcn = lambda x: delS*(log10(1.0+exp(x)))**gamS/(1.0+exp(10.0**(-x))) - log10(10**(alphaS*x)+1.0)
mstarS_fcn = lambda logM: (logepsS+logM1S) + fS_fcn(logM - logM1S) - fS_fcn(0.0)

################### interpolate metallicity from Mstar ##################
mstar_logarr, metal_logarr = genfromtxt('metallicity_mass_table.txt').T [[0,3]]
metal_interp = interpolate.interp1d(mstar_logarr, metal_logarr-8.69, fill_value="extrapolate") ### solar metalicity = 8.69

################## interpolate cooling rate from Temperature and metallicity
cool_table = genfromtxt('coolfunc.dat_orig.txt')
metal_cen_logarr = log10(array([1, 0, 0.316, 0.1, 0.0316, 0.01, 0.001]))
Tcool, Zcool, Lambdacool = array([[cool_table[i,0], metal_cen_logarr[j], cool_table[i,j+1]] for i in range(91) for j in range(7)]).T ## i is the cols for temp, j is for metal
cool_interp = interpolate.NearestNDInterpolator(array([Tcool, Zcool]).T, Lambdacool)
cool_interp1D = interpolate.interp1d(cool_table.T[0], cool_table.T[4], fill_value="extrapolate")

########## Tcrit, above which cannot cool, use actual cooling curve
#Tmax = lambda z, M, ilambdanet: 160.0/3 * M/mdot_fcn(z,M) *ilambdanet *rho_cz0*OmegaB*(1+z)**3 /(muH**2*mP)/kB*yr
#Cnfw_dolag = lambda Mvir, z: 9.59/(1.0+z)*(Mvir/h/1e14)**-0.102
#Edot_heat = lambda Mvir, Rvir, z: 3.0/2.0 * V_fcn(Mvir*Msun, Rvir)**2 *fbary*mdot_fcn(10**logM,z)*Msun/yr
Edot_heat = lambda  Mvir, Rvir, z: 3.0/2.0 *  (G*Mvir/Rvir/2.0) *fbary*mdot_fcn(Mvir/Msun,z)*Msun/yr

fx = lambda x: 1.0 - 1.0/x*log(1+x)
B = 3.0/1.2/(log(2.0) - 0.5)/2.0
ygas = lambda x: exp(-B*fx(x)) 
ygas1x2 = 4*pi*quad(lambda x: ygas(x) *x**2, 0, Cnfw)[0]
ygas2x2 = 4*pi*quad(lambda x: ygas(x)**2 *x**2, 0, Cnfw)[0]

ng_0 = lambda Mvir, Rvir: Mvir*fbary / (ygas1x2  * (Rvir/Cnfw)**3) / mP/muH
Edot_cool = lambda Mvir, Rvir, ilambda: ng_0(Mvir, Rvir)**2*ygas2x2*(Rvir/Cnfw)**3*ilambda

SFR_bins = logspace(-3,3,61)
mvir_bins = logspace (7, 15, 81)
mstar_bins = array([9.3, 9.6, 9.9, 10.2, 10.5, 10.8, 11.0, 15])## all gals
mstar_centers = array([9.45, 9.75, 10.05, 10.35, 10.64, 10.92, 11.2])

def createfile (N):
    
    icat_fn = folder+'npy/'+files_arr[N]
    out_fn = folder+'SFhalos_Suto_Weinberg/%s_M_V_Tcut_ENVonly.npy'%(files_arr[N][:-4])
    a=float(files_arr[N][6:13])
    z = 1/a-1.0
    if not os.path.isfile(out_fn):
        print icat_fn
        #############
        icat = load(icat_fn) # id, mvir, rvir, x, y, z
        idd, mvir, rvir, upid, vrms0 = icat[:,[0,1,2,6,7]].T#, ix, iy, iz
        pos = icat[:,3:6] ### /a to convert to proper distance
        icat = 0 #### release memory
        ##### comoving to proper distance
        pos *= a/h
        rvir *= a/h
        mvir /= h
        ############
        kdt = cKDTree(pos)

        ### cut by temperature (used max point of the cooling curve, ignore now)    
        Tvir_arr = Tvir_fcn(mvir*Msun, rvir*kpc)

        ########## cut by t_cool > t_heat using cooling table
        mstar = mstarC_fcn(log10(mvir))
        mstar[upid!=-1] = mstarS_fcn(log10(mvir[upid!=-1]))## subhalos
        
        ######## fixed lambda
        ##lambda_net = 10**cool_interp1D(log10(Tvir_arr))
        ##lambda_net = 10**-20.84 #Tpeak @ 10**5.3

        metal = metal_interp(mstar)
        ########## D Weinberg response added below 1 line
        metal -= 1.0
        lambda_net = 10.0**cool_interp(log10(Tvir_arr), metal)
        ####################
        #Tcut = Tmax(z, mvir, lambda_net)
        vrms = V_fcn(mvir*Msun, rvir*kpc)/1e5 #km/s,cut at 30km/s
        Edot_heat_arr = Edot_heat (mvir*Msun, rvir*kpc, z)
        Edot_cool_arr = Edot_cool (mvir*Msun, rvir*kpc, lambda_net)
        
        ### Cut out self-heating halos, then remove all halos that're near this halo, idx of halos have T > Tcut, then find neighbors around
        #idx_highT = where(((Tvir_arr-Tcut)>0) & (Tvir_arr>1e5))[0] 
        idx_highT = where((Edot_heat_arr>Edot_cool_arr) & (Tvir_arr>1e5))[0] 
        
        ########## 9/26/2016 ###########
        Tmax_flag = zeros(len(idd)).astype(int)
        Tmax_flag[idx_highT] = 10
        ##########################
        for impact_radius in (1.0, 2.0, 3.0):
            if len(idx_highT)>0:
                remove_lowVcut_array=[]
                for iidx_highT in idx_highT:     
                    ########## Tcut is a number
                    #ir = Tvir_arr[iidx_highT]*rvir[iidx_highT]/Tcut/1e3
                    ########## Tcut = Tcut (M)
                    #ir = (Tvir_arr[iidx_highT]/Tcut[iidx_highT])**0.5*rvir[iidx_highT]/1e3 ## /1e3 because rvir in kpc, pos has unit Mpc
                    #ir = amin([ir, impact_radius*rvir[iidx_highT]/1e3 ])
                    ir = impact_radius*rvir[iidx_highT]/1e3
                    iidx_hiNeibors=array(kdt.query_ball_point(pos[iidx_highT], ir))#
                    ### remove self from ENV impacted objects
                    iidx_hiNeibors = delete(iidx_hiNeibors, where(iidx_hiNeibors==iidx_highT))
                    ##########
                    remove_lowVcut_array.append(iidx_hiNeibors)

                if len(remove_lowVcut_array)>1:
                    remove_lowVcut_array = concatenate(array(remove_lowVcut_array))

                idx_ENV = unique(remove_lowVcut_array)#delete take cares of repeated indices.
                Tmax_flag[idx_ENV] += 1
        save(out_fn, [mvir, vrms, Tmax_flag, mstar])
    
    ######## add new part
    mvir, vrms, flag, mstar = load(out_fn)
    Mc = amax(mvir[flag<10])
    
    SFR = mdot_fcn(mvir, z)
    ################ add smoothed edge ############
    #f_smooth = 0.5*(tanh( (log10(Mc)-0.5-log10(mvir)) * 5.0 )+1)
    #SFR *= f_smooth
    ###############################################
    idx_arr = [range(len(mvir)), #0
            where(mvir<1e12)[0], ##1 exclude self heating with mass cut, don't use  
            where(flag < 10)[0], ##2 exclude self heating
            where(flag%10 <1)[0], ##3 exclude environment only!
            where(flag%10 <2)[0], 
            where(flag%10 <3)[0],
            where(vrms>20)[0], ##6 exclude by velocity
            where(vrms>30)[0], ##7
            where(vrms>40)[0], ##8
            where(vrms>50)[0], ##9
            where(vrms>60)[0], ##10
            where(vrms>70)[0], ##11
            where(vrms>80)[0], ##12
            where(vrms>90)[0], ##13
            where( (vrms>20) & (flag==0))[0], ##14
            where( (vrms>30) & (flag==0))[0], ##15
            where( (vrms>40) & (flag==0))[0],
            where( (vrms>50) & (flag==0))[0],
            where( (vrms>60) & (flag==0))[0],
            where( (vrms>70) & (flag==0))[0],
            where( (vrms>80) & (flag==0))[0],
            where( (vrms>90) & (flag==0))[0],
            where(flag == 10)[0], ##22
            where(flag == 11)[0], ##23
            where(flag == 12)[0], ##24
            where(flag == 13)[0]  ##25
            ]
    SFR_arr = array([histogram(SFR[idx], bins=SFR_bins)[0] for idx in idx_arr])
    Mhalo_arr = array([histogram(mvir[idx], bins=mvir_bins)[0] for idx in idx_arr])
    
    halo_means = array([ [sum(SFR[idx]), mean(SFR[idx]), median(SFR[idx]), sum(mvir[idx]), mean(mvir[idx]), median(mvir[idx]), len(idx)]  for idx in idx_arr])        
    
    SFR_Mstarbin = array([[mean(SFR[idx][ (mstar[idx]>mstar_bins[i])&(mstar[idx]<mstar_bins[i+1]) ]) for i in range(len(mstar_bins)-1)] for idx in idx_arr])
    
    halo_means_massweight = array([ [sum(SFR[idx]), sum(SFR[idx]*mvir[idx])/sum(mvir[idx]), median(SFR[idx]), sum(mvir[idx]), sum(mvir[idx]**2)/sum(mvir[idx]), median(mvir[idx]), len(idx)]  for idx in idx_arr]) 
    
    halo_means_SFRweight = array([ [sum(SFR[idx]), sum(SFR[idx]**2)/sum(SFR[idx]), median(SFR[idx]), sum(mvir[idx]), sum(mvir[idx]*SFR[idx])/sum(SFR[idx]), median(mvir[idx]), len(idx)]  for idx in idx_arr])
    
    return SFR_arr, Mhalo_arr, halo_means, SFR_Mstarbin, halo_means_massweight, halo_means_SFRweight

if operation:
    NN = len(files_arr)
    if machine == 'local':
        map(createfile, range(NN))
    else:
        pool=MPIPool()
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
        statsfolder = 'stats_Suto_Weinberg'#_smoothedge
        all_points = pool.map(createfile, range(NN))
        save(statsfolder + '/SFhalos_SFR_lambda%s.npy'%(prefix), [all_points[i][0] for i in range(NN)])
        save(statsfolder + '/SFhalos_mvir_lambda%s.npy'%(prefix), [all_points[i][1] for i in range(NN)])
        save(statsfolder + '/SFhalos_stats_lambda%s.npy'%(prefix), [all_points[i][2] for i in range(NN)]) 
        save(statsfolder + '/SFhalos_SFRmstarbins_lambda%s.npy'%(prefix), [all_points[i][3] for i in range(NN)])
        save(statsfolder + '/SFhalos_stats_Mweighted_lambda%s.npy'%(prefix), [all_points[i][4] for i in range(NN)])
        save(statsfolder + '/SFhalos_stats_SFRweighted_lambda%s.npy'%(prefix), [all_points[i][5] for i in range(NN)])
        pool.close()

