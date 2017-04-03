from scipy import *
import numpy as np
import os
import sys
from pylab import *
from scipy import optimize,interpolate
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
from scipy.integrate import quad

####### knobs #########
##### official plots ######
get_constants = 1
get_data = 1
plot_f = 0
plot_theory_SFRD = 0
plot_fit_data = 0
plot_downsize = 0
plot_Mquench = 0
plot_int_vs_ext = 1


#### tests ######
Mquench_ratio = 0
plot_SFR = 0
find_f_junk = 0
plot_Luminosity_distribution = 0
plot_massSFR_tacc = 0
#######################
OmegaM = 0.307115 ##### need to correct for the sims
fbary = 0.05/OmegaM
plot_dir = '/Users/jia/Documents/smbhs/CosmicDownsizing/plot/'
fs=22
##### mdot unit: Msun/yr

a_arr = array([ 0.05562,  0.06912,  0.08037,  0.09612,  0.113  ,  0.12987,
        0.15687,  0.18218,  0.23787,  0.26318,  0.2885 ,  0.31381,
        0.33912,  0.36443,  0.38975,  0.41506,  0.44037,  0.46569,
        0.491  ,  0.5315 ,  0.55681,  0.58212,  0.60743,  0.63275,
        0.65806,  0.68337,  0.72134,  0.75931,  0.79728,  0.83525,
        0.87322,  0.91119,  0.94915,  0.98712])
z_arr = 1/a_arr-1
#z_arr = z_arr[1:]
h=0.68

SFR_bins = logspace(-3,3,61)
mvir_bins = logspace (7, 15, 81)
mstar_bins = array([9.3, 9.6, 9.9, 10.2, 10.5, 10.8, 11.0, 15])## all gals
mstar_centers = array([9.45, 9.75, 10.05, 10.35, 10.64, 10.92, 11.2])
#mstar_bins = array([9.4, 9.8, 10.2, 10.6, 11.0, 15])## SF gals

if get_constants:
    h = 0.68
    Msun = 1.9891e33
    kpc = 3.08567758e18*1e3
    kB = 1.38064852e-16 # erg/K
    mP = 1.6726219e-24 # grams
    G = 6.67408e-8 # 6.67408e-11 m^3 kg^-1 s^-2
    yr = 31556926.0 # seconds

    OmegaM = 0.307115#1#Planck15 TT,TE,EE+lowP+lensing+ext
    OmegaV = 1.0-OmegaM
    OmegaB = 0.05
    #OmegaB_z = lambda z: OmegaB*(1+z)**3 / (OmegaM*(1+z)**3 + OmegaV)
    fbary = 0.05/OmegaM
    muH = 1.18#1/0.76
    mu_tot = 0.62#1/(0.76*2+0.24/4*3)
    Hcgs = lambda z: 68.0*sqrt(OmegaM*(1+z)**3+OmegaV)*3.24e-20
    rho_cz = lambda z: 0.375*Hcgs(z)**2/pi/G
    rho_cz0 = rho_cz(0)
    
    c_w, c_m, c_alpha, c_beta, c_gamma = 0.029, 0.097, -110.001, 2469.720, 16.885
    Cnfw = lambda Mvir, z: (Mvir/h) ** (c_w*z-c_m) * 10**(c_alpha/(z+c_gamma)+c_beta/(z+c_gamma)**2) ##### Mvir in unit of Msun    
    Cnfw_KS = lambda logM: 6*(10**logM/1e14/h)**-0.2
    Cnfw_dolag = lambda Mvir, z: 9.59/(1.0+z)*(Mvir/h/1e14)**-0.102
    dd = lambda z: OmegaM*(1+z)**3/(OmegaM*(1+z)**3+OmegaV)
    Delta_vir = lambda z: 18.0*pi**2+82.0*dd(z)-39.0*dd(z)**2#200.0#
    Rvir = lambda Mvir, z: (Mvir*Msun/(4.0/3.0*pi*Delta_vir(z)*rho_cz0*OmegaM*(1.0+z)**3))**0.3333
    ######dont use, wrong##Rvir = lambda M, z: (0.75/200.0*M*Msun/pi/rho_cz(z))**0.3333
    mdot_fcn = lambda M, z: 46.1*(M/1e12)**1.1*(1+1.11*z)*sqrt(OmegaM*(1+z)**3+(1-OmegaM))
    Tvir_fcn = lambda M, z: 0.58*G*M*Msun*mP/Rvir(M, z)/kB
    
    logM1C, logepsC, alphaC, delC, gamC = [11.493, -1.600, -2.138, 3.572, 0.487]
    fC_fcn = lambda x: delC*(log10(1.0+exp(x)))**gamC/(1.0+exp(10.0**(-x))) - log10(10**(alphaC*x)+1.0)
    mstarC_fcn = lambda logM: (logepsC+logM1C) + fC_fcn(logM - logM1C) - fC_fcn(0.0)
    mstar_logarr, metal_logarr = genfromtxt('metallicity_mass_table.txt').T [[0,3]]
    metal_interp = interpolate.interp1d(mstar_logarr, metal_logarr-8.69, fill_value="extrapolate") ### solar metalicity = 8.69
    cool_table = genfromtxt('coolfunc.dat_orig_unsorted.txt')
    cool_table_sorted = genfromtxt('coolfunc.dat_orig_sorted.txt')
    metal_cen_logarr_sorted = array([0.5, 0.0, -0.5, -1.0, -1.5, -2.0, -3.0, -99.0])[::-1]
    metal_cen_logarr = array([0, -99, -0.5, -1, -1.5, -2, -3, 0.5])
    
    Tcool, Zcool, Lambdacool = array([[cool_table[i,0], metal_cen_logarr[j], cool_table[i,j+1]] for i in range(91) for j in range(7)]).T ## i is the cols for temp, j is for metal
    ##cool_interp_rbf = interpolate.Rbf(Tcool, Zcool, Lambdacool,function='linear')
    ##cool_interp_linear = interpolate.LinearNDInterpolator(array([Tcool, Zcool]).T, Lambdacool,rescale=1)#
    #cool_interp_nearest = interpolate.NearestNDInterpolator(array([Tcool, Zcool]).T, Lambdacool)
    cool_interp2d = interpolate.interp2d(Tcool, Zcool, Lambdacool) 
    
    
    cool_interpRect = interpolate.RectBivariateSpline (cool_table_sorted[:,0], metal_cen_logarr_sorted[:-1], cool_table_sorted[:,1:-1], s=2)
    #cool_interp = cool_interpRect##cool_interp_rbf#cool_interp_nearest#cool_interp_linear#

    cool_interp1D = interpolate.interp1d(cool_table.T[0], cool_table.T[4], fill_value="extrapolate")
    
    iMetal = -0.5
    lambda_fcn0 = lambda logM, T: 10.0**cool_interp2d(log10(T), metal_interp(mstarC_fcn(logM)))
    lambda_fcn1 = lambda logM, T: 10.0**cool_interpRect(log10(T), metal_interp(mstarC_fcn(logM))-1.0)
    
    lambda_fcn_group = lambda logM, T: 10.0**cool_interp(log10(T),iMetal)
    #def lambda_fcn(logM, T):
        #ilambda0 = lambda_fcn0(logM, T)
        #ilambda1 = lambda_fcn_group(logM, T)
        ##print ilambda0, ilambda1
        #return amin([ilambda1,ilambda0])
    lambda_fcn = lambda_fcn0
    #lambda_fcn = lambda logM, T: 10**(-23.0-(log10(T)-6.0))
    
    V_fcn = lambda Mvir, Rvir: sqrt(G*Mvir/2.0/Rvir)
if get_data:
    SFRD = load('stats_Suto_Weinberg/SFhalos_stats_lambdaSubhalo.npy')[:,:,0]*(h/250)**3*fbary
    SFRD_smooth = load('stats_Suto_Weinberg_smoothedge/SFhalos_stats_lambdaSubhalo.npy')[:,:,0]*(h/250)**3*fbary#shape (34, 26, 7), z, cuts, 7 stats=[SFRD, SFR_mean, SFR_med, mvir_sum, mvir_mean, mvir_med, Nhalo]
    SFRD_tot = SFRD[:,0]
    SFRD_vrms20 = SFRD[:,6]
    SFRD_vrms30 = SFRD[:,7]
    SFRD_vrms40 = SFRD[:,8]
    SFRD_vrms50 = SFRD[:,9]
    SFRD_vrms80 = SFRD[:,12]
    SFRD_self = SFRD[:,2]
    SFRD_self10 = SFRD[:,22]
    SFRD_self11 = SFRD[:,23]
    SFRD_self12 = SFRD[:,24]
    SFRD_self13 = SFRD[:,25]
    SFRD_env1Rvir = SFRD[:,5]
    SFRD_env2Rvir = SFRD[:,4]
    SFRD_env3Rvir = SFRD[:,3]
    SRFD_env3_self_20km = SFRD[:,14]
    SRFD_env3_self_20km_smooth = SFRD_smooth[:,14]
    SRFD_env3_self_30km = SFRD[:,15]
    SRFD_env3_self_40km = SFRD[:,16]
    SRFD_env3_self_50km = SFRD[:,17]
    SRFD_env3_self_80km = SFRD[:,20]
    z_arr0 = 1+z_arr


if plot_int_vs_ext:
    
    SFRD_mdot, SFRD_ext, SFRD_data = load('SFRD_mdot_ext_data_sharp.npy')
    #loglog(z_arr, SFRD_mdot,label='SFRDmdot')
    #loglog(z_arr, SFRD_ext,label='SFRDext')
    #loglog(z_arr, SFRD_data,label='SFRDdata')
    #legend()
    #savefig('/Users/jia/Desktop/SFRD_text.png');close()
    
    f_ext2data = SFRD_data/SFRD_ext
    #f_mdot2ext = SFRD_ext/SFRD_mdot
    ratio_ext = (SFRD_mdot-SFRD_ext)/(SFRD_mdot-SFRD_data)
    ratio_int = (SFRD_ext-SFRD_data)/(SFRD_mdot-SFRD_data)
    
    for jj in range(2):
        f=figure(figsize=(8,4))
        ax1=f.add_subplot(111)
        if jj==0:
            ax1.fill_between(z_arr0, zeros(len(z_arr0)), ratio_int, color='salmon')
            ax1.fill_between(z_arr0, 1.0-ratio_ext, ones(len(z_arr0)), color='royalblue')
            ax1.set_ylim(0,1)
            ax1.set_title(r'$\boldsymbol{\rm Feedback\; Source}$',fontsize=fs)
            ax1.set_ylabel(r'$\boldsymbol{\rm Fraction}$',fontsize=fs)
            ax1.text(4,0.3,r'$\boldsymbol{\rm Internal}$',fontsize=fs)
            ax1.text(1.3,0.6,r'$\boldsymbol{\rm External}$',fontsize=fs,color='w')
            

        if jj!=0:
            ax1.plot(z_arr0, f_ext2data,color='salmon',lw=4,label=r'$\boldsymbol{f_{\rm int}=\psi_{\rm obs}/\psi_{\rm ext}}$')
            #ax2.plot(z_arr0, f_mdot2ext,'--',color='orange',lw=2)
            xx=linspace(0,6.5,100)
            aa = 0.286
            bb = 6.81
            cc = 0.922
            dd = 0.302
            ee = 0.032
            ff  = 0.750
            ## f_{\rm int}(z) = 0.268(1+0.75z)^6.81\exp(-\frac{z^{0.922}}{0.302})
            ax1.plot(1+xx,aa*(1+ff*xx)**bb*exp(-xx**cc/dd) + ee,'k--',lw=2, label=r'$\boldsymbol{\rm Fitting\; Function}$')
            ax1.set_ylim(0,0.6)
            ax1.set_ylabel(r'$\boldsymbol{f_{\rm int}}$', fontsize=fs)
            ax1.legend(frameon=0,fontsize=fs-2,loc=1)#
        ax1.set_xscale ('log')
        ax1.set_xlim(1.0,7.5)
        ax1.xaxis.set_major_formatter(FormatStrFormatter('%i'))
        ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax1.tick_params(which='both', labelsize=18, width=1.5)
        ax1.set_xticks(arange(1,8))
        a=ax1.get_xticks().tolist()
        anew = [str(ia-1) for ia in a]
        ax1.set_xticklabels(anew)
        [i.set_linewidth(1.5) for i in ax1.spines.itervalues()]
        
        ax1.set_xlabel(r'$\boldsymbol{z}$', fontsize=fs)

        plt.tight_layout()
        #savefig('plot_official/plot_%s.pdf'%(['feedback_ratio','fint'][jj]))
        savefig('plot_official/plot_%s_sharp.png'%(['feedback_ratio','fint'][jj]))
        close()
    
if plot_theory_SFRD:
    

    f = figure(figsize=(16,12))
    ax1 = f.add_subplot(2,2,1)
    ax2 = f.add_subplot(2,2,2)
    ax3 = f.add_subplot(2,2,3)
    ax4 = f.add_subplot(2,2,4)
    
    #f=figure()
    #ax1=f.add_subplot(111)
    for iax in [ax1, ax2, ax3, ax4]:
        iax.plot(z_arr0, SFRD_tot, 'k--', lw=2, label=r'$\boldsymbol{{\rm Prop. \;to \;Mass \;Accretion}}$')
    ax1.plot(z_arr0, SFRD_vrms50, '-', color='dodgerblue', lw=4, label=r'$\boldsymbol{\rm Photoheating\,(\sigma_v\!<\!50 km/s)}$')
    ax1.plot(z_arr0, SFRD_self, '-', color='forestgreen', lw=2, label=r'$\boldsymbol{{\rm Self\mbox{--}heating}  \,(M_h>M_c)}$')
    ax1.plot(z_arr0, SFRD_env3Rvir, '-', color='firebrick', lw=3, label=r'$\boldsymbol{{\rm Hot\;Environment}\,(d\!<\!3r_{\rm vir})}$')
    
    ax2.plot(z_arr0, SFRD_vrms20, '-', color='mediumblue', lw=3, label=r'$\boldsymbol{\rm Photoheating\,(\sigma_v\!<\!20 km/s)}$')
    ax2.plot(z_arr0, SFRD_vrms50, '-', color='dodgerblue', lw=4, label=r'$\boldsymbol{\rm Photoheating\,(\sigma_v\!<\!50 km/s)}$')
    ax2.plot(z_arr0, SFRD_vrms80, '-', color='cyan', lw=2, label=r'$\boldsymbol{\rm Photoheating\,(\sigma_v\!<\!80 km/s)}$')
    
    ax3.plot(z_arr0, SFRD_self, '-', color='forestgreen', lw=2, label=r'$\boldsymbol{{\rm Self\mbox{--}heating}\,(M_h>M_c)}$')
    ax3.plot(z_arr0, SFRD_tot-SFRD_self10, '-', color='limegreen', lw=4, label=r'$\boldsymbol{{\rm Self\mbox{--}heating\,(no\; env.)}}$')
    
    ax4.plot(z_arr0, SFRD_env3Rvir, '-', color='firebrick', lw=3, label=r'$\boldsymbol{{\rm Hot\;Environment}\,(d\!<\!3r_{\rm vir})}$')
    ax4.plot(z_arr0, SFRD_env2Rvir, '-', color='orange', lw=4, label=r'$\boldsymbol{{\rm Hot\;Environment}\,(d\!<\!2r_{\rm vir})}$')
    ax4.plot(z_arr0, SFRD_env1Rvir, '-', color='brown', lw=2, label=r'$\boldsymbol{{\rm Hot\;Environment}\,(d\!<\!1r_{\rm vir})}$')
    ax4.plot(z_arr0, SFRD_env3Rvir+SFRD_self11+SFRD_self12+SFRD_self13, 'x', color='grey', mew=2, label=r'$\boldsymbol{{\rm Hot\;Environment\,(no\;self,}\;3r_{\rm vir})}$')
    boxes = [r'$\boldsymbol{\rm (A) \;Summary}$', 
             r'$\boldsymbol{\rm (B) \;Photoheating}$',
             r'$\boldsymbol{\rm (C) \;Self\mbox{--}heating}$',
             r'$\boldsymbol{\rm (D) \;Hot\;Environment}$']
    colors = ['beige','lightblue','lightgreen','lightpink']#'mistyrose']#['gray','dodgerblue', 'forestgreen', 'firebrick']
    ii=0
    for ax1 in [ax1, ax2, ax3, ax4]:
        ax1.text(0.08, 0.85, boxes[ii],
        verticalalignment='bottom', horizontalalignment='left',
        transform=ax1.transAxes,
        color='k', fontsize=fs, bbox={'facecolor':colors[ii],'alpha':0.8, 'edgecolor':'k','pad':10, 'linewidth':2})
        ii+=1
        
        ax1.set_xscale ('log')
        ax1.set_yscale ('log')
        ax1.set_xlim(1,7.5)
        ax1.set_ylim(2e-2,3)
        ax1.legend(frameon=0,fontsize=fs-2,loc=4)#
        ax1.set_xlabel(r'$\boldsymbol{z}$', fontsize=fs)
        ax1.set_ylabel(r'$\boldsymbol{ \psi\; [M_\odot \,{\rm \, yr}^{-1} \,{\rm \; Mpc}^{-3}]}$',fontsize=fs)
        ax1.xaxis.set_major_formatter(FormatStrFormatter('%i'))
        ax1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax1.set_xticks(arange(1,8))
        a=ax1.get_xticks().tolist()
        anew = [str(ia-1) for ia in a]
        ax1.set_xticklabels(anew)
        [i.set_linewidth(1.5) for i in ax1.spines.itervalues()]
        ax1.tick_params(which='both', labelsize=18, width=1.5)

    plt.tight_layout()
    savefig('plot_official/plot_sims.pdf')
    close()


if plot_Mquench:
    iCnfw = 5.0
    Edot_heat = lambda  Mvir, Rvir, z: 3.0/2.0 *  (G*Mvir/Rvir/2.0) *fbary*mdot_fcn(Mvir/Msun,z)*Msun/yr
    fx = lambda x: 1.0 - 1.0/x*log(1+x)
    B = 3.0/1.2/(log(2.0) - 0.5)/2.0
    ygas = lambda x: exp(-B*fx(x)) 
    
    ygas1x2 = 4*pi*quad(lambda x: ygas(x) *x**2, 0.0, iCnfw)[0]
    ygas2x2 = 4*pi*quad(lambda x: ygas(x)**2 *x**2, 0.0, iCnfw)[0]
            
    ng_0 = lambda Mvir, Rvir: Mvir*fbary / (ygas1x2  * (Rvir/iCnfw)**3) / mP/muH
    #ng_0 = lambda z: Delta_vir(z) * rho_cz0 * OmegaB*(1+z)**3 * (4.0*pi/3.0)*iCnfw**3 / (4.0*pi*quad(lambda x: ygas(x) * x**2, 0, iCnfw)[0])/ mP/muH
    
    Edot_cool = lambda Mvir, Rvir, ilambda, z: ng_0(Mvir, Rvir)**2*ygas2x2*(Rvir/iCnfw)**3*ilambda
    def root_fcn(logM, z, metal0=0, reducemetal=0):
        #print logM, z
        iTvir = Tvir_fcn(10**logM, z)
        if metal0:
            ilambda = 10.0**cool_interp2d(log10(iTvir),-99.0)
        elif reducemetal:
            ilambda = lambda_fcn1(logM, iTvir)
        else:
            ilambda = lambda_fcn0(logM, iTvir)
        iRvir = Rvir(10**logM, z)
        iheat = Edot_heat(10**logM*Msun, iRvir, z)
        icool = Edot_cool(10**logM*Msun, iRvir, ilambda,z)
        return log10(iheat) - log10(icool)
    zz_arr = logspace(0,1,100)-1
    M_arr = zeros(len(zz_arr))
    M_arr_Z0 =  zeros(len(zz_arr))
    M_arr_reducemetal =  zeros(len(zz_arr))
    for jj in range(len(zz_arr)):
        iz=zz_arr[jj]
        iroot = optimize.bisect(root_fcn, 10.0, 15.0, args=(iz,))
        iroot_err = root_fcn(iroot, iz)
        iroot_Z0 = optimize.bisect(root_fcn, 10.0, 15.0, args=(iz,1.0))
        iroot_reducemetal = optimize.bisect(root_fcn, 10.0, 15.0, args=(iz,0.0,1.0))
        M_arr[jj] = iroot
        M_arr_Z0[jj] = iroot_Z0
        M_arr_reducemetal[jj] = iroot_reducemetal
     
    f=figure(figsize=(8,4))
    ax1=f.add_subplot(111)
    ax1.plot(zz_arr+1, M_arr_reducemetal, '-', color='firebrick', lw=4,label=r'$\boldsymbol{\rm Metal \;Cooling}$')
    #ax1.plot(zz_arr+1, M_arr, '-', color='firebrick', lw=4,label=r'$\boldsymbol{\rm Metal \;Cooling}$')
    #ax1.plot(zz_arr+1, M_arr_reducemetal, '--', color='green', lw=2, label=r'$\boldsymbol{\rm Reduced \;Metal \;Cooling}$')
    ax1.plot(zz_arr+1, M_arr_Z0, '--', color='royalblue', lw=2, label=r'$\boldsymbol{\rm No \;Metal \;Cooling}$')
    

    ax1.set_ylabel(r'$\boldsymbol{ \log M_c\; [M_\odot]}$',fontsize=fs)    
    ax1.set_xscale('log')
    ax1.set_xlim(1,7.5)
    ax1.legend(frameon=0,fontsize=fs-2,loc=0)
    ax1.set_xlabel(r'$\boldsymbol{z}$', fontsize=fs)               
    ax1.set_xticks(arange(1,8))
    a=ax1.get_xticks().tolist()
    anew = [str(ia-1) for ia in a]
    ax1.set_xticklabels(anew)
    [i.set_linewidth(1.5) for i in ax1.spines.itervalues()]
    ax1.tick_params(which='both', labelsize=18, width=1.5)
    ax1.set_yticks(arange(12.0, 13.6,0.4))
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))    
    plt.tight_layout()
    savefig('plot_official/plot_Mquench.png')
    close()


if Mquench_ratio:   
    #t_ff = lambda logM, z: 2.0*pi*Rvir(10**logM, z)**1.5/sqrt(G*10**logM*Msun)
    
    Edot_heat = lambda logM, z: 3.0/2.0 * V_fcn(10**logM*Msun, Rvir(10**logM, z))**2 *fbary*mdot_fcn(10**logM,z)*Msun/yr
    ########## Utot
    #N_fcn = lambda logM, z:80*rho_cz0*OmegaB*(1+z)**3 / mP/mu_tot * 4.0/3.0*Rvir(10**logM, z)**3
    #Utot = lambda logM,z: 3/2.0 * N_fcn(logM, z)*kB*Tvir_fcn(10**logM, z)
    #t_heat = lambda logM,z: Utot(logM,z)/Edot_heat(logM,z)
    ##############
    gamma = lambda Cnfw: 1.15+0.01*(Cnfw-6.5)
    mfcn = lambda x: log(1.0+x) - x/(1.0+x)
    eta0 = lambda Cnfw: 0.00676*(Cnfw-6.5)**2 + 0.206*(Cnfw - 6.5) + 2.48
    ###eta0 = lambda logM, z: G*mu_tot*mP*10**logM*Msun/3.0/Rvir(10**logM, z)/kB/Tvir_fcn(10**logM, z)
    ygas_KS = lambda x, iCnfw: (1.0-3.0/eta0(iCnfw)*(1.0-1.0/gamma(iCnfw))*(iCnfw/mfcn(iCnfw)) * (1-log(1+x)/x)) **(1.0/(gamma(iCnfw) - 1.0))
    ########## plot y_gas profile ########
    ######## suto profile
    m1 = log(2.0) - 0.5
    fx = lambda x: 1.0 - 1.0/x*log(1+x)
    igamma = 1.2
    B = 3.0/igamma/m1/2
    rho_x = lambda x, iCnfw: exp(-B*fx(x)) 
    ygas = lambda x, iCnfw: exp(-B*fx(x)) 
    
    #########
    def n2_ratio_fcn (iCnfw):
        integran = lambda x, a, b: ygas(x, iCnfw)**a * x**b
        #n_mean = ygas(iCnfw, iCnfw) / 80.0 * 200.0
        n_mean = quad(integran, 0, iCnfw, args=(1,2))[0] /  quad(integran, 0, iCnfw, args=(0,2))[0]
        ratio = quad(integran, 0, iCnfw, args=(2,2))[0] / quad(integran, 0, iCnfw, args=(0,2))[0] / n_mean**2
        #print iCnfw, ratio
        return ratio
    iCnfw_arr = linspace(2,18,50)
    ratio_arr = array([n2_ratio_fcn(iCnfw) for iCnfw in iCnfw_arr])
    
    def norm_ratio_fcn (iCnfw):
        integran = lambda x: ygas(x, iCnfw) * x**2
        #rho0_tot = 200.0 * rho_cz0 * OmegaB*(1+z)**3 * iCnfw**3 / (3*quad(integran, 0, iCnfw)[0])
        #rho0_boundary = 80.0 *rho_cz0 * OmegaB*(1+z)**3 / ygas(iCnfw, iCnfw)
        #ratio = rho0_tot / rho0_boundary
        ratio = (200.0  * iCnfw**3 / (3*quad(integran, 0, iCnfw)[0])) / (80.0 / ygas(iCnfw, iCnfw))
        return ratio
    iCnfw_arr = linspace(2,18,50)
    norm_ratio_arr = array([norm_ratio_fcn(iCnfw) for iCnfw in iCnfw_arr])
        
    
    ####################################### 
    ##Cnfw_KS = lambda logM: 6*(10**logM/1e14/h)**-0.2
    def Edot_cool(logM, z, int_num): 
        ##iCnfw = Cnfw_KS(10**logM,z)#5.0#
        #iCnfw = Cnfw(10**logM, z)
        iCnfw = Cnfw_dolag(10**logM, z)
        igamma = gamma(iCnfw)
        iRvir = Rvir(10**logM, z)
        iRs = iRvir / iCnfw
        iTvir = Tvir_fcn(10**logM, z)
        ilambda = lambda_fcn(logM, iTvir)
        #iT0 = iTvir*ygas (iCnfw, iCnfw)**(1-igamma)
        #ilambda = lambda x: lambda_fcn(logM, iT0 * ygas (x, iCnfw)**(igamma-1))
        
        rho0 = 80.0*rho_cz0*OmegaB*(1+z)**3/ygas(iCnfw, iCnfw)##rho_cz0*OmegaB*(1+z)**3
        n_gas = lambda x: (rho0/ mP/muH)*ygas (x, iCnfw)  
        
        if int_num == '<n^2>':
            integrand = lambda r: n_gas(r/iRs)**2 * r**2 * ilambda
        elif int_num == 'n_vir^2':            
            integrand = lambda r: n_gas(iCnfw)**2 * r**2 * ilambda   
        elif int_num == '<n>^2':  
            irho_ratio = (200/80.0)**2
            integrand = lambda r: n_gas(iCnfw)**2 * r**2 * ilambda * irho_ratio
        elif int_num == '<n>^2_boost':
            irho_ratio = (200/80.0)**2
            boost_factor = n2_ratio_fcn (iCnfw) 
            integrand = lambda r: n_gas(iCnfw)**2 * r**2 * ilambda * irho_ratio * boost_factor
        elif int_num == '<n>^2_boost_rho':
            irho_ratio = (200/80.0)**2
            irho0_ratio = 1.0/norm_ratio_fcn (iCnfw)**2
            boost_factor = n2_ratio_fcn (iCnfw)
            integrand = lambda r: n_gas(iCnfw)**2 * r**2 * ilambda * irho_ratio * boost_factor * irho0_ratio

        out = 4.0*pi*quad (integrand, 0, iCnfw*iRs, limit=1000)[0]
        return out
    
    root_fcn = lambda logM, z, int_num: log10(Edot_cool(logM, z, int_num))-log10(Edot_heat(logM, z))
    #root_fcn = lambda logM, z: log10(t_ff(logM, z))-log10(t_heat(logM, z))
    
    zz_arr = logspace(0,1,20)-1
    M_arr = zeros(len(zz_arr))
    M_nvir_arr = zeros(len(zz_arr))
    M_nmean_arr = zeros(len(zz_arr))
    M_nmeanboost_arr = zeros(len(zz_arr))
    M_nmeanboost_rho_arr = zeros(len(zz_arr))
    
    for jj in range(len(zz_arr)):
        iz=zz_arr[jj]
        M_arr[jj] = optimize.bisect(root_fcn, 9, 15.0, args=(iz, '<n^2>'))
        M_nvir_arr[jj] = optimize.bisect(root_fcn, 9, 15.0, args=(iz, 'n_vir^2'))
        M_nmean_arr[jj] = optimize.bisect(root_fcn, 9, 15.0, args=(iz, '<n>^2'))
        M_nmeanboost_arr[jj] = optimize.bisect(root_fcn, 9, 15.0, args=(iz, '<n>^2_boost'))
        M_nmeanboost_rho_arr[jj] = optimize.bisect(root_fcn, 9, 15.0, args=(iz, '<n>^2_boost_rho'))
        
        
        #print '%.2f\t%.2f\t(%e)'%(iz, iroot, iroot_err)
    ##################### plot mass profile
    f = figure(figsize=(10,8))
    ax1 = f.add_subplot(221)
    ax2 = f.add_subplot(222)
    ax3 = f.add_subplot(223)
    ax4 = f.add_subplot(224)
    x_arr = logspace(-3,1,100)
    for logM in arange(11,16):
        iCnfw = Cnfw_KS(logM)
        ygas_arr = array([ygas_KS (x, iCnfw) for x in x_arr])
        ax1.loglog(x_arr, ygas_arr,label='logM=%s c=%.2f'%(logM,iCnfw))
    ax1.loglog(x_arr, rho_x(x_arr, iCnfw),'--',label='Suto')
    ax1.legend(frameon=0,loc=0,fontsize=12)
    ax1.set_ylim(1e-3, 10)
    ax1.set_xlabel('x')
    ax1.set_ylabel(r'$\rho_{gas}$')
    ax1.set_title('KS rho')
    
    ax2.plot(iCnfw_arr, ratio_arr,label=r'$ \frac {\int n^2 dV}{\langle n \rangle ^2V}$')
    ax2.set_xlabel(r'$C_{nfw}$')
    ax2.set_ylabel(r'$ \frac {\int n^2 dV}{\langle n \rangle ^2V}$')
    
    ax3.set_xlabel(r'$C_{nfw}$')
    ax3.plot(iCnfw_arr, norm_ratio_arr)
    ax3.set_ylabel(r'$\frac{\rho0_{tot}}{\rho0_{boundary}}$')
    
    ax4.plot(zz_arr+1, M_arr, '-', color='firebrick',lw=2,label=r'$\langle n^2 \rangle, \rho_0\propto\;\rho_{vir}$')
    ax4.plot(zz_arr+1, M_nmeanboost_rho_arr, 'x', color='b',lw=2,label=r'$\langle n \rangle ^2, \rho_0\propto\;\rho_{vir}\times Boost$')
    ax4.plot(zz_arr+1, M_nmeanboost_arr, '--', color='g',lw=2,label=r'$\langle n \rangle ^2, \rho_0\propto\;M_{tot}\times Boost$')
    ax4.plot(zz_arr+1, M_nmean_arr, '--', color='c',lw=2,label=r'$$\langle n \rangle ^2, \rho_0\propto\;M_{tot}$')
    ax4.plot(zz_arr+1, M_nvir_arr, '--', color='purple',lw=2,label=r'$n_{vir} ^2, \rho_0\propto\;M_{tot}$')
    
    
    
    ax4.set_ylabel(r'$\boldsymbol{ \log M_c\; [M_\odot]}$',fontsize=fs)    
    ax4.set_xscale('log')
    ax4.set_xlim(1,8.5)
    ax4.set_ylim(10,14.5)
    ax4.legend(frameon=0,loc=0,fontsize=10)
    ax4.set_xlabel(r'$\boldsymbol{z}$', fontsize=fs)               
    ax4.set_xticks(arange(1,10))
    a=ax4.get_xticks().tolist()
    anew = [str(ia-1) for ia in a]
    ax4.set_xticklabels(anew)
    [i.set_linewidth(1.5) for i in ax1.spines.itervalues()]
    ax4.tick_params(which='both', labelsize=18, width=1.5)
    ax4.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))    
    
    plt.tight_layout()
    savefig('rho_g_suto_Cnfw_dolag.png');close()
     
    #f=figure(figsize=(8,4))
    #ax1=f.add_subplot(111)
    #ax1.plot(zz_arr+1, M_arr, '-', color='firebrick', lw=4)

    #ax1.set_ylabel(r'$\boldsymbol{ \log M_c\; [M_\odot]}$',fontsize=fs)    
    #ax1.set_xscale('log')
    #ax1.set_xlim(1,8.5)
    ##ax1.set_yscale('log')
    ##ax1.set_ylim(0.5e12,6e12)
    ##ax1.set_yticks(arange(11.8,12.6,0.2))
        
    #ax1.set_xlabel(r'$\boldsymbol{z}$', fontsize=fs)               
    #ax1.set_xticks(arange(1,10))
    #a=ax1.get_xticks().tolist()
    #anew = [str(ia-1) for ia in a]
    #ax1.set_xticklabels(anew)
    #[i.set_linewidth(1.5) for i in ax1.spines.itervalues()]
    #ax1.tick_params(which='both', labelsize=18, width=1.5)
    #ax1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))    
    #plt.tight_layout()
    #savefig('plot_Mquench_Nvir_test.png')
    #close()

if plot_f:
    z_bins = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2, 2.5]
    Mstar, zcenter, SFRkarim, z1,z2,M1,M2,errP, errM= genfromtxt('KarimSFRall.txt').T
    mstar_centers = array([9.45, 9.75, 10.05, 10.35, 10.64, 10.92, 11.2])
    SFRinterp = interpolate.interp2d(zcenter, Mstar, SFRkarim, fill_value='nan')
    simSFR=load('stats/SFhalos_SFRmstarbins_lambdaSubhalo.npy')
    fbary= 0.1628
    

    zkarim = array([0.28, 0.48, 0.68, 0.9, 1.1, 1.4, 1.8, 2.3, 2.75])
    newzcenter= array([zkarim[argmin(abs(izcenter-zkarim))] for izcenter in zcenter])

    karim_errM_arr = array([mean((errM/SFRkarim)[newzcenter==iz]) for iz in zkarim])
    karim_errP_arr = array([mean((errP/SFRkarim)[newzcenter==iz]) for iz in zkarim])
    
    errPSFR = interpolate.interp1d(zkarim[:4], karim_errP_arr[:4],bounds_error=0,fill_value="extrapolate")(z_arr)
    errMSFR = interpolate.interp1d(zkarim[:4], karim_errM_arr[:4],fill_value="extrapolate")(z_arr)
    errPSFR[z_arr>0.9]=0.3
    errMSFR[z_arr>0.9]=0.3
    errPSFR = sqrt(errPSFR**2+0.58**2)
    errMSFR = sqrt(errMSFR**2+0.58**2)
    
    zidx=arange(34)
    ffA_arr = zeros((len(zidx), 7))
    ffB_arr = zeros((len(zidx), 7))
    
    for ii in zidx:
        iz = z_arr[ii]
        SFRtheo = SFRinterp (iz, mstar_centers).flatten()
        ff0 = SFRtheo / (simSFR[ii,0]*fbary)
        ffA_arr [ii] = SFRtheo / (simSFR[ii,0]*fbary) ###### all
        ffB_arr [ii] = SFRtheo / (simSFR[ii,14]*fbary) ##### env3, self, 20km/s
    meanffB_arr=mean(ffB_arr,axis=1)
    ff_arr = zeros(len(zidx))
    idx_lim = where((z_arr<1)&(z_arr>0.27))[0]
    
    ff_data = mean(ffB_arr,axis=1)
    ff_arr = 0.3*tanh( 2.5*(z_arr-0.4))+0.45
    
    plot (z_arr[z_arr<0.9], ff_data[z_arr<0.9], label='data')
    plot (z_arr, ff_arr, label='fit')

    #fill_between(z_arr, ff_arr*(1-errMSFR), ff_arr*(1+errPSFR),alpha=0.3)

    legend(frameon=0,loc=2)
    
    xlim(0,3)
    #xlim(0.31,0.9)
    ylim(0.1,1.0)
    xlabel('z')
    ylabel('f')
    #savefig('plot_official/test_f.png')
    close()
    
    #ff_arr=0.4*tanh( 2.0*(z_arr-0.4))+0.45#
    #ff_arr3 = 0.9*(1-exp(-z_arr/.7))+0.1
    #ff_arr = 0.4*tanh( 2.0*(z_arr-0.4))+0.45

if plot_fit_data:
    #ff_arr = 0.3*tanh( 2.5*(z_arr-0.4))+0.45
    ff = 1
    ff_arr = ff*ones(len(z_arr))
    data_arr = genfromtxt('data.txt').T     
    z_data, SFRD_data, errP, errM, IR, author = data_arr
    ### IR: 0=UV, 1=IR, 3=Halpha, 4=UV+IR, 14=1.4GHz
    #z_data, SFRD_data, errP, errM, IR, author = data_arr[:,IR>0]
    z_data0 = z_data+1
    
    fs=22    
    colors = ['dodgerblue','firebrick','limegreen','cyan','orange']
    #symbols = ['o', 'd', 's', 'x', '*', '^', 'v', '+', 'x','p','h','D','H'] 
    symbols = ['o', 'd', 's', 'x', '*', '^', 'v', '<', '>','+','D','.','|']
    
    f=figure(figsize=(8,6))
    ax1=f.add_subplot(111)
    
    
    ax1.plot(z_arr0, SFRD_tot, 'k--', lw=2, label=r'$\boldsymbol{\rm Proportional \;to \;Mass \;Accretion}$')
    #ax1.fill_between(z_arr0, SFRD_tot*(1-errMSFR), SFRD_tot*(1+errPSFR), color='gray',alpha=0.2)
    
    centerrr=SRFD_env3_self_20km
    ax1.plot(z_arr0, centerrr*ff_arr, '-',color='royalblue', lw=4, label=r'$\boldsymbol{ {\rm External\; Feedback\; Only} \;(f_{\rm int}=%s)}$'%(ff))
    ########
    ax1.plot(z_arr0, SRFD_env3_self_20km_smooth, '-',color='cyan', lw=4, label=r'$\boldsymbol{ {\rm External\; Feedback\; Only} \;(f_{\rm int}=%s, smooth)}$'%(ff))
    
    #ax1.fill_between(z_arr0, centerrr*ff_arr*(1-errMSFR), centerrr*ff_arr*(1+errPSFR), color='deepskyblue',alpha=0.3)#, 'b-', lw=2, label=r'$\boldsymbol{\rm After\;cut (30km/s, f=tanh\,fcn)}$')
    ############## plot Madau bestfit #######
    #best_fit_madau = lambda z: 0.015*(1+z)**2.7/(1+((1+z)/2.9)**5.6)
    #best_fit_madau_test = lambda A, B, C, D, z: A*(1+z)**B/(1+((1+z)/C)**D)
    #chisq_madau = lambda ABCD: sum((best_fit_madau_test(ABCD[0],ABCD[1],ABCD[2],ABCD[3],z_data)-SFRD_data)**2/(errM+errP)**2)
    #A,B,C,D=optimize.minimize(chisq_madau, (0.015, 2.7, 2.9, 5.6)).x
    A,B,C,D = 0.01594907,  2.94688778,  2.70011684,  5.9136489
    best_fit_madau = lambda z: A*(1+z)**B/(1+((1+z)/C)**D)
    ax1.plot(z_arr0, best_fit_madau(z_arr0-1.0),'-',color='salmon',lw=2,label=r'$\boldsymbol{\rm Best\; Fit\; to\; Data }$')
    
    ########################################

    iii=-1
    for iIR in unique(IR):
        iii+=1
        authors = unique(author[(author<100) & (IR==iIR)])
        jjj=-1
        for iauthor in authors:
            jjj+=1
            idx = where((IR==iIR)&(author==iauthor))[0]
            ax1.errorbar(z_data0[idx], SFRD_data[idx], yerr=[errP[idx], errM[idx]], ecolor=colors[iii], mec=colors[iii], mfc=colors[iii], mew=2, ms=8, ls='None', fmt=symbols[jjj], capsize=0)
            ax1.set_xscale ('log')
    
    
    ax1.legend(frameon=0,fontsize=fs-3,loc=2)#
    
    ax1.set_yscale ('log')
    ax1.set_xlim(1,7.5)
    ax1.set_ylim(8e-3,3.5)#1.5)
    
    ax1.set_xlabel(r'$\boldsymbol{z}$', fontsize=fs)
    ax1.set_ylabel(r'$\boldsymbol{ \psi\; [M_\odot \,{\rm \, yr}^{-1} \,{\rm \; Mpc}^{-3}]}$',fontsize=fs)
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%i'))
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax1.set_xticks(arange(1,8))
    a=ax1.get_xticks().tolist()
    anew = [str(ia-1) for ia in a]
    ax1.set_xticklabels(anew)
    [i.set_linewidth(1.5) for i in ax1.spines.itervalues()]
    ax1.tick_params(which='both', labelsize=18, width=1.5)
    plt.tight_layout()
    #savefig('plot_official/plot_data20_smooth.pdf')
    savefig('plot_official/plot_data20_smoothtest.png')
    close()
    #save('SFRD_mdot_ext_data.npy',[SFRD_tot,centerrr, best_fit_madau(z_arr0-1.0)])
    
if plot_downsize:
    ########
    ff_arr = ones(len(z_arr0))
    #########
    SFR_bins = logspace(-3,3,61)
    mvir_bins = logspace (7, 15, 81)

    from matplotlib.ticker import OldScalarFormatter, ScalarFormatter
    idxselect = 14 # 14 is velocity 20km/s
    #SFR_hist = load('stats_Suto/SFhalos_SFR_lambdaSubhalo.npy')[:,idxselect]
    #mvir_hist = load('stats_Suto/SFhalos_mvir_lambdaSubhalo.npy')[:,idxselect]
    #SFR_typical = SFR_bins[argmax(SFR_hist,axis=1)]
    #mvir_typical = mvir_bins[argmax(mvir_hist,axis=1)]
    folder='stats_Suto_Weinberg_smoothedge/'
    SFR_mean = load(folder+'SFhalos_stats_lambdaSubhalo.npy')[:,idxselect,1] #mean is 1
    Mhalo_mean = load(folder+'SFhalos_stats_lambdaSubhalo.npy')[:,idxselect,4] #mean is 4
    SFR_all = load(folder+'SFhalos_stats_lambdaSubhalo.npy')[:,0, 1] #mean is 1
    Mhalo_all = load(folder+'SFhalos_stats_lambdaSubhalo.npy')[:,0, 4] #mean is 4
    SFR_mean_Mweighted = load(folder+'SFhalos_stats_Mweighted_lambdaSubhalo.npy')[:,idxselect,1] #mean is 1
    Mhalo_mean_Mweighted = load(folder+'SFhalos_stats_Mweighted_lambdaSubhalo.npy')[:,idxselect,4] #mean is 4
    
    SFR_mean_SFRweighted = load(folder+'SFhalos_stats_SFRweighted_lambdaSubhalo.npy')[:,idxselect,1] #mean is 1
    Mhalo_mean_SFRweighted = load(folder+'SFhalos_stats_SFRweighted_lambdaSubhalo.npy')[:,idxselect,4] #mean is 4
    
    
    f=figure(figsize=(8,5))
    ax1=f.add_subplot(111)
    ax2 = ax1.twinx()

    lns1=ax1.plot(z_arr0, SFR_mean*fbary*ff_arr, '--', color='firebrick', lw=4,label=r'$\boldsymbol{\rm Number\mbox{-}weighted\; \langle \dot{M}_\star \rangle}$')
    #ax1.fill_between(z_arr0, SFR_mean*fbary*ff_arr*(1-errMSFR), SFR_mean*fbary*ff_arr*(1+errPSFR), color='firebrick', alpha=0.3)
    lns2=ax2.plot(z_arr0, Mhalo_mean, '-', color='mediumblue', lw=4, label=r'$\boldsymbol{ {\rm Number\mbox{-}weighted}\; \langle M_h \rangle}$')
        
    ########## mass weighted
    #ax1.plot(z_arr0, SFR_mean_Mweighted*fbary*ff_arr, '--', color='firebrick', lw=2)
    #ax1.fill_between(z_arr0, SFR_mean_Mweighted *fbary*ff_arr*(1-errMSFR), SFR_mean_Mweighted*fbary*ff_arr*(1+errPSFR), color='firebrick', alpha=0.3)
    #ax2.plot(z_arr0, Mhalo_mean_Mweighted, '-', color='mediumblue', lw=2, label='mass weighted')

    lns3=ax1.plot(z_arr0, SFR_mean_SFRweighted*fbary*ff_arr, '--', color='firebrick', lw=1, label=r'$\boldsymbol{\rm SFR\mbox{-}weighted\; \langle \dot{M}_\star \rangle}$')
    #ax1.fill_between(z_arr0, SFR_mean_SFRweighted *fbary*ff_arr*(1-errMSFR), SFR_mean_SFRweighted*fbary*ff_arr*(1+errPSFR), color='firebrick', alpha=0.3)
    lns4=ax2.plot(z_arr0, Mhalo_mean_SFRweighted, '-', color='mediumblue', lw=1, label=r'$\boldsymbol{ {\rm SFR\mbox{-}weighted}\; \langle M_h \rangle}$')

    #
    #ax2.legend(frameon=0,fontsize=fs-2,loc=2)#
    lns = lns1+lns2+lns3+lns4
    labs = [l.get_label() for l in lns]
    #ax.legend(lns, labs, loc=0)
    ax1.legend(lns, labs, ncol=2, fontsize=13,loc=9,frameon=1)#, ,loc=2)#
    
    ax1.set_ylabel(r'$\boldsymbol{ \langle \dot{M}_\star\rangle \; [M_\odot \,{\rm \, yr}^{-1} ]}$',fontsize=fs, color='firebrick')
    ax2.set_ylabel(r'$\boldsymbol{ \langle M_h\rangle\; [M_\odot]}$',fontsize=fs,  rotation=270, labelpad=30, color='mediumblue')
 
    for tl in ax1.get_yticklabels():
        tl.set_color('firebrick')
    for tl in ax2.get_yticklabels():
        tl.set_color('mediumblue')
    
    ax1.set_ylim(0.1, 200)    
    ax2.set_ylim(3e9, 3e13)
    
    for iax in [ax1,ax2]:
        iax.set_xscale('log')
        
        iax.set_xlim(1,7.5)
        iax.set_xlabel(r'$\boldsymbol{z}$', fontsize=fs)
        
        if iax==ax2:
            iax.set_xticks(arange(1,8))
            a=iax.get_xticks().tolist()
            anew = [str(ia-1) for ia in a]
            iax.set_xticklabels(anew)
            [i.set_linewidth(1.5) for i in iax.spines.itervalues()]
        iax.tick_params(which='both', labelsize=18, width=1.5)
        iax.set_yscale ('log')
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    plt.tight_layout()
    #ax1.set_title('smoothe edge')
    savefig('plot_official/plot_downsize20.pdf')
    close()

    
    
    
    






################## test plots ###############


#mdot_fcn = lambda z, M: 46.1*(M/1e12)**1.1*(1+1.11*z)*sqrt(OmegaM*(1+z)**3+(1-OmegaM))

if find_f_junk:
    Mhalo_params_arr = [[12.520, 10.916, 0.457, 0.566, 1.53],
            [12.725, 11.038, 0.466, 0.610, 1.95],
            [12.722, 11.100, 0.470, 0.393, 2.51]]
    ##log10M1, log10M0, beta, sig, gamma
    redshift_edges=[[0, 0.48], [0.48,0.74], [0.74, 1.30]]
    Mstar_arr = linspace(5, 13, 100)
    def Mstar2Mhalo (log10Mstar, i=0):
        '''input log10 (Mstar_arr/Msun), return Mhalo with unit: log10(Mhalo/Mstar)
        i = 0,1,2 is the redshift bin
        '''
        z0,z1 = redshift_edges[i]
        log10M1, log10M0, beta, sig, gamma = Mhalo_params_arr[i]
            #print log10M1, log10M0, beta, sig, gamma 
        Mhalo = log10M1+beta*(log10Mstar-log10M0)+10.0**(sig*(log10Mstar-log10M0))/(1+10.0**(-gamma*(log10Mstar-log10M0)))-0.5
        return Mhalo
    Mhalo_interp0 = interpolate.interp1d(Mstar2Mhalo(Mstar_arr,i=0), Mstar_arr,fill_value='extrapolate')
    Mhalo_interp1 = interpolate.interp1d(Mstar2Mhalo(Mstar_arr,i=1), Mstar_arr,fill_value='extrapolate')
    Mhalo_interp2 = interpolate.interp1d(Mstar2Mhalo(Mstar_arr,i=2), Mstar_arr,fill_value='extrapolate')
    
    def return_points(ifile):        
        print ifile
        a = float(ifile[6:13])    
        z = 1/a-1.0
        mvir,  sigma, flag = load(folder+ifile)
        SFR = mdot_fcn(z, mvir)
        
        test_bin = z - array([0, 0.48, 0.74])
        i = where (test_bin>0)[0][-1]
        Mhalo_interp = [Mhalo_interp0, Mhalo_interp1, Mhalo_interp2][i]
        mstar = Mhalo_interp(log10(mvir))
        
        #idx_arr = [where((flag%10 ==0 )&(sigma>isigma)& (mvir<1e12))[0] for isigma in range(20,90,10)]
        idx_arr = [where((flag == 0 )&(sigma>isigma))[0] for isigma in range(20,90,10)]
        idx_arr.append(arange(len(mstar)))

        SFR_arr = array([histogram(SFR[idx], bins=SFR_bins)[0] for idx in idx_arr])
        Mhalo_arr = array([histogram(mvir[idx], bins=mvir_bins)[0] for idx in idx_arr])
        halo_means = array([ [sum(SFR[idx]), mean(SFR[idx]), median(SFR[idx]), sum(mvir[idx]), mean(mvir[idx]), median(mvir[idx]), len(idx)]  for idx in idx_arr])        
        SFR_Mstarbin = array([[mean(SFR[idx][ (mstar[idx]>mstar_bins[i])&(mstar[idx]<mstar_bins[i+1]) ]) for i in range(len(mstar_bins)-1)] for idx in idx_arr])
        
        return SFR_arr, Mhalo_arr, halo_means, SFR_Mstarbin
    
    all_points = map(return_points, files)
    save('SFhalos_SFR_noMcut_lambda.npy', [all_points[i][0] for i in range(len(files))])
    save('SFhalos_mvir_noMcut_lambda.npy', [all_points[i][1] for i in range(len(files))])
    save('SFhalos_stats_noMcut_lambda.npy', [all_points[i][2] for i in range(len(files))]) # redshift, sigma cuts, stats [22,7,8]
    save('SFhalos_SFRmstarbins_noMcut_lambda.npy', [all_points[i][3] for i in range(len(files))])
    
    
if plot_SFR:
    def return_points(ifile):
        print ifile
        a = float(ifile[6:13])    
        z = 1/a-1.0
        mvir,  sigma, flag = load(folder+ifile)
        ###### flag 0: not excluded by any Tmax
        ###### flag 13: excluded by self Tmax
        ###### flat 1/2/3: exluded by Hot\;Environment, 1, 2, 3 virial radiis
        SFR = mdot_fcn(z, mvir)
        
        #if z>=3:
            #sigma_crit = 50.0
        #else:
            #sigma_crit = 50.0*sqrt((1+z)/3.0)
        ##SFR =  mdot_fcn(z, mvir)*exp(-0.69*(sigma_crit/sigma)**2)
        ##SFR = SFR[where( (flag==0) & (sigma>sigma_crit) )[0]]
        #mvir = mvir[where( (flag==0) & (sigma>sigma_crit) )[0]]
        #SFR_arr = histogram(SFR, bins=SFR_bins)[0]
        #Mhalo_arr = histogram(mvir, bins=mvir_bins)[0]
        #halo_means = array([sum(SFR), mean(SFR), median(SFR), sum(mvir), mean(mvir), median(mvir)] )
        
        #idx_arr = [range(len(mvir)),
                   #where(sigma>30)[0],
                   #where(sigma>50)[0],
                   #where(sigma>80)[0],
                   #where(flag<13)[0], ## 4, exclude self Tmax halo
                   #where(flag<3)[0], ## 5, exclude self Tmax 1 virial
                   #where(flag<2)[0], ## 6, 2 virial
                   #where(flag==0)[0], ## 7, 3 virial
                   #where((flag==0)&(sigma>30))[0]]
        idx_arr = [range(len(mvir)),
            where(sigma>30)[0],
            where(sigma>50)[0],
            where(sigma>80)[0],
            where(mvir<1e12)[0], ## exclude large self   
            where(flag != 10)[0],
            where(flag%10 <1)[0], ## exclude self
            where(flag%10 <2)[0], 
            where(flag%10 <3)[0],
            where((sigma>50 )&(mvir<1e12))[0],
            where((flag%10 ==0 )&(sigma>50))[0],
            where((flag%10 ==0 )&(mvir<1e12))[0],
            where((flag%10 ==0 )&(sigma>50)& (mvir<1e12))[0]]
        ## sigma_arr is from 10,20,..100
        
        SFR_arr = array([histogram(SFR[idx], bins=SFR_bins)[0] for idx in idx_arr])
        Mhalo_arr = array([histogram(mvir[idx], bins=mvir_bins)[0] for idx in idx_arr])
        halo_means = array([ [sum(SFR[idx]), mean(SFR[idx]), median(SFR[idx]), sum(mvir[idx]), mean(mvir[idx]), median(mvir[idx]), len(idx)]  for idx in idx_arr])
        
        return SFR_arr, Mhalo_arr, halo_means#halo_means#
 
    #all_points = map(return_points, files)
    #save('SFhalos_SFR_M11_ENVonly.npy', [all_points[i][0] for i in range(len(files))])
    #save('SFhalos_mvir_M11_ENVonly.npy', [all_points[i][1] for i in range(len(files))])
    #save('SFhalos_stats_M11_ENVonly.npy', [all_points[i][2] for i in range(len(files))])
    

    #SFhalos_SFR  = load('SFhalos_SFR_M11.npy')
    #SFhalos_mvir = load('SFhalos_mvir_M11.npy')             
    SFhalos_stats= load('SFhalos_stats_M11_ENVonly.npy')   

    legs = [['All', '*'], 
            ['v>30', '--'],
            ['v>50', '--'],
            ['v>80', '--'],
            ['mvir<1e12', '.'],
            ['Tmax_self', '.'],
            ['3Rvir', 'x'],
            ['2Rvir', 'x'],
            ['1Rvir', 'x'],
            ['v>50, mvir<1e12', '-'],
            ['v>50, 3Rvir', '-'],
            ['mvir<1e12, 3Rvir', '-'],
            ['v>50, mvir<1e12, 3Rvir', '-']]
    
    #legs = [['All', '*'], 
            #['v>30', '--'],
            #['v>50', '--'],
            #['v>80', '--'],
            #['Tmax', '.'],
            #['1Rvir', 'x'],
            #['2Rvir', 'x'],
            #['3Rvir', 'x'],
            #['v>30, Tmax, 3Rvir', '-']]
    
    points_labels = [[SFhalos_stats[:,:,0] * (h/250)**3*fbary, r'$SFRD [Msun/yr/Mpc^3] /f$'],
                [SFhalos_stats[:,:,3] * (h/250)**3, r'$M_{halo} [Msun/Mpc^3]$'],
                [SFhalos_stats[:,:,1] * fbary, r'$SFR^{mean} [Msun/yr/f]$'],
                [SFhalos_stats[:,:,4] , r'$M_{halo}^{mean} [Msun]$'],
                [SFhalos_stats[:,:,2] * fbary, r'$SFR^{median} [Msun/yr/f]$'],
                [SFhalos_stats[:,:,5] , r'$M_{halo}^{median} [Msun]$']]

    f=figure(figsize=(12,10))
    for i in range(6):
        ax = f.add_subplot(3,2,i+1)
        iy, ilabel = points_labels[i]
        for j in (0, 9, 10, 11, 12):
            ax.plot(z_arr, iy[:,j], legs[j][1], label=legs[j][0])
        ax.set_xlabel('z')
        ax.set_ylabel(ilabel)
        ax.set_yscale('log')
        ax.set_xlim(0,10)
        if i==0:
            ax.set_ylim(1e-2, 1)
            legend(frameon=0, fontsize=12, loc=0, ncol=1)
            data_arr = genfromtxt('data.txt').T     
            z_data, SFRD_data, delSFRD_data, IR, author = data_arr
            z_data, SFRD_data, delSFRD_data, IR, author = data_arr[:,IR==1]
            ax.errorbar(z_data,  10**SFRD_data, 0.5*(10**(SFRD_data+delSFRD_data)-10**(SFRD_data-delSFRD_data)),fmt='+')
    plt.tight_layout()
    #show()
    ##savefig(plot_dir+'MassSFR_allcuts_M11.pdf')
    savefig(plot_dir+'MassSFR_allcuts_M11_ENVonly2.png')
    close()
            
############### minimizer to find f, sigma ########
#mdep=11
#eqer=0
#SFhalos_stats = load('SFhalos_stats_M%i_SigmaArr_cutM12.npy'%(mdep))
#SFRD = log10(SFhalos_stats[:,:,0]* (h/250)**3*fbary)[:,:-1]

#data_arr = genfromtxt('data.txt').T     
#z_data, SFRD_data, delSFRD_data, IR, author = data_arr
#z_data, SFRD_data, delSFRD_data, IR, author = data_arr[:,IR==1]
#if eqer:
    #delSFRD_data = ones(len(delSFRD_data))*mean(delSFRD_data)
#sigma_arr = arange(10,100,10)
#from scipy import optimize,interpolate
### models depends on f, z, sigma
#input_params = array([[[isigma, iz] for isigma in sigma_arr ]for iz in z_arr]).reshape(22*9,2)
#interp2d = interpolate.LinearNDInterpolator(input_params, SFRD.flatten())
#model_fcn = lambda sigma, z, f: log10(f) + interp2d(sigma, z)
#def chisq_fun(fsigma):
    #f,sigma = fsigma
    #imodel = model_fcn(sigma, z_data, f)
    #chisq = (imodel - SFRD_data)**2/delSFRD_data**2
    ##print chisq
    #return sum(chisq)

#for initg in arange(10,100,10):
    #x=optimize.minimize(chisq_fun, (0.2, initg), bounds=[(0,40), (0.001, 1.0)], method = 'Nelder-Mead')
    #print initg, x.x, x.fun/z_data.shape[0]

##errorbar(z_data,  SFRD_data, delSFRD_data, fmt='o' )
##plot(model_fcn(z_arr, x.x[0], x.x[1]))
##xlim(0,4)
##show()

#logf=log10(0.2)
#plot(z_arr, logf+interp2d(10, z_arr),'-',label='10km/s f=0.2')
#plot(z_arr, logf+interp2d(50, z_arr),'-',label='50km/s f=0.2')
#plot(z_arr, logf+interp2d(80, z_arr),'-',label='80km/s f=0.2')
#plot(z_arr, log10(x.x[0])+interp2d(x.x[1], z_arr),'--',label='f=%.2f, %.2fkm/s'%(x.x[0], x.x[1]))
##plot(z_arr, log10(0.25)+interp2d(14.7, z_arr),'--',label='f=0.25, 14.7km/s')

#errorbar(z_data,  SFRD_data, delSFRD_data, fmt='+' )
#legend(loc=0)
#xlabel('z')
#ylabel('SFRD')
#title('M^%.1f, Mvir_self <1e12'%(mdep/10.0))
#savefig('/Users/jia/Desktop/down_fitdataM%i_IRonly%s_cutM12.png'%(mdep,['','_equalerror'][eqer]))
#close()

################################################

if plot_Luminosity_distribution:
    ititle = 't_cool = t_acc'##'T_max = %s'%(printSci(Tcut))
    printSci = lambda a: '%ie%i'%(int(a/10**(int(log10(a)))),int(log10(a)) )
    #z_arr[3::3]=[ 4.48907674,  2.08699142,  1.3242301 ,  0.81241504,  0.59314311, 0.23081468,  0.02089778]
    data = [load(folder+ifile) for ifile in files[3::3]]
    Mdot_arr = [mdot_fcn(z_arr[i], data[i]) for i in range(len(data))]
    SFRD_arr = array([sum(imdot_arr) for imdot_arr in Mdot_arr])
    SFR_mean_arr = array([mean(imdot_arr) for imdot_arr in Mdot_arr])

    SFmass_mean_arr = array([mean(idata) for idata in data])
    #std_arr = array([std(idata) for idata in data])
    SFnum_halos = array([len(idata) for idata in data])
    SFmass_sum_arr = array([sum(idata) for idata in data])

    f = figure(figsize=(8,6))
    ax=f.add_subplot(111)
    for i in range(len(z_arr[3::3])):
        x,y=histogram(Mdot_arr[i],range=(0,100),bins=50)
        y=y[:-1]+0.5*(y[1]-y[0])
        ax.plot(y,x,label='%.2f'%(z_arr[3::3][i]))
    legend(title='z',ncol=2,frameon=0)
    xlabel('Mdot [Msun/yr]')
    ylabel('# of SF Halos')
    title(ititle)
    ax.set_yscale('log')
    plt.tight_layout()
    #plt.ticklabel_format(style='sci', axis='y')#, scilimits=(0,0)
    savefig(plot_dir+'Lum_dist_tacc.png')##'Lum_dist_T%s.png'%(printSci(Tcut))
    close()
    #data = [load(folder+ifile) for ifile in files]
    #Mdot_arr = [mdot_fcn(z_arr[i], data[i]) for i in range(len(data))]
    #SFRD_arr = array([sum(imdot_arr) for imdot_arr in Mdot_arr])
    #SFR_mean_arr = array([mean(imdot_arr) for imdot_arr in Mdot_arr])

    #SFmass_mean_arr = array([mean(idata) for idata in data])
    #SFnum_halos = array([len(idata) for idata in data])
    #SFmass_sum_arr = array([sum(idata) for idata in data])
    
    #data2 = genfromtxt('/Users/jia/Documents/smbhs/CosmicDownsizing/down_3vir_Tcut%s.out'%(printSci(Tcut))).T
    #ititle = '$r=min(3R_{vir}, T_{vir}*R_{vir}/T_{max}), T_{max}=%s$'%(printSci(Tcut))
    #z_arr2,Tot_halos,Vcut_halos,Tcut_halos,SF_halos,sum_arr,mean_arr,std_arr, mass_SF_sum = data2
    #iTot_halos=Tot_halos[argsort(z_arr2)[::-1]]
    
    

if plot_massSFR_tacc:

    f=figure(figsize=(10,8))
    
    for icut in (30,50,80):
        folder = '/Users/jia/Documents/smbhs/CosmicDownsizing/SFhalos_tacc_Vcut%i/'%(icut)
        files = os.listdir(folder)    
        ifn = 'points_tacc_Vcut%i.npy'%(icut)
        if not os.path.isfile(ifn):
            
            
            data = [load(folder+ifile) for ifile in files]
            Mdot_arr = [mdot_fcn(z_arr[i], data[i]) for i in range(len(data))]
            SFRD_arr = array([sum(imdot_arr) for imdot_arr in Mdot_arr])
            SFR_mean_arr = array([mean(imdot_arr) for imdot_arr in Mdot_arr])

            SFmass_mean_arr = array([mean(idata) for idata in data])
            SFnum_halos = array([len(idata) for idata in data])
            SFmass_sum_arr = array([sum(idata) for idata in data])
            
            data2 = genfromtxt('/Users/jia/Documents/smbhs/CosmicDownsizing/down_3vir_tacc_Vcut%i.out'%(icut)).T
            z_arr2,Tot_halos,Vcut_halos,Tcut_halos,SF_halos,sum_arr,mean_arr,std_arr, mass_SF_sum = data2
            iTot_halos=Tot_halos[argsort(z_arr2)[::-1]]
        else:
            SFRD_arr, SFR_mean_arr, SFmass_sum_arr, SFmass_mean_arr, iTot_halos, SFnum_halos= load(ifn)
       
        points_arr = [[SFRD_arr, 'SFRD [Mdot/yr/box]'],
                    [SFR_mean_arr, '<SFR> [Mdot/yr]'],
                    [SFmass_sum_arr, 'Mhalo_SF_all [Msun]'],
                    [SFmass_mean_arr, '<Mhalo_SF> [Msun]'],
                    [iTot_halos, '# tot. halos'],
                    [SFnum_halos, '# SF halos']]
        if not os.path.isfile(ifn):
            save('points_tacc_Vcut%i.npy'%(icut),array([ipoints_arr[0] for ipoints_arr in points_arr]))

        for i in range(6):
            ax = f.add_subplot(3,2,i+1)
            iy, iylabel = points_arr[i]
            ax.plot(z_arr, iy, '-', label='Vmin = %i'%(icut))
            #ax.set_title('t_cool = t_acc, Vmin = 50')
            ax.set_xlabel('z')
            ax.set_ylabel(iylabel)
            ax.set_yscale('log')
            ax.set_xlim(0,10)
            #ax.set_xscale('log')
            if i==0:
                legend(frameon=0, fontsize=12,loc=0)
    plt.tight_layout()
    #show()
    savefig(plot_dir+'MassSFR_r3vir_tacc_Vcut.png')
    close()
