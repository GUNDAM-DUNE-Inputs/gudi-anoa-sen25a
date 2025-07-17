from numpy import *
import matplotlib.pyplot as plt
from scipy import interpolate, integrate
from time import time
import ROOT
import sys, os

path_to_Honda_folder = 'HondaHomestake/'
out_folder = 'InterpolatedFluxes'

def ExtendToZero(Enu, h):
    a = (h[1] - h[0])/(Enu[1]-Enu[0])
    b = h[0] - a*Enu[0]
    return array([b])

# Function to implement log-interpolation (E)
def log_interp1d(x, y, kind='linear'):
    logx = log10(x)
    logy = log10(y)
    lin_interp = interpolate.interp1d(logx, logy, kind=kind, bounds_error=None, fill_value="extrapolate")
    log_interp = lambda z: power(10.0, lin_interp(log10(z)))
    return log_interp

# Area preserving interpolation (theta)
def interp(avgs):
    X = arange(-1, 1.1, .1)#range(len(avgs)+1)
    F = [0]
    for i in arange(1, len(avgs)+1, 1):
        F.append(F[i - 1] + avgs[i - 1] * (X[i] - X[i - 1]))
    spl = interpolate.CubicSpline(X, F)
    xnew=cos(linspace(0, 180, 181)*pi/180)
    ynew = spl(xnew, nu=1)
    return(ynew[::-1])

#Honda Tables
solmax = 'hms-nu-20-01-n3650.d'
solmin = 'solmin_file_name'

def InterpolateFluxE(solcycle, SubGeV=False):
    if SubGeV==True: mr=21
    else: mr = 101

    flux_mu = []
    flux_mub = []
    flux_e = []
    flux_eb = []
    for i in arange(0, 20, 1):                               
        var = loadtxt(path_to_Honda_folder+solcycle, unpack=True, skiprows=2+103*(19-i), max_rows=mr) #starting from [-1,-0.9]
        flux_mu.append(concatenate((ExtendToZero(var[0], var[1]), var[1]), axis=0))
        flux_mub.append(concatenate((ExtendToZero(var[0], var[2]), var[2]), axis=0))
        flux_e.append(concatenate((ExtendToZero(var[0], var[3]), var[3]), axis=0))
        flux_eb.append(concatenate((ExtendToZero(var[0], var[4]), var[4]), axis=0))
    x = concatenate((array([1e-3]), var[0]), axis=0)

    if SubGeV==True: xnew = linspace(0.1, 1, 901)
    else: 
        xnew = concatenate((linspace(0.001, 100, 100001)[:-1], logspace(2, 4, 101)), axis=0)
        
    
    os.system('mkdir -p InterpolatedFluxes')
    orig_stdout = sys.stdout
    if solcycle==solmin:
        if SubGeV==True: f = open('InterpolatedFluxes/Honda_Interpolated_noOsc_solmin_SubGeV.txt', 'w')
        else: f = open('InterpolatedFluxes/Honda_Interpolated_noOsc_solmin.txt', 'w')
    else:
        if SubGeV==True: f = open('InterpolatedFluxes/Honda_Interpolated_noOsc_SubGeV.txt', 'w')
        else: f = open('InterpolatedFluxes/Honda_Interpolated_noOsc.txt', 'w')
    sys.stdout = f

    for i in range(20):
        print('\t\t cos(theta_z) in ['+str((i-10)/10)+','+str((i-9)/10)+']')
        print('E(GeV)  phi_e       phi_m       phi_eb      phi_mb\t (m^2 sec sr GeV)^-1')
        fe = log_interp1d(x, flux_e[i])
        fm = log_interp1d(x, flux_mu[i])
        feb = log_interp1d(x, flux_eb[i])
        fmb = log_interp1d(x, flux_mub[i])

        for e in xnew:
            print(str(format(e, '.7f'))+'    '+ str(format(fe(e), '.7f'))+'    '+ str(format(fm(e), '.7f'))+'    '+ str(format(feb(e), '.7f'))+'    '+ str(format(fmb(e), '.7f')))
        print('\n')

    sys.stdout = orig_stdout
    f.close()

def InterpolateFluxTheta(solcycle, SubGeV=False):
    if SubGeV == True:
        mr=901
        if solcycle=='solmin':file= "InterpolatedFluxes/Honda_Interpolated_noOsc_solmin_SubGeV.txt"
        else: file= "InterpolatedFluxes/Honda_Interpolated_noOsc_SubGeV.txt"
    else: 
        mr = len(concatenate((linspace(0.001, 100, 100001)[:-1], logspace(2, 4, 101)), axis=0))
        if solcycle=='solmin':file= "InterpolatedFluxes/Honda_Interpolated_noOsc_solmin.txt"
        else: file= "InterpolatedFluxes/Honda_Interpolated_noOsc.txt"

    Honda = []
    #print('Loaded file: %s'%file)
    for bin in range(20):
        Honda.append(loadtxt(file, skiprows=2+((mr+4)*bin), unpack=True, max_rows=mr))
    Honda = array(Honda)

    phi_e = []
    phi_m = []
    phi_eb = []
    phi_mb = []

    for e in range(mr):
        phi_e.append(interp(Honda[:, 1, e]))
        phi_m.append(interp(Honda[:, 2, e]))
        phi_eb.append(interp(Honda[:, 3, e]))
        phi_mb.append(interp(Honda[:, 4, e]))
    phi_e = array(phi_e)
    phi_m = array(phi_m)
    phi_eb = array(phi_eb)
    phi_mb = array(phi_mb)

    if SubGeV == True:
        if solcycle=='solmin':
            savetxt('InterpolatedFluxes/phi_e_Interpolated_noOsc_solmin_SubGeV.txt', phi_e)
            savetxt('InterpolatedFluxes/phi_eb_Interpolated_noOsc_solmin_SubGeV.txt', phi_eb)
            savetxt('InterpolatedFluxes/phi_m_Interpolated_noOsc_solmin_SubGeV.txt', phi_m)
            savetxt('InterpolatedFluxes/phi_mb_Interpolated_noOsc_solmin_SubGeV.txt', phi_mb)
        else:
            savetxt('InterpolatedFluxes/phi_e_Interpolated_noOsc_SubGeV.txt', phi_e)
            savetxt('InterpolatedFluxes/phi_eb_Interpolated_noOsc_SubGeV.txt', phi_eb)
            savetxt('InterpolatedFluxes/phi_m_Interpolated_noOsc_SubGeV.txt', phi_m)
            savetxt('InterpolatedFluxes/phi_mb_Interpolated_noOsc_SubGeV.txt', phi_mb)

    else:
        if solcycle=='solmin':
            savetxt('InterpolatedFluxes/phi_e_Interpolated_noOsc_solmin.txt', phi_e)
            savetxt('InterpolatedFluxes/phi_eb_Interpolated_noOsc_solmin.txt', phi_eb)
            savetxt('InterpolatedFluxes/phi_m_Interpolated_noOsc_solmin.txt', phi_m)
            savetxt('InterpolatedFluxes/phi_mb_Interpolated_noOsc_solmin.txt', phi_mb)
        else:
            savetxt('InterpolatedFluxes/phi_e_Interpolated_noOsc.txt', phi_e)
            savetxt('InterpolatedFluxes/phi_eb_Interpolated_noOsc.txt', phi_eb)
            savetxt('InterpolatedFluxes/phi_m_Interpolated_noOsc.txt', phi_m)
            savetxt('InterpolatedFluxes/phi_mb_Interpolated_noOsc.txt', phi_mb)

print('Starting Flux Interpolation...')
#InterpolateFluxE(solmax)
InterpolateFluxTheta('solmax')

#Prepare for regular grid interpolator
flavors=['e', 'mu', 'ebar', 'mubar']
flv_label={'e':'e', 'mu':'m', 'ebar':'eb', 'mubar':'mb'}
E = concatenate((linspace(0.001, 100, 100001)[:-1], logspace(2, 4, 101)), axis=0)
cosZ = cos(linspace(0, 180, 181)*pi/180)
points = [E, cosZ[::-1]]
for flv in flavors:
    Phi=loadtxt('%s/phi_%s_Interpolated_noOsc.txt'%(out_folder, flv_label[flv]))
    f = interpolate.RegularGridInterpolator(points, Phi, method='linear', bounds_error=False, fill_value=None)
    save('%s/Phi_%s'%(out_folder, flv), f)

os.system('rm InterpolatedFluxes/*.txt')
print('Concluded.')