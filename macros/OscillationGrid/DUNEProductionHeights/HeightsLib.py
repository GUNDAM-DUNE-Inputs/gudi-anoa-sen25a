from numpy import *
import matplotlib.pyplot as plt
from scipy import interpolate, integrate
from time import time
import ROOT
import sys, os
import scienceplots # type: ignore

plt.style.use('science')
font = {'size'   : 16}
plt.rc('font', **font)
plt.rcParams['mathtext.fontset'] = 'stix'

path_to_CDF_folder = 'HondaHomestake'


'''Needed Functions'''

# Reads Honda tables and returns E, CDF⁻1(0.3), CDF⁻1(0.5) and CDF⁻1(0.7)
def Read_h(flavor, cosZ_indexbin, PLOT=False):
    #flavor = mu, mubar, e, ebar
    cosDict = {0:'[0.9, 1.0]', 1:'[0.8, 0.9]', 2:'[0.7, 0.8]', 3:'[0.6, 0.7]', 4:'[0.5, 0.6]', 5:'[0.4, 0.5]', 6:'[0.3, 0.4]', 7:'[0.2, 0.3]', 8:'[0.1, 0.2]',
               9:'[0.0, 0.1]', 10:'[-0.1, 0.0]', 11:'[-0.2, -0.1]', 12:'[-0.3, -0.2]', 13:'[-0.4, -0.3]', 14:'[-0.5, -0.4]', 15:'[-0.6, -0.5]', 16:'[-0.7, -0.6]', 
               17:'[-0.8, -0.7]', 18:'[-0.9, -0.8]', 19:'[-1.0, -0.9]'}

    mr = 101 # 21 or 101
    h_CDF = loadtxt('%s/hms-ally-aa-nu%s.d'%(path_to_CDF_folder, flavor), skiprows=1+102*cosZ_indexbin, max_rows=mr)

    # E[60] = 105.9 GeV
    E = h_CDF[:,0] # energy
    h_0p5 = h_CDF[:, 11]*1e-3 # median
    h_0p3 = h_CDF[:, 7]*1e-3 # CDF = 0.3
    h_0p7 = h_CDF[:, 15]*1e-3 # CDF = 0.7
    #h_0p99 = h_CDF[:, -1]*1e-3 # CDF = 0.7

    if PLOT==True:
        plt.figure(figsize=(5, 3.5))
        plt.plot(E, h_0p3, label=r'CDF=0.3')
        plt.plot(E, h_0p5, label=r'CDF=0.5')
        plt.plot(E, h_0p7, label=r'CDF=0.7')
        #plt.plot(E, h_0p99, label=r'CDF=0.99')
        plt.xlabel(r'$E$ [GeV]')
        plt.ylabel(r'$h$ [km]')
        plt.xscale('log')
        plt.legend()
        plt.title(r'$\cos\theta_z$ in $%s$'%cosDict[cosZ_indexbin])
        plt.show()
    return E, h_0p3, h_0p5, h_0p7
# E.g.: Read_h('mu', 10, True)

# Reads production height from specific quantile and returns E, h (km)
def Read_quant(flavor, cosZ_indexbin, Enu_indexbin, quant):
    Quant = {0.01:1, 0.05:2, 0.10:3, 0.15:4, 0.20:5, 0.25:6, 0.30:7, 0.35:8, 0.40:9, 0.45:10, 0.50:11, 0.55:12, 0.60:13, 0.65:14, 0.70:15, 0.75:16, 0.80:17, 0.85:18, 0.90:19, 0.95:20, 0.99:21}
    #cosDict = {0:'[0.9, 1.0]', 1:'[0.8, 0.9]', 2:'[0.7, 0.8]', 3:'[0.6, 0.7]', 4:'[0.5, 0.6]', 5:'[0.4, 0.5]', 6:'[0.3, 0.4]', 7:'[0.2, 0.3]', 8:'[0.1, 0.2]',
    #           9:'[0.0, 0.1]', 10:'[-0.1, 0.0]', 11:'[-0.2, -0.1]', 12:'[-0.3, -0.2]', 13:'[-0.4, -0.3]', 14:'[-0.5, -0.4]', 15:'[-0.6, -0.5]', 16:'[-0.7, -0.6]', 
    #           17:'[-0.8, -0.7]', 18:'[-0.9, -0.8]', 19:'[-1.0, -0.9]'}
    h_CDF = loadtxt('%s/hms-ally-aa-nu%s.d'%(path_to_CDF_folder,flavor), skiprows=1+102*cosZ_indexbin, max_rows=101)
    E = h_CDF[Enu_indexbin,0] # energy
    point = h_CDF[Enu_indexbin, Quant[quant]]*1e-3
    return E, point
# E.g.: Read_quant('mu', 10, 0, 0.5)

# Linear extrapolate to get b from ax+b=y
def ExtendToZero(Enu, h):
    a = (h[1] - h[0])/(Enu[1]-Enu[0])
    b = h[0] - a*Enu[0]
    return b

# Interpolates Honda points in a [E, cosZ] grid and saves functions for each percentile.
# In practice we will use that to calculate the height for a given percentile at an arbitrary [E, cosZ]
def InterpolateQuantile(Q):
    out_folder = 'Interpolated_CDF'
    os.system('mkdir -p %s'%out_folder)
    Points_mu = []
    Points_e = []
    Points_mub = []
    Points_eb = []
    E = Read_h('mu', 0)[0] # energy
    for c in range(20):
        points_mu = []
        points_e = []
        points_mub = []
        points_eb = []
        for e in range(101):
            points_mu.append(Read_quant('mu', c, e, Q)[1])
            points_e.append(Read_quant('e', c, e, Q)[1])
            points_mub.append(Read_quant('mubar', c, e, Q)[1])
            points_eb.append(Read_quant('ebar', c, e, Q)[1])
        b = ExtendToZero(E, points_mu)
        Points_mu.append([b]+points_mu)
        b = ExtendToZero(E, points_e)
        Points_e.append([b]+points_e)
        b = ExtendToZero(E, points_mub)
        Points_mub.append([b]+points_mub)
        b = ExtendToZero(E, points_eb)
        Points_eb.append([b]+points_eb)
    Points_mu = array(Points_mu).reshape((20, 102))
    Points_e = array(Points_e).reshape((20, 102))
    Points_mub = array(Points_mub).reshape((20, 102))
    Points_eb = array(Points_eb).reshape((20, 102))
    cosZ = linspace(-0.95, 0.95, 20)[::-1] # bin centers
    E = concatenate((array([0]), E), axis=0)
    points = [cosZ, E]
    f_mu = interpolate.RegularGridInterpolator(points, Points_mu, method='cubic', bounds_error=False, fill_value=None)
    f_e = interpolate.RegularGridInterpolator(points, Points_e, method='cubic', bounds_error=False, fill_value=None)
    f_mub = interpolate.RegularGridInterpolator(points, Points_mub, method='cubic', bounds_error=False, fill_value=None)
    f_eb = interpolate.RegularGridInterpolator(points, Points_eb, method='cubic', bounds_error=False, fill_value=None)
    save('%s/InterpQuant_mu_%.2f'%(out_folder, Q), f_mu)
    save('%s/InterpQuant_e_%.2f'%(out_folder, Q), f_e)
    save('%s/InterpQuant_mubar_%.2f'%(out_folder, Q), f_mub)
    save('%s/InterpQuant_ebar_%.2f'%(out_folder, Q), f_eb)
# E.g.: InterpolateQuantile(0.01)

# Makes a CDF in the same way as Honda but for arbitrary E, cosZ.
# Calls (if needed) and uses the output of InterpolateQuantile().
def MakeCDF(flavor, cosZ, Enu):
    # Check if interpolated functions exist. If not, make them:
    if os.path.exists(os.path.join(os.getcwd(), 'Interpolated_CDF', 'InterpQuant_mu_0.99.npy'))==False:
        print('Interpolated quantiles not found. Calculating...')
        for q in [0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99]:
            InterpolateQuantile(q)
            print('Calculated q=%.2f'%q)
    CDF = []
    for q in [0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99]:
        f = load('Interpolated_CDF/InterpQuant_%s_%.2f.npy'%(flavor, q), allow_pickle=True).item()
        CDF.append(f([cosZ, Enu])[0])
    return array(CDF)
# E.g.: MakeCDF('mu', -0.5, 1.)

# Makes interpolated CDF and PDF for a given [E, cosZ]. Returns h_Honda, CDF_Honda, h_CDF, CDF, h_PDF, PDF
def MakeDistributions(flavor, cosZ, Enu, hpoints=1000):
    CDF_Honda = [0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99] # percentiles
    h_Honda = MakeCDF(flavor, cosZ, Enu) # CDF points
    
    spl = interpolate.CubicSpline(CDF_Honda, h_Honda) # interpolation
    CDF = linspace(0, 1, hpoints) # CDF between 0 and 1 with hpoints
    h_CDF = spl(CDF) # calculate height for each CDF point

    # It may be that h_CDF has negative points. Let's correct that:
    if h_CDF[0]<0:
        CDF_Honda = [0]+CDF_Honda
        h_Honda = concatenate((array([0]), h_Honda), axis=0)
        spl = interpolate.CubicSpline(CDF_Honda, h_Honda)
        h_CDF = spl(CDF)
    
    #depricated
    #lp = where(h_CDF > 0, h_CDF, inf).argmin() #index of the least positive h
    #aux_h_CDF = linspace(0, h_CDF[lp], lp+1)
    #h_CDF = [aux_h_CDF[i] if hh<=0 else hh for i, hh in enumerate(h_CDF)]

    # Get PDF
    spl2 = interpolate.CubicSpline(h_CDF, CDF)
    h_PDF = linspace(h_CDF[0], h_CDF[-1], hpoints)
    PDF = spl2(h_PDF, nu=1)

    return h_Honda, CDF_Honda, h_CDF, CDF, h_PDF, PDF
# E.g.: MakeDistributions('mu', -0.5, 0.1)

# Finds index of a value closets to 'a' in a given array 'arr'
def FindEntry(a, arr):
    l = []
    for entry in arr:
        l.append(abs(a-entry))
    return array(l).argmin()

# Makes 1D histogram for a given flavor, cosZ and Enu
def MakeHistogram(flavor, cz, Enu, hbin_width=2, PLOT=False):

    # Get PDF
    h_Honda, CDF_Honda, h_CDF, CDF, h_PDF, PDF = MakeDistributions(flavor, cz, Enu)
    #_, _, _, _, h_PDF, PDF = MakeDistributions(flavor, cz, Enu)
    dh = h_PDF[1] - h_PDF[0]

    # Add 0 before
    if h_PDF[0] > 0:
        pre_start_x = arange(0, h_PDF[0]+dh, dh)
    else: pre_start_x=array([])
    pre_start_y = zeros(len(pre_start_x))

    # Add points until 70km
    post_end_x = arange(h_PDF[-1], 70+dh, dh)
    post_end_y = zeros(len(post_end_x))

    # Concatenate
    full_pdf = concatenate((pre_start_y, PDF, post_end_y), axis=None)
    full_h = concatenate((pre_start_x, h_PDF, post_end_x), axis=None)

    bins = arange(0, 70+hbin_width, hbin_width)
    int_i=[FindEntry(b, full_h) for b in bins]
    integrals=[sum(full_pdf[int_i[i]:int_i[i+1]])*dh for i in range(len(int_i)-1)]
    integrals=integrals/sum(integrals)
    histo=array(integrals)/hbin_width
    if abs(sum(histo)-1/hbin_width)>1e-3:print('Histogram doesn\'t sum to unity')

    if PLOT==True:
        plt.plot(h_CDF, CDF)
        plt.plot(h_Honda, CDF_Honda, 'or')
        plt.figure()
        plt.plot(h_PDF, PDF)
        plt.plot(pre_start_x, pre_start_y)
        plt.plot(post_end_x, post_end_y)
        plt.plot(full_h, full_pdf)
        plt.step(bins[1:], histo)
        plt.show()
    return histo
# E.g.: MakeHistogram('mu', -0.5, 1., hbin_width=2, PLOT=True)

# Makes a 3D histogram for a given Energy, Angle array, averaged over N subdivisions
def GetAverageHist(flavor, Energies, Angles, N, hbin_width=2, savename=''):
    print('Calculating nu_%s heights histogram for:'%flavor)
    print('Energy bins: ', Energies)
    print('Angle bins: ', Angles)
    print('with subdivision N = %d\n'%N)
    div=50
    Ntot=(len(Energies)-1)*(len(Angles)-1)*N*N
    counter=0
    HistogramAverage = []
    if N>1:
        f = load('InterpolatedFluxes/Phi_%s.npy'%(flavor), allow_pickle=True).item()
        for i in range(len(Energies)-1):
            Energy_subrange = linspace(Energies[i], Energies[i+1], N+1)[:-1]
            Energy_subrange = Energy_subrange+(Energy_subrange[1]-Energy_subrange[0])/2
            for z in range(len(Angles)-1):
                Angle_subrange = linspace(Angles[z], Angles[z+1], N+1)[:-1]
                Angle_subrange = Angle_subrange+(Angle_subrange[1]-Angle_subrange[0])/2
                Histograms= []
                WeightsSum = 0
                for ee in Energy_subrange:
                    for zz in Angle_subrange:
                        counter+=1
                        w = f([ee, zz])[0]
                        WeightsSum+=w
                        Histograms.append(w*MakeHistogram(flavor, zz, ee, hbin_width=hbin_width))
                        sys.stdout.write('\r')
                        sys.stdout.write("[%-50s] %.2f%%" % ('='*int(counter/Ntot*div), (counter*100/Ntot)))
                        sys.stdout.flush()
                HistogramAverage.append(sum(Histograms, axis=0)/WeightsSum)
    elif N==1:
        for i in range(len(Energies)-1):
            Energy_subrange = [(Energies[i]+Energies[i+1])/2]
            for z in range(len(Angles)-1):
                Angle_subrange = [(Angles[z]+Angles[z+1])/2]
                Histograms= []
                for ee in Energy_subrange:
                    for zz in Angle_subrange:
                        counter+=1
                        Histograms.append(MakeHistogram(flavor, zz, ee, hbin_width=hbin_width))
                        sys.stdout.write('\r')
                        sys.stdout.write("[%-50s] %.2f%%" % ('='*int(counter/Ntot*div), (counter*100/Ntot)))
                        sys.stdout.flush()
                HistogramAverage.append(sum(Histograms, axis=0))
    else:print('Error: N must be positive!')
    hbins = len(HistogramAverage[0])
    HistogramAverage=array(HistogramAverage).reshape(len(Energies)-1, len(Angles)-1, hbins)
    if savename=='': savename='Histogram_%s_E%d_cZ%d_N%d_hbin%.1f'%(flavor, len(Energies)-1, len(Angles)-1, N, hbin_width)
    save('HeightHistograms/%s'%savename, HistogramAverage)
    save('Binning/%s_BinsX'%savename, Energies)
    save('Binning/%s_BinsY'%savename, Angles)
    save('Binning/%s_BinsZ'%savename, arange(0, 70+hbin_width, hbin_width))
    print('\nDone calculating. Histogram saved to file %s.npy with dimensions '%savename, HistogramAverage.shape)
    print('\n')
    return HistogramAverage
# E.g.: GetAverageHist(flavor='mu', Energies=logspace(-1, 2, 31), Angles=linspace(-1, 1, 11), N=1, hbin_width=5)

# Returns a finer binning from a coarse binning.
def GetFineBinning(arr, N):
    if N==1:
        return arr
    FineBinning = []
    for i in range(len(arr)-1):
        FineBinning.append(linspace(arr[i], arr[i+1], N+1)[:-1])
    FineBinning=array(FineBinning).flatten()
    FineBinning=concatenate((FineBinning, array([arr[-1]])), axis=0)
    return FineBinning
# E.g.: GetFineBinning(linspace(0, 100, 11), 10)

# Makes a template to use the same binning when making an oscillogram.
def MakeOscillogramTemplate(EnergiesCoarse, AnglesCoarse, N_E, N_Z):
    EnergiesFine=GetFineBinning(EnergiesCoarse, N_E)
    AnglesFine=GetFineBinning(AnglesCoarse, N_Z)
    HistNameCoarse="OscillogramTemplate_Coarse"
    HistNameFine="OscillogramTemplate_Fine"
    HistLabels=";TrueNeutrinoEnergy;TrueCosineZ"
    hCoarse = ROOT.TH2D(HistNameCoarse, HistLabels, len(EnergiesCoarse)-1, EnergiesCoarse, len(AnglesCoarse)-1, AnglesCoarse)
    hFine = ROOT.TH2D(HistNameFine, HistLabels, len(EnergiesFine)-1, EnergiesFine, len(AnglesFine)-1, AnglesFine)
    OutputFile = ROOT.TFile("OscillogramTemplate_new.root", "RECREATE")
    hCoarse.Write()
    hFine.Write()
    OutputFile.Close()
# E.g.: MakeOscillogramTemplate(logspace(-1, 2, 31), linspace(-1, 1, 11), N_E=10, N_Z=10)

# Saves to ROOT file in the format needed bu CUDAProb3. Note that although a specific flavor is needed in filename, the root file will contain all flavors.
def SaveRootTH3D(filename):
    Energies, Angles, Heights = load('Binning/%s_BinsX.npy'%filename), load('Binning/%s_BinsY.npy'%filename), load('Binning/%s_BinsZ.npy'%filename)
    Energies, Angles, Heights = array(Energies, dtype=double), array(Angles, dtype=double), array(Heights, dtype=double)
    flavors = ['mu', 'e', 'mubar', 'ebar']
    flavor_index=filename.index('_E')
    Histograms = []
    for f in flavors:
        flavorfile='Histogram_%s'%f+filename[flavor_index:]
        Histograms.append(load('HeightHistograms/%s.npy'%flavorfile))

    # Fill TH3D
    H4 = []
    for j, f in enumerate(flavors):
        HistName="ProductionHeight_nu%s"%f
        HistLabels=";True Neutrino Energy (GeV);True Cosine Zenith;Production Height (km)"
        h3 = ROOT.TH3D(HistName, HistLabels, len(Energies)-1, Energies, len(Angles)-1, Angles, len(Heights)-1, Heights)
        for enu in range(len(Energies)-1):
            for cz in range(len(Angles)-1):
                for h in range(len(Heights)-1):
                    h3.Fill(Energies[enu], Angles[cz], Heights[h], Histograms[j][enu, cz, h])#change
        H4.append(h3)

    # Open file and write histogram
    savename = filename[flavor_index:]
    OutputFile = ROOT.TFile("ProdHeightDist%s.root"%savename, "RECREATE")
    for h in H4:
        h.Write()
    OutputFile.Close()
# E.g.: SaveRootTH3D('Histogram_mu_E30_cZ10_N1_hbin5.0')

# Reads ROOT TH3D amd returns the probabilities histogram, binsX, binsY and binsZ.
def ReadRootTH3D(filename, flavor):
    file = ROOT.TFile.Open("%s.root"%filename)
    h3 = file.Get("ProductionHeight_nu%s"%flavor)
    
    # Get bins:
    NbinsX, NbinsY, NbinsZ = h3.GetNbinsX(), h3.GetNbinsY(), h3.GetNbinsZ()
    binsX, binsY, binsZ = [], [], []
    for x in range(NbinsX):
        binsX.append(h3.GetXaxis().GetXbins().At(x))
    for y in range(NbinsY):
        binsY.append(h3.GetYaxis().GetXbins().At(y))
    for z in range(NbinsZ):
        binsZ.append(h3.GetZaxis().GetXbins().At(z))
    
    # Get probabilities:
    H3 = []
    for x in range(NbinsX):
        H2 = []
        for y in range(NbinsY):
            for z in range(NbinsZ):
                H2.append(h3.GetBinContent(x+1, y+1, z+1)) # Bins are indexed by 1 and not 0. Took a long - LONG - time to figure that out.
        H3.append(array(H2).reshape(NbinsY, NbinsZ))
    H3=array(H3)
    
    file.Close()

    print('Histogram has shape: ', H3.shape)
    return H3, binsX, binsY, binsZ
# E.g.: ReadRootTH3D('ProdHeightDist_E30_cZ10_N1_hbin5.0', 'mu')

''' Plotting '''

def PlotProjCZ(filename):
    if 'bar' in filename:
        if 'mu' in filename: flavor='mubar'
        else: flavor='ebar'
    elif 'mu' in filename: flavor='mu'
    else: flavor='e'
    LatexDict = {'mu':r'\nu_\mu', 'mubar':r'\bar{\nu}_\mu', 'e':r'\nu_e', 'ebar':r'\bar{\nu}_e'}
    binsY, binsZ = load('Binning/%s_BinsY.npy'%filename), load('Binning/%s_BinsZ.npy'%filename)
    H = sum(load('HeightHistograms/%s.npy'%filename), axis=0)
    X, Y = meshgrid(binsY, binsZ)
    fig, ax = plt.subplots(figsize=(11, 8))
    c0 = ax.pcolor(X, Y, H.T, cmap='cividis')
    ax.set_xlabel(r'$\cos\theta_z$')
    ax.set_ylabel(r'$h$ [km]')
    ax.set_title(r'$%s$ - Projection in $\cos\theta_z \times h$'%LatexDict[flavor])
    cbar0 = fig.colorbar(c0, ax=ax)
    cbar0.set_label(r'$\sum_{E_\nu}\quad p(h|E_\nu$, $\theta_z)$ [km$^{-1}$]', rotation=270, labelpad=27)
    plt.tight_layout()
    plt.savefig('Plots/ProjectionCZ_%s.png'%filename)
    plt.savefig('Plots_pdf/ProjectionCZ_%s.pdf'%filename)
    plt.show()
# E.g.: PlotProjCZ('Histogram_mu_E30_cZ10_N1_hbin5.0')

def PlotProjE(filename):
    if 'bar' in filename:
        if 'mu' in filename: flavor='mubar'
        else: flavor='ebar'
    elif 'mu' in filename: flavor='mu'
    else: flavor='e'
    LatexDict = {'mu':r'\nu_\mu', 'mubar':r'\bar{\nu}_\mu', 'e':r'\nu_e', 'ebar':r'\bar{\nu}_e'}
    binsX, binsZ = load('Binning/%s_BinsX.npy'%filename), load('Binning/%s_BinsZ.npy'%filename)
    H = sum(load('HeightHistograms/%s.npy'%filename), axis=1)
    X, Y = meshgrid(binsX, binsZ)
    
    fig, ax = plt.subplots(figsize=(11, 8))
    c1 = ax.pcolor(X, Y, H.T, cmap='cividis', shading='flat')
    ax.set_xlabel(r'$E_\nu$ [GeV]')
    ax.set_xscale('log')
    ax.set_ylabel(r'$h$ [km]')
    cbar1 = fig.colorbar(c1, ax=ax)
    cbar1.set_label(r'$\sum_{\theta_z}\quad p(h|E_\nu$, $\theta_z)$ [km$^{-1}$]', rotation=270, labelpad=27)

    plt.title(r'$%s$ - Projection in $E_\nu \times h$'%LatexDict[flavor])
    plt.tight_layout()
    plt.savefig('Plots/ProjectionE_%s.png'%filename)
    plt.savefig('Plots_pdf/ProjectionE_%s.pdf'%filename)
    plt.show()
# E.g.: PlotProjE('Histogram_mu_E30_cZ10_N1_hbin5.0')

def PlotSliceCZ(filename, cosZ):
    if 'bar' in filename:
        if 'mu' in filename: flavor='mubar'
        else: flavor='ebar'
    elif 'mu' in filename: flavor='mu'
    else: flavor='e'
    LatexDict = {'mu':r'\nu_\mu', 'mubar':r'\bar{\nu}_\mu', 'e':r'\nu_e', 'ebar':r'\bar{\nu}_e'}
    
    binsX, binsY, binsZ = load('Binning/%s_BinsX.npy'%filename), load('Binning/%s_BinsY.npy'%filename), load('Binning/%s_BinsZ.npy'%filename)
    H_Honda = load('HeightHistograms/%s.npy'%filename)
    czb=FindEntry(cosZ, binsY)
    X, Y = meshgrid(binsX, binsZ)
    
    fig, ax = plt.subplots(figsize=(11, 8))
    
    c1 = ax.pcolor(X, Y, H_Honda[:, czb, :].T, cmap='cividis', shading='flat')
    ax.set_xlabel(r'$E_\nu$ [GeV]')
    ax.set_xscale('log')
    ax.set_ylabel(r'$h$ [km]')
    cbar1 = fig.colorbar(c1, ax=ax)
    cbar1.set_label(r'$p(h|E_\nu$, $\theta_z)$ [km$^{-1}$]', rotation=270, labelpad=27)

    plt.title(r'$%s$ - $\cos\theta_z$ in [%.3f, %.3f]'%(LatexDict[flavor], binsY[czb], binsY[czb+1]))
    plt.tight_layout()
    plt.savefig('Plots/SliceCZ_%.2f_%s.png'%(cosZ, filename))
    plt.savefig('Plots_pdf/SliceCZ_%.2f_%s.pdf'%(cosZ, filename))
    plt.show()
# E.g.: PlotSliceCZ('Histogram_mu_E30_cZ10_N1_hbin5.0', 0.)

def PlotSliceE(filename, E):
    if 'bar' in filename:
        if 'mu' in filename: flavor='mubar'
        else: flavor='ebar'
    elif 'mu' in filename: flavor='mu'
    else: flavor='e'
    LatexDict = {'mu':r'\nu_\mu', 'mubar':r'\bar{\nu}_\mu', 'e':r'\nu_e', 'ebar':r'\bar{\nu}_e'}
    
    binsX, binsY, binsZ = binsX, binsY, binsZ = load('Binning/%s_BinsX.npy'%filename), load('Binning/%s_BinsY.npy'%filename), load('Binning/%s_BinsZ.npy'%filename)

    H = load('HeightHistograms/%s.npy'%filename)
    eb=FindEntry(E, binsX)
    X, Y = meshgrid(binsY, binsZ)
    
    fig, ax = plt.subplots(figsize=(11, 8))
    
    c1 = ax.pcolor(X, Y, H[eb, :, :].T, cmap='cividis', shading='flat')
    ax.set_xlabel(r'$\cos\theta_z$')
    ax.set_ylabel(r'$h$ [km]')
    cbar1 = fig.colorbar(c1, ax=ax)
    cbar1.set_label(r'$p(h|E_\nu$, $\theta_z)$ [km$^{-1}$]', rotation=270, labelpad=27)
    
    plt.title(r'$%s$ - $E_\nu$ in [%.2f, %.2f]'%(LatexDict[flavor], binsX[eb], binsX[eb+1]))
    plt.tight_layout()
    plt.savefig('Plots/SliceE_%.2f_%s.png'%(E, filename))
    plt.savefig('Plots_pdf/SliceE_%.2f_%s.pdf'%(E, filename))
    plt.show()
# E.g.: PlotSliceE('Histogram_mu_E30_cZ10_N1_hbin5.0', 1.)