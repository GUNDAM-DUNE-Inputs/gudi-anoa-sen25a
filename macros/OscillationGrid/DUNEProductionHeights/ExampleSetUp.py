from HeightsLib import *

Energies = logspace(-1, 2, 31)
Angles = linspace(-1, 1, 11)

print('Calculating Height Histograms...')

GetAverageHist('mu', Energies, Angles, N=1, hbin_width=5)
GetAverageHist('e', Energies, Angles, N=1, hbin_width=5)
GetAverageHist('mubar', Energies, Angles, N=1, hbin_width=5)
GetAverageHist('ebar', Energies, Angles, N=1, hbin_width=5)

print('Concluded.')

print('Saving to root file...')

SaveRootTH3D('Histogram_mu_E30_cZ10_N1_hbin5.0')

MakeOscillogramTemplate(Energies, Angles, 1, 1)

print('Done.')

print('Plotting example:')

PlotProjCZ('Histogram_mu_E30_cZ10_N1_hbin5.0')
#PlotProjE('Histogram_mu_E30_cZ10_N1_hbin5.0')
#PlotSliceCZ('Histogram_mu_E30_cZ10_N1_hbin5.0', 0.)
#PlotSliceE('Histogram_mu_E30_cZ10_N1_hbin5.0', 1.)

print('Done.')