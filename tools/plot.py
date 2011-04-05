import healpy
import glob
import numpy as np
import matplotlib.pyplot as plt
from planck import ps

tag = '70GHz 9x9'
short = '70'
import os
plotsfolder = 'plots/'+short+'/'
try:
    os.mkdir(plotsfolder)
except:
    pass

spar = ['S%d' % i for i in range(1,6+1)]
for m in spar[-1:]:
    print(m)
    a=healpy.read_map('out/'+short+'/binmap_%s.fits' % m)
    a = ps.smooth(a, 30)
    healpy.mollview(healpy.ud_grade(a,64)*1e3,min=-.1,max=.1,title='%s %s nside 64' % (m, tag),unit='mK')
    plt.savefig(plotsfolder + '%s_%s.png' % (m, tag.replace(' ','_')))

print('RCOND')
r=healpy.read_map('out/'+short+'/rcondmap.fits')
healpy.mollview(np.log10(r),min=-5,max=0,title='LOG10(RCOND) ' + tag,xsize=800)
plt.savefig(plotsfolder + 'rcond_%s.png' % (tag.replace(' ','_')))

ext={'Q':1,'U':2}
wmap = 'wmap/wmap_band_iqumap_r9_7yr_V_v4.fits'
for m in ['Q','U']:
    print(m)
    a=healpy.read_map('out/'+short+'/binmap_%s.fits' % m)
    a = ps.smooth(a, 30)
    healpy.mollview(healpy.ud_grade(a,64)*1e3,min=-.1,max=.1,title='IQUSS %s %s nside 64' % (m, tag),unit='mK')
    plt.savefig(plotsfolder + 'IQUSS_%s_%s.png' % (m, tag.replace(' ','_')))

    dpc = glob.glob('/project/projectdirs/planck/data/mission/DPC_maps/dx6/lfi/LFI_*_%s.fits' % (tag.split()[1]))[0]
    a=healpy.read_map(dpc, ext[m])
    a = ps.smooth(a, 30)
    healpy.mollview(healpy.ud_grade(a,64)*1e3,min=-.1,max=.1,title='DPC %s %s nside 64' % (m, tag),unit='mK')
    plt.savefig(plotsfolder + 'DPC_%s_%s.png' % (m, tag.replace(' ','_')))

    a=healpy.read_map(wmap, ext[m])
    a = ps.smooth(a, 30)
    healpy.mollview(healpy.ud_grade(a,64),min=-.1,max=.1,title='WMAP %s %s nside 64' % (m, tag),unit='mK')
    plt.savefig(plotsfolder + 'WMAP_%s_%s.png' % (m, tag.replace(' ','_')))
