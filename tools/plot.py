import healpy
import glob
import numpy as np
import matplotlib.pyplot as plt

tag = '70GHz 18-23'
for m in ['S1','S2']:
    print(m)
    a=healpy.read_map('binmap_%s.fits' % m)
    healpy.mollview(healpy.ud_grade(a,64)*1e3,min=-.1,max=.1,title='%s %s nside 64' % (m, tag),unit='mK')
    plt.savefig('%s_%s.png' % (m, tag.replace(' ','_')))

print('RCOND')
r=healpy.read_map('rcondmap.fits')
healpy.mollview(np.log10(r),min=-5,max=0,title='LOG10(RCOND) ' + tag,xsize=800)
plt.savefig('rcond_%s.png' % (tag.replace(' ','_')))

ext={'Q':1,'U':2}
wmap = 'wmap_band_iqumap_r9_7yr_V_v4.fits'
for m in ['Q','U']:
    print(m)
    a=healpy.read_map('binmap_%s.fits' % m)
    healpy.mollview(healpy.ud_grade(a,64)*1e3,min=-.1,max=.1,title='IQUSS %s %s nside 64' % (m, tag),unit='mK')
    plt.savefig('IQUSS_%s_%s.png' % (m, tag.replace(' ','_')))

    dpc = glob.glob('/u/zonca/dx6/LFI_*_%s.fits' % (tag.split()[1]))[0]
    a=healpy.read_map(dpc, ext[m])
    healpy.mollview(healpy.ud_grade(a,64)*1e3,min=-.1,max=.1,title='DPC %s %s nside 64' % (m, tag),unit='mK')
    plt.savefig('DPC_%s_%s.png' % (m, tag.replace(' ','_')))

    a=healpy.read_map(wmap, ext[m])
    healpy.mollview(healpy.ud_grade(a,64),min=-.1,max=.1,title='WMAP %s %s nside 64' % (m, tag),unit='mK')
    plt.savefig('WMAP_%s_%s.png' % (m, tag.replace(' ','_')))

