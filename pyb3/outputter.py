import os
import glob
import sys
import astropy.coordinates
import astropy.units as u
import astropy.time
from datetime import datetime, timedelta
import json
import numpy as np

# ----------------------------------------- EPOCH -----------------------------------------
ds50epoch = astropy.time.Time( datetime.strptime( '1949-12-31T00:00:00', "%Y-%m-%dT%H:%M:%S") )

def ds50ToATime( flt ):
    try: return astropy.time.Time( ds50epoch.jd + flt , format='jd' )
    except Exception as e: 
        print('ds50ToATime: cannot parse input: {}'.format(str(e)))
        return None

def ds50ToDateTime( flt ):
    '''
    given the days-since-1950 float, return a datetime
     "YYDDDHHMMSS.SSS"
     '''
    return ds50ToATime( flt ).datetime


def B3_float_field( val, left, right ):
    '''
    given a float, we need a structured output that puts the decimal at a specific column
    total field length = left + right
    further, B3's deal with negatives with a character mapping to keep the field-lengths constant (see below)
    '''
    if val < 0: neg = True
    else: neg = False
    l,r = str(np.abs(val)).split('.')
    l = l[-left:].rjust(left,'0')
    r = r[:right].ljust(right,'0')
    if not neg: return l + r
    if l[0] == '0' : return '-' + l[1:] + r
    if l[0] == '1' : return 'J' + l[1:] + r
    if l[0] == '2' : return 'K' + l[1:] + r
    if l[0] == '3' : return 'L' + l[1:] + r
    if l[0] == '4' : return 'M' + l[1:] + r
    if l[0] == '5' : return 'N' + l[1:] + r
    if l[0] == '6' : return 'O' + l[1:] + r
    if l[0] == '7' : return 'P' + l[1:] + r
    if l[0] == '8' : return 'Q' + l[1:] + r
    if l[0] == '9' : return 'R' + l[1:] + r

def makeDate( datetm ): return datetm.strftime('%y%j%H%M%S%f')[:14]

def makeCommon( obdata, datetm=None, classification='U' ):
    ts      = list(' '*76)  # init the string
    ts[0]   = 'U'
    ts[1:6] = '{:05d}'.format( int(obdata['XA_OBS_SATNUM'] ) )
    ts[6:9] = '{:03d}'.format( int(obdata['XA_OBS_SENNUM'] ) )
    # we can pass in a datetime, or use what's in the struct (mostly this will be used for non A.S. data)
    if datetm == None:
        timedatetime = ds50ToDateTime( obdata['XA_OBS_DS50UTC'])
        ts[9:23] = makeDate( timedatetime )
    else:
        ts[9:23] = makeDate( datetm )
    return ts

def makeEl( el ): return B3_float_field( el, 2, 4)

def makeRA( dec ):
    '''
    Each hour is 360∘/24=15∘. 
    Each minute of time is 15∘/60=15′, i.e. 15 arcminutes. 
    Each second of time is 15′/60=15′′, i.e. 15 arcseconds.
    '''
    dec = (dec + 360) % 360  # lock it to [0,360] and the values below should not overflow
    hours  = int( dec / 15 )
    frac   = dec - (15 * hours)
    minut  = int( frac / 0.25 )
    frac -= 0.25 * minut 
    secs   =  frac * 86400./360.
    secsS  = '{:04.1f}'.format( secs ).replace('.','')
    frac -= secs * 0.25/60
    return "{:02d}{:02d}{}".format( hours, minut, secsS )

def makeRange( rangeval ):
    '''
    return the range field and the exponent location
    '''
    exp = int(np.log10( rangeval ))
    expval = exp-1  # this is the exponent mapping (valid ranges are 99.99999 to 9,999,999, so 0=10^1)
    if expval < 0 or expval > 5: raise
    return str(rangeval).replace('.','').ljust(7,'0')[0:7], str(expval)

def fortran9p3( flt ): 
    '''
    this is used for type 9 EFG sensor locations
    '''
    if flt < 0: neg = True
    else : neg = False
    l,r = str(np.abs(flt)).split('.')
    l = l[-6:].rjust(6,'0')
    r = r[-3:].ljust(3,'0')
    if neg : return '-' + l[1:] + r
    return '+' + l[1:] + r


# In[4]:


# ------------------------------------  TYPE 1 ---------------------------------------
def maketype1( obdata, datetm=None ):
    ts = makeCommon( obdata, datetm=datetm )
    ts[23:29] = makeEl( obdata['XA_OBS_ELORDEC'])
    ts[30:37] = B3_float_field( obdata['XA_OBS_AZORRA'], 3,4)
    ts[74] = '1'
    ts[75] = '0'
    return ''.join(ts)

# ------------------------------------  TYPE 2 ---------------------------------------
def maketype2( data, datetm=None ):
    '''
    4 - Elevation, azimuth, range, range rate, elevation rate, azimuth rate, rate acceleration
    types 2,3,4 are a subset of this...

    '''
    ts = makeCommon(data, datetm=datetm) 
    ts[23:29] = makeEl( data['XA_OBS_ELORDEC'])
    ts[30:37] = B3_float_field( data['XA_OBS_AZORRA'],3,4)
    rgval, rgexp = makeRange(data['XA_OBS_RANGE'])  # this carves out the exponent...
    ts[38:45] = rgval
    ts[45]    = rgexp
    ts[74] = '2'
    ts[75] = '0'
    return ''.join(ts)


# ------------------------------------  TYPE 3 ---------------------------------------
def maketype3( data, datetm=None ):
    '''
    4 - Elevation, azimuth, range, range rate, elevation rate, azimuth rate, rate acceleration
    types 2,3,4 are a subset of this...

    '''
    ts = makeCommon(data, datetm=datetm) 
    ts[23:29] = makeEl( data['XA_OBS_ELORDEC'])
    ts[30:37] = B3_float_field( data['XA_OBS_AZORRA'],3,4)
    rgval, rgexp = makeRange(data['XA_OBS_RANGE'])  # this carves out the exponent...
    ts[38:45] = rgval
    ts[45]    = rgexp
    ts[47:54] = B3_float_field( data['XA_OBS_RANGERATE'],2,5)  
    ts[74] = '3'
    ts[75] = '0'
    return ''.join(ts)

# ------------------------------------  TYPE 4 ---------------------------------------
def maketype4( data, datetm=None ):
    '''
    4 - Elevation, azimuth, range, range rate, elevation rate, azimuth rate, rate acceleration
    types 2,3,4 are a subset of this...

    '''
    ts = makeCommon(data, datetm=datetm) 
    ts[23:29] = makeEl( data['XA_OBS_ELORDEC'])
    ts[30:37] = B3_float_field( data['XA_OBS_AZORRA'],3,4)
    rgval, rgexp = makeRange(data['XA_OBS_RANGE'])  # this carves out the exponent...
    ts[38:45] = rgval
    ts[45]    = rgexp
    ts[47:54] = B3_float_field( data['XA_OBS_RANGERATE'],2,5)  
    ts[55:60] = B3_float_field( data['XA_OBS_ELRATE'],1,4) 
    ts[61:66] = B3_float_field( data['XA_OBS_AZRATE'],1,4) 
    ts[67:72] = B3_float_field( data['XA_OBS_RANGEACCEL'],1,4)
    ts[74] = '4'
    ts[75] = '0'
    return ''.join(ts)

# ------------------------------------  TYPE 5 ---------------------------------------
def maketype5( data, datetm=None ):
    ts = makeCommon( data, datetm=datetm )
    ts[23:29] = B3_float_field( data['XA_OBS_ELORDEC'], 2, 4)
    ts[30:37] = makeRA( data['XA_OBS_AZORRA'])
    ts[74] = '5'
    ts[75] = '0'
    return ''.join(ts)

# ------------------------------------  TYPE 6 ---------------------------------------
def maketype6( data, datetm=None ):
    ts = makeCommon( data, datetm=datetm )
    rgval, rgexp = makeRange(data['XA_OBS_RANGE'])  # this carves out the exponent...
    ts[38:45] = rgval
    ts[45]    = rgexp
    ts[74] = '6'
    ts[75] = '0'
    return ''.join(ts)

# ------------------------------------  TYPE 9 ---------------------------------------
def maketype9( data, datetm=None ):
    ts = makeCommon( data, datetm=datetm )
    ts[23:29] = makeEl( data['XA_OBS_ELORDEC'])
    ts[30:37] = makeRA( data['XA_OBS_AZORRA'])
    ts[38:45] = '0000000'
    ts[46:55] = fortran9p3( data['XA_OBS_POSX'] ) 
    ts[55:64] = fortran9p3( data['XA_OBS_POSY'] ) 
    ts[64:73] = fortran9p3( data['XA_OBS_POSZ'] ) 
    ts[74] = '9'
    ts[75] = '0'
    return ''.join(ts)

def b3_dispatcher( data ):
    if data['XA_OBS_OBSTYPE'] == 1: return maketype1( data  )
    if data['XA_OBS_OBSTYPE'] == 2: return maketype2( data  )
    if data['XA_OBS_OBSTYPE'] == 3: return maketype3( data  )
    if data['XA_OBS_OBSTYPE'] == 4: return maketype4( data  )
    if data['XA_OBS_OBSTYPE'] == 5: return maketype5( data  )
    if data['XA_OBS_OBSTYPE'] == 6: return maketype6( data  )
    if data['XA_OBS_OBSTYPE'] == 9: return maketype9( data  )


if __name__ == "__main__":
    Q = B3s[0].toAstrostdDict()
    b3_dispatcher(Q)

    for B in B3s:
        print(B.origline)
        print( b3_dispatcher( B.toAstrostdDict() ) )
