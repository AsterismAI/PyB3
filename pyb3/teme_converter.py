import numpy as np
import astropy.coordinates
import astropy.units as u
import astropy.time

def GCRS_RADEC_to_TEME_RADEC( ra, dec ):
    L = np.ones( len(ra) ) * 1e20 * u.km
    ans = astropy.coordinates.spherical_to_cartesian( L, np.radians(dec), np.radians(ra) )
    gcrf = np.array(ans).T
    gcrs = astropy.coordinates.GCRS( obstime=date,
                                    x = gcrf[:,0] * u.km,
                                    y = gcrf[:,1] * u.km,
                                    z = gcrf[:,2] * u.km,
                                    representation_type='cartesian')
    # teme conversion
    teme = gcrs.transform_to( astropy.coordinates.TEME(obstime=date))
    temexyz = teme.cartesian.xyz.to_value(u.km).T
    temeradec = astropy.coordinates.cartesian_to_spherical( temexyz[:,0] * u.km,
                                                            temexyz[:,1] * u.km,
                                                            temexyz[:,2] * u.km )
    
    # 0: dist, 1: latitude, 2: longitude
    return temeradec[2].to_value(u.deg), temeradec[1].to_value(u.deg)
