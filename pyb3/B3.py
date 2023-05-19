# ###############################################################################
# MIT License
# 
# Copyright (c) 2023 Kerry Wood
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ###############################################################################

import astropy.time
import json
from datetime import datetime
import outputter


# ----------------------------------------- EPOCH -----------------------------------------
ds50epoch = astropy.time.Time( datetime.strptime( '1949-12-31T00:00:00', "%Y-%m-%dT%H:%M:%S") )


# these fields are pulled out of the AstroStandards code and pushed into JSON
# spit out of the A.S. parser
as_fields = json.loads('["XA_OBS_ASTAT", "XA_OBS_AZORRA", "XA_OBS_AZRATE", "XA_OBS_DS50UTC", "XA_OBS_ELORDEC", "XA_OBS_ELRATE", "XA_OBS_OBSTYPE", "XA_OBS_POSX", "XA_OBS_POSY", "XA_OBS_POSZ", "XA_OBS_RANGE", "XA_OBS_RANGEACCEL", "XA_OBS_RANGERATE", "XA_OBS_SATNUM", "XA_OBS_SECCLASS", "XA_OBS_SENNUM", "XA_OBS_SIGMAEL1", "XA_OBS_SIGMAEL10", "XA_OBS_SIGMAEL11", "XA_OBS_SIGMAEL12", "XA_OBS_SIGMAEL13", "XA_OBS_SIGMAEL14", "XA_OBS_SIGMAEL15", "XA_OBS_SIGMAEL16", "XA_OBS_SIGMAEL17", "XA_OBS_SIGMAEL18", "XA_OBS_SIGMAEL19", "XA_OBS_SIGMAEL2", "XA_OBS_SIGMAEL20", "XA_OBS_SIGMAEL21", "XA_OBS_SIGMAEL3", "XA_OBS_SIGMAEL4", "XA_OBS_SIGMAEL5", "XA_OBS_SIGMAEL6", "XA_OBS_SIGMAEL7", "XA_OBS_SIGMAEL8", "XA_OBS_SIGMAEL9", "XA_OBS_SIGMATYPE", "XA_OBS_SITETAG", "XA_OBS_SPADOCTAG", "XA_OBS_TRACKIND", "XA_OBS_VELX", "XA_OBS_VELY", "XA_OBS_VELZ"]')

# classification map (astrostandards stores the classification as a number)
classmap = {'U':1.0, 'C':2.0, 'S':3.0}

# float map
charmap = {
    'J' : '-1',
    'K' : '-2',
    'L' : '-3',
    'M' : '-4',
    'N' : '-5',
    'O' : '-6',
    'P' : '-7',
    'Q' : '-8',
    'R' : '-9'
}

# -----------------------------------------------------------------------------------------------------
class B3:
    jd1950 = 2433281.5
    origline = None
    def __init__( self, L=None ):
        if L == None: 
            raise Exception('B3 was not passed in input line')
            return
        self.default = -1.0
        self.origline = L
        self.parse( self.origline )

    def setdate( self ):
        self.isot = "{:04d} {:03d} {:02d} {:02d} {:02d} {}".format(
                            self.year,
                            self.doy,
                            self.hour,
                            self.minute,
                            self.second, 
                            str(self.millis).ljust(6,'0'))
        self.datetime = datetime.strptime( self.isot, '%Y %j %H %M %S %f' )
        return self.datetime

    def parse( self, L ):
        # common fields
        self.classification = L[0]
        self.satid   = int( L[1:6] )
        self.sensid  = int( L[6:9] )
        self.year    = int( L[9:11] )
        if self.year < 50 : self.year += 2000
        else: self.year += 1900
        self.doy     = int( L[11:14] )
        self.hour    = int( L[14:16] )
        self.minute  = int( L[16:18] )
        self.second  = int( L[18:20] )
        self.millis  = int( L[20:23])
        self.alldate = L[9:23]

        self.obstype = int( L[74] )

        # elevation or declination
        if self.obstype in set([1,2,3,4,5,8,9]):
            # deal with the char mapping
            poschar = L[23]
            prefix  = L[23]
            rest    = L[24:29]
            if poschar in charmap: prefix = charmap[poschar]
            self.eledec = float( prefix + rest ) / 10000.

        # azimuth or right ascension
        self.azra = 0.
        if self.obstype in set([1,2,3,4,8]): self.azra = float( L[30:37] ) / 10000
        if self.obstype in set([5,9]): 
            self.azra = float(L[30:32]) + float(L[32:34])/60. + float(L[34:37])/36000
            self.azra *= 360./24.

        try: self.rgexp = float( L[45] )
        except: self.rgexp = self.default

        try: self.range = (float(L[38:45]) / 100000) * (10 ** self.rgexp)
        except: self.range = self.default

        # other types of obs
        if self.obstype in set([8,9]):
            self.rngrate = 0.
            self.ecfx    = float(L[46:55]) / 1000.
            self.ecfy    = float(L[55:64]) / 1000.
            self.ecfz    = float(L[64:73]) / 1000.
        else: 
            self.ecfx = self.ecfy = self.ecfz = 0.

        try: self.rngrate = float(L[47:54]) / 100000.
        except: self.rngrate = self.default

        try: self.elrate = float(L[55:60]) / 10000
        except: self.elrate = self.default

        try: self.azrate = float(L[61:66]) / 10000
        except: self.azrate = self.default

        try: self.rangeacc = float( L[67:72] ) / 10000 
        except: self.rangeacc = self.default

        self.equinox = None
        if L[75] == ' ' or L[75] == '0' : self.equinox = 'TEME'
        if L[75] == '1' : self.equinox = 'YEAR'
        if L[75] == '2' : self.equinox = 'J2K'
        if L[75] == '3' : self.equinox = '1950'


        try: self.track_position = int( L[76] )
        except: self.track_position = None

        # sometimes, astat value is included
        try: self.astat = int( L[79] )
        except: self.astat = None

        # sometimes the sensor applies a different satellite number
        try: self.site_tag    = int(L[80:84])
        except: self.site_tag = None

        # SPADOC applied tag number
        try: self.spadoc_tag    = int(L[85:90])
        except: self.spadoc_tag = None

        self.datetime = self.setdate()


    def todict ( self ):
        td = {'obstype' : self.obstype,
              'classification' : self.classification,
              'sensid'  : self.sensid,
              'satid'   : self.satid,
              'date'    : self.datetime.isoformat(),
              'el_dec'  : self.eledec,
              'az_ra'   : self.azra,
              'range'   : self.range,
              'range_rate' : self.rngrate,
              'ecfx'    : self.ecfx,
              'ecfy'    : self.ecfy,
              'ecfz'    : self.ecfz,
              'elrate'  : self.elrate,
              'azrate'  : self.azrate,
              'rangeacc' : self.rangeacc,
              'equinox' : self.equinox,
              'site_tag' : self.site_tag,
              'spadoc_tag' : self.spadoc_tag,
              'track_position' : self.track_position
              }
        return td

    def toAstrostdDict( self ):
        rv = {}
        for v in as_fields: rv[v] = 0.0
        rv['XA_OBS_SECCLASS'] = classmap[ self.classification ]
        rv['XA_OBS_DS50UTC'] = astropy.time.Time( self.datetime ).jd - self.jd1950
        rv['XA_OBS_SATNUM']  = self.satid
        rv['XA_OBS_OBSTYPE'] = self.obstype
        rv['XA_OBS_SENNUM']  = self.sensid
        if self.obstype in set([1,2,3,4,5,8,9]): rv['XA_OBS_ELORDEC'] = self.eledec
        if self.obstype in set([1,2,3,4,5,8,9] ): rv['XA_OBS_AZORRA'] = self.azra
        if self.obstype in set([2,3,4,6]) : rv['XA_OBS_RANGE'] = self.range
        if self.obstype in set([3,4,0]) : rv['XA_OBS_RANGERATE'] = self.rngrate
        if self.obstype == 4:
            rv['XA_OBS_AZRATE'] = self.azrate
            rv['XA_OBS_ELRATE'] = self.elrate
        if self.obstype in set([8,9]):
            rv['XA_OBS_POSX'] = self.ecfx
            rv['XA_OBS_POSY'] = self.ecfy
            rv['XA_OBS_POSZ'] = self.ecfz
        if self.site_tag :        rv['XA_OBS_SITETAG'] = self.site_tag
        if self.spadoc_tag :      rv['XA_OBS_SPADOCTAG'] = self.spadoc_tag
        if self.track_position:   rv['XA_OBS_TRACKIND'] = self.track_position
        if self.astat:            rv['XA_OBS_ASTAT'] = self.astat
        return rv

    def __repr__( self ): 
        return str( self.todict() )

    def prettyp( self ):
        import json
        return json.dumps( self.todict(), indent=5 )

    def toB3( self ):
        return outputter.b3_dispatcher( self.toAstrostdDict() )



# =====================================================================================================
if __name__ == "__main__" : 
    pass

