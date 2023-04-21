# Copyright (C) 2016, 2020
# Joint Institute for VLBI ERIC, Dwingeloo, The Netherlands
#
# This library is free software; you can redistribute it and/or modify it
# under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
# License for more details.
#
# You should have received a copy of the GNU Library General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
#
from __future__ import print_function
import datetime
import math
import os
import sys
import time as tm
import re

try:
    # Python 2
    from StringIO import StringIO
except:
    # Python 3
    from io import StringIO

import numpy as np

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

from casavlbitools import key

os.environ['TZ'] = 'UTC'
tm.tzset()

class IdiHDU(pyfits.PrimaryHDU):
    @classmethod
    def match_header(cls, header):
        try:
            keyword = header.cards[0].keyword
        except:
            keyword = header.ascard[0].key
            pass
        return (keyword == 'SIMPLE' and 'GROUPS' in header and
                header['GROUPS'] == True and 'NAXIS' in header and
                header['NAXIS'] == 0)

pyfits.register_hdu(IdiHDU)

class IdiData:
    def __init__(self, idifiles):
        hdulist = pyfits.open(idifiles[0])

        # Fundamentals
        header = hdulist['ARRAY_GEOMETRY'].header
        self.obscode = header['OBSCODE']
        self.n_stokes = header['NO_STKD']
        self.stk_1 = header['STK_1']
        self.n_band = header['NO_BAND']
        self.n_chan = header['NO_CHAN']
        self.ref_freq = header['REF_FREQ']
        self.chan_bw = header['CHAN_BW']
        self.ref_pixl = header['REF_PIXL']
        self.rdate = header['RDATE']

        # Create antenna mapping
        tbhdu = hdulist['ARRAY_GEOMETRY']
        assert tbhdu.header['EXTVER'] == 1
        antennas = [s.strip() for s in tbhdu.data['ANNAME']]
        self.antenna_map = dict(zip(antennas, tbhdu.data['NOSTA']))

        # Create list of frequencies
        tbhdu = hdulist['FREQUENCY']
        assert tbhdu.data['FREQID'][0] == 1
        freqs = tbhdu.data['BANDFREQ'][0]
        if self.n_band == 1:
            freqs = [freqs]
            pass
        self.freqs = np.array([x + self.ref_freq for x in freqs])
        assert len(self.freqs) == self.n_band

        # Create source mapping
        tbhdu = hdulist['SOURCE']
        if 'ID_NO.' in tbhdu.columns.names:
            source_id_col = 'ID_NO.'
        else:
            source_id_col = 'SOURCE_ID'
            pass
        self.source_map = dict(zip(tbhdu.data['SOURCE'],
                                   tbhdu.data[source_id_col]))

        # Reference date in convenient form
        tupletime = tm.strptime(self.rdate, "%Y-%m-%d")
        self.reftime = tm.mktime(tupletime)

        # Create an index
        self.idx = []
        source_id = -1
        self.last_time = float("-inf")
        self.first_time = float("inf")
        for idifile in idifiles:
            hdulist = pyfits.open(idifile)
            tbhdu = hdulist['UV_DATA']
            if 'SOURCE' in tbhdu.columns.names:
                source_id_col = 'SOURCE'
            else:
                source_id_col = 'SOURCE_ID'
                pass
            for data in tbhdu.data:
                jd = data['DATE']
                time = (jd - 2440587.5 + data['TIME']) * 86400
                if time > self.last_time:
                    self.last_time = time
                    pass
                if time < self.first_time:
                    self.first_time = time
                    pass
                if data[source_id_col] != source_id:
                    source_id = data[source_id_col]
                    self.idx.append((time, source_id))
                    pass
                continue
            continue
        self.idx.sort()
        return

    def find_source(self, time):
        if time < self.first_time or time > self.last_time:
            return -1
        source_id = -1
        i = 0
        try:
            while time >= self.idx[i][0]:
                source_id = self.idx[i][1]
                i += 1
                continue
        except:
            pass
        return source_id


class TSysTable:
    def __init__(self):
        self.time = []
        self.time_interval = []
        self.source_id = []
        self.antenna_no = []
        self.array = []
        self.freqid = []
        self.tsys_1 = []
        self.tsys_2 = []
        self.tant = []
        return

class GainCurveTable:
    def __init__(self):
        self.antenna_no = []
        self.array = []
        self.freqid = []
        self.nterm = []
        self.y_typ = []
        self.gain = []
        self.sens_1 = []
        self.sens_2 = []
        return


def update_map(pols, spws, spwmap, index):
    idx = 0
    if not isinstance(index, (list, tuple)):
        index = [index]
        pass
    for labels in index:
        for label in labels.split('|'):
            pol = label[0]
            rng = label[1:].split(':')
            if pol != 'X':
                if not pol in pols:
                    pols.append(pol)
                    pass
                if len(rng) == 1:
                    rng.append(rng[0])
                    pass
                rng = [int(x) - 1 for x in rng]
                for spw in range(rng[0], rng[1] + 1):
                    if not spw in spws:
                        spws.append(spw)
                        pass
                    spwmap[(pol, spw)] = idx
                    continue
                pass
            continue
        idx += 1
        continue
    spws = sorted(spws)
    return

def find_antenna(keys, ignore):
    for key in keys[1:]:
        if not type(key[1]) is bool:
            continue
        if key[0] in ignore:
            continue
        return key[0]
    return None

def skip_values(infp):
    for line in infp:
        if line.startswith('!'):
            continue
        if line.strip().endswith('/'):
            break
        continue
    return

def get_timetuple(ts):
    # ts as string with these possible formats:
    # hh.hh      
    # hh:mm.mm   
    # hh:mm:ss.ss
    # NOTE: Regexs below will match any number of decimals on the last quantity (e.g. 19.8222222 and 19.8 both work)
    if re.match(r"[0-9]{2}\.[0-9]+", ts):
        # hh.hh 
        tm_hour = int(ts.split('.')[0])
        tm_min = math.modf(60*float(ts.split('.')[1]))
        tm_sec = int(60 * tm_min[0])
        tm_min = int(tm_min[1])
    elif re.match(r"[0-9]{2}:[0-9]{2}\.[0-9]+", ts):
        # hh:mm.mm 
        tm_hour = int(ts.split(':')[0])
        tm_min = math.modf(float(ts.split(':')[1]))
        tm_sec = int(60 * tm_min[0])
        tm_min = int(tm_min[1])
    elif re.match(r"[0-9]{2}:[0-9]{2}:[0-9]{2}$", ts):
        # hh:mm:ss
        tm_hour = int(ts.split(':')[0])
        tm_min = int(ts.split(':')[1])
        tm_sec = int(ts.split(':')[2])
    elif re.match(r"[0-9]{2}:[0-9]{2}:[0-9]{2}\.[0-9]+", ts):
        # hh:mm:ss.ss
        tm_hour = int(ts.split(':')[0])
        tm_min = int(ts.split(':')[1])
        tm_sec = float(ts.split(':')[2])
    return tm_hour, tm_min, tm_sec

def process_values(infp, keys, pols, idi, data):
    year = tm.gmtime(idi.reftime).tm_year
    antenna_name = find_antenna(keys[0], ['SRC/SYS'])
    if not antenna_name:
        print('Antenna missing from TSYS group')
        skip_values(infp)
        return
    try:
        antenna = idi.antenna_map[antenna_name]
    except:
        print('Antenna %s not present in FITS-IDI file' % antenna_name)
        skip_values(infp)
        return
    keys = dict(keys[0])
    spws = []
    spwmap = {}
    update_map(pols, spws, spwmap, keys['INDEX'])
    if 'INDEX2' in keys:
        update_map(pols, spws, spwmap, keys['INDEX2'])
        pass
    if len(spws) != idi.n_band:
        print('INDEX for antenna %s does not match FITS-IDI file'
              % antenna_name, file=sys.stderr)
        pass
    spws = range(idi.n_band)
    timeoff = 0
    if 'TIMEOFF' in keys:
        timeoff = float(keys['TIMEOFF'])
    for line in infp:
        if line.startswith('!'):
            continue
        fields = line.split()
        if len(fields) > 1:
            tm_year = year
            tm_yday = int(fields[0])
            # Get timestamp data depending on data format
            tm_hour, tm_min, tm_sec = get_timetuple(fields[1])
            t = "%dy%03dd%02dh%02dm%02ds" % \
                (tm_year, tm_yday, tm_hour, tm_min, tm_sec)
            t = tm.mktime(tm.strptime(t, "%Yy%jd%Hh%Mm%Ss"))
            days = (t + timeoff - idi.reftime) / 86400
            source = idi.find_source(t)
            values = fields[2:]
            tsys = {'R': [], 'L': []}
            for spw in spws:
                for pol in ['R', 'L']:
                    try:
                        value = float(values[spwmap[(pol, spw)]])
                        if value == 999.9:
                            value = float(-999.9)
                            pass
                    except:
                        value = float(-999.9)
                        pass
                    tsys[pol].append(value)
                    continue
                continue
            if source != -1:
                data.time.append(days)
                data.time_interval.append(0.0)
                data.source_id.append(source)
                data.antenna_no.append(antenna)
                data.array.append(1)
                data.freqid.append(1)
                data.tsys_1.append(tsys['R'])
                data.tsys_2.append(tsys['L'])
                data.tant.append(idi.n_band * [float(-999.9)])
                pass
            pass
        if line.strip().endswith('/'):
            break
        continue
    return

gain_keys = [ 'EQUAT', 'ALTAZ', 'ELEV', 'GCNRAO', 'TABLE', 'RCP', 'LCP' ]

def process_gc_values(infp, keys, pols, idi, data):
    antenna_name = find_antenna(keys[0], gain_keys)
    if not antenna_name:
        print('Antenna missing from GAIN group')
        skip_values(infp)
        return
    try:
        antenna = idi.antenna_map[antenna_name]
    except:
        print('Antenna %s not present in FITS-IDI file' % antenna_name)
        skip_values(infp)
        return
    keys = dict(keys[0])

    dpfu = {}
    try:
        dpfu['R'] = keys['DPFU'][0]
        dpfu['L'] = keys['DPFU'][1]
        pols.add('R')
        pols.add('L')
    except:
        dpfu['R'] = dpfu['L'] = keys['DPFU']
        pols.add('X')
        pass
    try:
        value = keys['POLY'][0]
    except:
        keys['POLY'] = [keys['POLY']]
        pass

    y_typ = 0
    if 'ELEV' in keys:
        y_typ = 1
    elif 'EQUAT' in keys:
        y_typ = 1
    elif 'ALTAZ' in keys:
        y_typ = 2
    else:
        print('Unknown gain curve type for antenna %s' % antenna_name)
        return

    poly = keys['POLY']
    data.antenna_no.append(antenna)
    data.array.append(1)
    data.freqid.append(1)
    data.y_typ.append(y_typ)
    data.nterm.append(len(poly))
    data.gain.append(poly)
    data.sens_1.append(dpfu['R'])
    data.sens_2.append(dpfu['L'])
    return

def append_tsys(antabfile, idifiles, replace=False):
    # Make sure we're dealing with a list
    if not isinstance(idifiles, list):
        idifiles = [idifiles]
        pass

    # Check if we already have a SYSTEM_TEMPERATURE table
    try:
        hdulist = pyfits.open(idifiles[0])
        hdu = hdulist['SYSTEM_TEMPERATURE']
        if not replace:
            raise RuntimeError('SYSTEM_TEMPERATURE table already present in FITS-IDI file')
        update = True
    except KeyError:
        update = False
        pass

    idi = IdiData(idifiles)
    data = TSysTable()

    pols = []
    keys = StringIO()
    fp = open(antabfile, 'r')
    for line in fp:
        if line.startswith('!'):
            continue
        keys.write(line)
        if line.strip().endswith('/'):
            keys.seek(0)
            try:
                tsys = key.read_keyfile(keys)
            except RuntimeError:
                print("\n", keys.getvalue(), file=sys.stderr)
                raise RuntimeError('error parsing ANTAB file')
            if tsys and tsys[0] and tsys[0][0][0] == 'TSYS':
                process_values(fp, tsys, pols, idi, data)
                pass
            keys = StringIO()
            continue
        continue

    time = data.time
    time_interval = data.time_interval
    source_id = data.source_id
    antenna_no = data.antenna_no
    array = data.array
    freqid = data.freqid
    tsys_1 = data.tsys_1
    tsys_2 = data.tsys_2
    tant = data.tant

    if len(pols) == 1 and pols[0] == 'L':
        tsys_1 = tsys_2
        pass

    cols = []
    col = pyfits.Column(name='TIME', format='1D', unit='DAYS', array=time)
    cols.append(col)
    col = pyfits.Column(name='TIME_INTERVAL', format='1E', unit='DAYS',
                        array=time_interval)
    cols.append(col)
    col = pyfits.Column(name='SOURCE_ID', format='1J', array=source_id)
    cols.append(col)
    col = pyfits.Column(name='ANTENNA_NO', format='1J', array=antenna_no)
    cols.append(col)
    col = pyfits.Column(name='ARRAY', format='1J', array=array)
    cols.append(col)
    col = pyfits.Column(name='FREQID', format='1J', array=freqid)
    cols.append(col)
    format = '%dE' % idi.n_band
    col = pyfits.Column(name='TSYS_1', format=format, unit='K', array=tsys_1)
    cols.append(col)
    col = pyfits.Column(name='TANT_1', format=format, unit='K', array=tant)
    cols.append(col)
    if len(pols) > 1:
        col = pyfits.Column(name='TSYS_2', format=format, unit='K',
                            array=tsys_2)
        cols.append(col)
        col = pyfits.Column(name='TANT_2', format=format, unit='K', array=tant)
        cols.append(col)
        pass
    coldefs = pyfits.ColDefs(cols)

    header = pyfits.Header()
    try:
        # Attempt to use the new interfaces available in PyFITS 3.1.x
        # and later.  The old interfaces are still available, but
        # generate warnings that say they're deprecated.
        header['EXTNAME'] = 'SYSTEM_TEMPERATURE'
        header['EXTVER'] = 1
        header['TABREV'] = 1
        header['OBSCODE'] = idi.obscode
        header['NO_STKD'] = idi.n_stokes
        header['STK_1'] = idi.stk_1
        header['NO_BAND'] = idi.n_band
        header['NO_CHAN'] = idi.n_chan
        header['REF_FREQ'] = idi.ref_freq
        header['CHAN_BW'] = idi.chan_bw
        header['REF_PIXL'] = idi.ref_pixl
        # Repeat the reference data even though the FITS-IDI standard
        # doesn't seem to require it.
        header['RDATE'] = idi.rdate
        header['NO_POL'] = len(pols)
        tbhdu = pyfits.BinTableHDU.from_columns(coldefs, header)
    except:
        # Fall back on the old interfaces available in PyFits 3.0.x
        # and earlier.
        header.update('EXTNAME', 'SYSTEM_TEMPERATURE')
        header.update('EXTVER', 1)
        header.update('TABREV', 1)
        header.update('OBSCODE', idi.obscode)
        header.update('NO_STKD', idi.n_stokes)
        header.update('STK_1', idi.stk_1)
        header.update('NO_BAND', idi.n_band)
        header.update('NO_CHAN', idi.n_chan)
        header.update('REF_FREQ', idi.ref_freq)
        header.update('CHAN_BW', idi.chan_bw)
        header.update('REF_PIXL', idi.ref_pixl)
        # Repeat the reference data even though the FITS-IDI standard doesn't
        # seem to require it.
        header.update('RDATE', idi.rdate)
        header.update('NO_POL', len(pols))
        rows = len(tsys_1)
        tbhdu = pyfits.core.new_table(coldefs, header, rows, False,
                                      'BinTableHDU')
        pass

    if update:
        pyfits.update(idifiles[0], tbhdu.data, tbhdu.header, 'SYSTEM_TEMPERATURE')
    else:
        hdulist = pyfits.open(idifiles[0], mode='append')
        hdulist.append(tbhdu)
        hdulist.close()
        pass
    return

def append_gc(antabfile, idifile, replace=False):
    # Check if we already have a GAIN_CURVE table
    try:
        hdulist = pyfits.open(idifile)
        hdu = hdulist['GAIN_CURVE']
        if not replace:
            raise RuntimeError('GAIN_CURVE table already present in FITS-IDI file')
        update = True
    except KeyError:
        update = False
        pass

    idi = IdiData([idifile])
    data = GainCurveTable()

    pols = set()
    keys = StringIO()
    fp = open(antabfile, 'r')
    for line in fp:
        if line.startswith('!'):
            continue
        keys.write(line)
        if line.strip().endswith('/'):
            keys.seek(0)
            gain = key.read_keyfile(keys)
            # Standard ANTAB
            if gain and gain[0] and gain[0][0][0] == 'GAIN':
                process_gc_values(fp, gain, pols, idi, data)
            elif gain and gain[0] and gain[0][0][0] == 'TSYS':
                skip_values(fp)
                pass
            keys = StringIO()
            continue
        continue

    n_tab = max(data.nterm)
    nrows = len(data.nterm)

    antenna_no = data.antenna_no
    array = data.array
    freqid = data.freqid
    _type = nrows * [idi.n_band * [2]]
    nterm = np.repeat(data.nterm, idi.n_band).reshape(nrows, idi.n_band)
    x_typ = nrows * [idi.n_band * [0]]
    y_typ = np.repeat(data.y_typ, idi.n_band).reshape(nrows, idi.n_band)
    x_val = nrows * [idi.n_band * [float('nan')]]
    y_val = nrows * [(n_tab * idi.n_band) * [float('nan')]]
    sens_1 = np.repeat(data.sens_1, idi.n_band).reshape(nrows, idi.n_band)
    if len(pols) > 1:
        sens_2 = np.repeat(data.sens_2, idi.n_band).reshape(nrows, idi.n_band)
        pass

    gain = np.zeros((nrows, n_tab))
    for row in range(nrows):
        for elem in range(len(data.gain[row])):
            gain[(row, elem)] = data.gain[row][elem]
            continue
        continue
    gain = np.repeat(gain, idi.n_band, axis=0)
    gain = gain.reshape(nrows, n_tab * idi.n_band)

    cols = []
    col = pyfits.Column(name='ANTENNA_NO', format='1J', array=antenna_no)
    cols.append(col)
    col = pyfits.Column(name='ARRAY', format='1J', array=array)
    cols.append(col)
    col = pyfits.Column(name='FREQID', format='1J', array=freqid)
    cols.append(col)
    format = '%dJ' % idi.n_band
    col = pyfits.Column(name='TYPE_1', format=format, array=_type)
    cols.append(col)
    col = pyfits.Column(name='NTERM_1', format=format, array=nterm)
    cols.append(col)
    col = pyfits.Column(name='X_TYP_1', format=format, array=x_typ)
    cols.append(col)
    col = pyfits.Column(name='Y_TYP_1', format=format, array=y_typ)
    cols.append(col)
    format = '%dE' % idi.n_band
    col = pyfits.Column(name='X_VAL_1', format=format, array=x_val)
    cols.append(col)
    format = '%dE' % (n_tab * idi.n_band)
    col = pyfits.Column(name='Y_VAL_1', format=format, array=y_val)
    cols.append(col)
    col = pyfits.Column(name='GAIN_1', format=format, array=gain)
    cols.append(col)
    format = '%dE' % idi.n_band
    col = pyfits.Column(name='SENS_1', format=format, unit='K/Jy',
                        array=sens_1)
    cols.append(col)
    if len(pols) > 1:
        format = '%dJ' % idi.n_band
        col = pyfits.Column(name='TYPE_2', format=format, array=_type)
        cols.append(col)
        col = pyfits.Column(name='NTERM_2', format=format, array=nterm)
        cols.append(col)
        col = pyfits.Column(name='X_TYP_2', format=format, array=x_typ)
        cols.append(col)
        col = pyfits.Column(name='Y_TYP_2', format=format, array=y_typ)
        cols.append(col)
        format = '%dE' % idi.n_band
        col = pyfits.Column(name='X_VAL_2', format=format, array=x_val)
        cols.append(col)
        format = '%dE' % (n_tab * idi.n_band)
        col = pyfits.Column(name='Y_VAL_2', format=format, array=y_val)
        cols.append(col)
        col = pyfits.Column(name='GAIN_2', format=format, array=gain)
        cols.append(col)
        format = '%dE' % idi.n_band
        col = pyfits.Column(name='SENS_2', format=format, unit='K/Jy',
                            array=sens_2)
        cols.append(col)
        pass
    coldefs = pyfits.ColDefs(cols)

    header = pyfits.Header()
    try:
        # Attempt to use the new interfaces available in PyFITS 3.1.x
        # and later.  The old interfaces are still available, but
        # generate warnings that say they're deprecated.
        header['EXTNAME'] = 'GAIN_CURVE'
        header['EXTVER'] = 1
        header['TABREV'] = 1
        header['OBSCODE'] = idi.obscode
        header['NO_STKD'] = idi.n_stokes
        header['STK_1'] = idi.stk_1
        header['NO_BAND'] = idi.n_band
        header['NO_CHAN'] = idi.n_chan
        header['REF_FREQ'] = idi.ref_freq
        header['CHAN_BW'] = idi.chan_bw
        header['REF_PIXL'] = idi.ref_pixl
        # Repeat the reference data even though the FITS-IDI standard
        # doesn't seem to require it.
        header['RDATE'] = idi.rdate
        header['NO_POL'] = 2 if len(pols) > 1 else 1
        header['NO_TABS'] = n_tab
        tbhdu = pyfits.BinTableHDU.from_columns(coldefs, header)
    except:
        # Fall back on the old interfaces available in PyFits 3.0.x
        # and earlier.
        header.update('EXTNAME', 'GAIN_CURVE')
        header.update('EXTVER', 1)
        header.update('TABREV', 1)
        header.update('OBSCODE', idi.obscode)
        header.update('NO_STKD', idi.n_stokes)
        header.update('STK_1', idi.stk_1)
        header.update('NO_BAND', idi.n_band)
        header.update('NO_CHAN', idi.n_chan)
        header.update('REF_FREQ', idi.ref_freq)
        header.update('CHAN_BW', idi.chan_bw)
        header.update('REF_PIXL', idi.ref_pixl)
        # Repeat the reference data even though the FITS-IDI standard doesn't
        # seem to require it.
        header.update('RDATE', idi.rdate)
        header.update('NO_POL', 2 if len(pols) > 1 else 1)
        header.update('NO_TABS', n_tab)
        rows = len(nterm)
        tbhdu = pyfits.core.new_table(coldefs, header, rows, False,
                                      'BinTableHDU')
        pass
    if update:
        pyfits.update(idifile, tbhdu.data, tbhdu.header, 'GAIN_CURVE')
    else:
        hdulist = pyfits.open(idifile, mode='append')
        hdulist.append(tbhdu)
        hdulist.close()
        pass
    return

def find_first_dobs(idifiles):
    first_dobs = datetime.datetime(datetime.MAXYEAR, 12, 31)
    for match in idifiles:
        hdulist = pyfits.open(match)
        dobs = hdulist['UV_DATA'].header['DATE-OBS']
        dobs = datetime.datetime.strptime(dobs, '%Y-%m-%d')
        if dobs < first_dobs:
            first_dobs = dobs
    return first_dobs

def convert_flags(infile, idifiles, outfp=sys.stdout, outfile=None):
    # Make sure we're dealing with a list
    if not isinstance(idifiles, list):
        idifiles = [idifiles]
        pass

    # Create new output file if requested to do so
    if outfile:
        outfp = open(outfile, 'w')
        pass

    first_date = find_first_dobs(idifiles)

    infp = open(infile, 'r')
    flags = key.read_keyfile(infp)
    infp.close()

    flags = [dict(x) for x in flags]

    for flag in flags:
        antenna = flag['ANT_NAME'].upper()
        if 'TIMERANG' in flag:
            timerange = flag['TIMERANG']
            if timerange == [0.0, 0.0, 0.0, 0.0, 400.0, 0.0, 0.0, 0.0]:
                timerange = ""
            elif isinstance(timerange, float) or isinstance(timerange, int):
                timerange = ""
            else:
            year = datetime.datetime(first_date.year, 1, 1)
            date1 = year + datetime.timedelta(timerange[0] - 1)
            date1 = date1.strftime("%Y/%m/%d")
            date2 = year + datetime.timedelta(timerange[4] - 1)
            date2 = date2.strftime("%Y/%m/%d")
            timerange = "%s/%02d:%02d:%02d~%s/%02d:%02d:%02d" % \
                (date1, timerange[1], timerange[2], timerange[3],
                 date2, timerange[5], timerange[6], timerange[7])
        else:
            timerange = ""
            
        if 'BIF' in flag:
            spw = "%d~%d" % (flag['BIF'] -1, flag['EIF'] - 1)
        else:
            spw = ""
            pass
        if 'BCHAN' in flag:
            spw = spw + ":%d~%d" % (flag['BCHAN']-1, flag['ECHAN']-1)
            pass
        cmd = "antenna='%s'" % antenna
        if timerange:
            cmd = cmd + " timerange='%s'" % timerange
            pass
        if spw:
            cmd = cmd + " spw='%s'" % spw
            pass
        outfp.write(cmd + "\n")
    return
