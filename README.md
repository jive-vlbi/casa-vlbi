# casa-vlbi
Scripts to assist VLBI data reduction in CASA

The functionality of the append_tsys.py, gc.py and flag.py scripts is
also available in Python modules.  The following interfaces are
available:

```python
casvlbitools.fitsidi.append_tsys(andtabfile, idifiles)
```

Append Tsys measurements from the file named by `antabfile` (which
should be in AIPS ANTAB format or VLBA VLOG format) to the FITS-IDI
file named by `idifiles`.  If an observation is split into multiple
FITS-IDI files, `idifile` should be a list of file names.

```python
casvlbitools.fitsidi.append_gc(andtabfile, idifile)
```

Append gain curves and sensitivity information from the file named by
`antabfile` (which should be in AIPS ANTAB format) to the FITS-IDI
file named by `idifile`.

```python
casavlbitools.casa.convert_gaincurve(antabfile, gc, min_elevation=0.0, max_elevation=90.0)
```

Create a gaincurve table with name `gc` from gain curves from the file
named by `antablfile` (which should be in AIPS ANTAB format or VLBA
gains format).  The gain curves are sampled between `min_elevation`
and `max_elevation`, which should be specified in degrees, and then
refitted.

```python
casavlbitools.fitsidi.convert_flags(infile, idifiles, outfp=sys.stdout, outfile=None)
```

Convert the flag file named by `infile` (which should be in AIPS UVFLG
format) into a format that can be used by CASA.  The name of the
FITS-IDI files corresponding to the observsation should be passed in
`idifiles`.  By default the converted output is written to stdout.
Output is written to a file named `outpufile` instead if that
parameter is provided.
