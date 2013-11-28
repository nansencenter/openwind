# README  

OpenWind depends on Nansat (https://github.com/nansencenter/nansat). 

# Installation

- Add the top level folder openwind/ to your PYTHONPATH
- Add the subfolder openwind/openwind/ to your PATH (for command line usage)

# Command line usage:

```
./sar_wind.py -s SAR_image -w winddir -f figure_filename -n netCDF_filename -r resize_factor
```

- SAR_image is a file readable by Nansat, containing NRCS in VV polarisation

- winddir is a file readable by Nansat, containing wind direction (U10 and V10), or an integer indicating constant wind direction (0 from North, 90 from East etc)
If winddir the string 'archive', OpenWind looks for matching win dfields in local file archive (not yet implemented)
If winddir is not given, OpenWind tries to download NCEP GFS model wind for the time of the SAR image.

To see explanation of all features:
```
./sar_wind.py -h
```


# Python usage:
```
>>> from openwind import sar_wind

>>> s = sar_wind.SARwind(s, w) # To calculate SAR wind
```

- s is a filename (string) or a corresponding Nansat object containing NRCS in VV polarisation. The SAR Nansat object should preferrably be resized to pixels at least 200 m size, but presently only non-reprojected SAR images are supported.

- w is a filename (string) or a corresponding Nansat object containing wind direction (U10 and V10), or an integer indicating constant wind direction (0 from North, 90 from East etc). If winddir is not given, OpenWind tries to download NCEP GFS model wind for the time of the SAR image.


The following explicit usage allows e.g. (not yet!) reprojecting the SARWind object before calculation of wind:
```
>>> from openwind import sar_wind

>>> s = sar_wind.SARwind(s, winddir=None)

>>> s.reproject(some_domain)

>>> s.calculate_wind(winddir)
```

To plot the result and save as figure:
```
>>> plt = s.plot() # to plot SAR wind overlaid wind vectors.

>>> plt.savefig(filename, bbox_inches='tight', dpi=300) # Save to file
```

See code and comments therein for more features.

# Notes:
- GDAL might need to be compiled with the option --with-jasper to be able to read the downloaded NCEP GFS GRIB2-files
- Due to a bug in Nansat (#43), the following workaround is needed before resizing ASAR images:

```
>>> s_tmp = Nansat(ASAR_filename)

>>> s = Nansat(s_tmp.vrt.fileName)
```

# Acknowledgments

Thanks to the Royal Netherlands Meteorological Institute
(http://www.knmi.nl/scatterometer/cmod5/) for making the cmod5 software
available.

# References

Hersbach, H., Comparison of C-Band Scatterometer CMOD5.N Equivalent Neutral
Winds with ECMWF J. Atm. Oceanic Technol., 2010, 27, 721-736,
doi:10.1175/2009JTECHO698.1.

Portabella, M. and A.C.M. Stoffelen, On Scatterometer Ocean Stress J. Atm.
Oceanic Technol., 2009, 26, 2, 368-382, doi:10.1175/2008JTECHO578.1.

Verhoef, A., M. Portabella, A. Stoffelen and H. Hersbach, CMOD5.n - the CMOD5
GMF for neutral winds Document external project: 2008,
SAF/OSI/CDOP/KNMI/TEC/TN/165, EUMETSAT, 2008.

Hersbach, H., A. Stoffelen and S. de Haan, An Improved C-band scatterometer
ocean geophysical model function: CMOD5 J. Geophys. Res., 2007, 112, C03006,
doi:10.1029/2006JC003743.

Portabella, M., Wind field retrieval from satellite radar systems. Thesis:
University of Barcelona, 2002, Barcelona, Spain, 207p.

Haan, S. de and A. Stoffelen, CMOD5 Document external project: 2001,
SAF/OSI/KNMI/TEC/TN/140, EUMETSAT, 2001.
