# README  

OpenWind depends on Nansat (https://github.com/nansencenter/nansat). 

# Installation

- Add the top level folder openwind/ to your PYTHONPATH
- Add the subfolder openwind/openwind/ to your PATH (for command line usage)

# Command line usage:

```
./sar_wind.py -s SAR_filename -w winddir -f figure_filename -n netCDF_filename -p pixelsize
```

- SAR_filename is a file readable by Nansat, containing Normalised Radar Cross Section (NRCS)

- winddir is either:
  - a file readable by Nansat, containing wind direction (U10 and V10)
  - an integer indicating constant wind direction (0 from North, 90 from East etc)
  - name (string) of a mapper with functionality to find a wind file (online or on local disk) matching the SAR image time [DEFAULT: 'ncep_wind_online']

- pixelsize is given in meters (500 m by default). Use pixelsize='fullres' for no resizing (usually not recommended, unless image is small or cropped).

To see explanation of all features:
```
./sar_wind.py -h
```


# Python usage:
```
>>> from openwind import SARWind

>>> s = SARWind(SAR_image, winddir, pixelsize)
```

See above (command line usage) for specification of the input parameters.
SAR image may here also be a Nansat SAR object (containing NRCS, incidence angle and look direction)


To plot the result and save as figure:
```
>>> s.plot() # to plot SAR wind overlaid wind vectors.

>>> s.plot(filename, show=False) # Save to file
```


# Notes:
- GDAL might need to be compiled with the option --with-jasper to be able to read the downloaded NCEP GFS GRIB2-files

# Acknowledgments

OpenWind is developed with support from the Norwegian Space Centre.

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
