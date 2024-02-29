# The platform flag is required for making it work on Apple Silicon build
# for the nansat image
FROM --platform=linux/amd64 nansencenter/nansat
LABEL maintainer="Artem Moiseev <artem.moiseev@nersc.no>"

RUN apt update && apt install -y unzip
# Install additional python packages
# Including ASF python API for downloading SAR scenes <asf_search> and Copernicus <cdsapi
RUN pip install jupyterlab matplotlib ipykernel \
 && pip install folium cdsapi asf_search \
 && pip install pyresample dateparser==1.1.1 pytest ipytest pythesint==1.6.6 pandas "xarray[complete]"
# Install Sentinel1denoise package
RUN pip install https://github.com/nansencenter/sentinel1denoised/archive/v1.3.1.tar.gz
# Update meadata vocabularies
RUN python -c 'import pythesint as pti; pti.update_all_vocabularies()'

WORKDIR /src
# Add GDAL related env vars
ARG GDAL_ENABLE_DEPRECATED_DRIVER_DODS='YES'
    