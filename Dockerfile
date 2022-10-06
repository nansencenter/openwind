FROM nansencenter/nansat

RUN apt update && apt install -y unzip
# Install additional python packages
# Including ASF python API for downloading SAR scenes <asf_search> and Copernicus <cdsapi
RUN conda install jupyterlab matplotlib \
 && conda install -c conda-forge folium cdsapi asf_search \
 && pip install pythesint==1.6.6 pyresample dateparser==1.1.1 pytest ipytest
# Update meadata vocabularies
RUN python -c 'import pythesint as pti; pti.update_all_vocabularies()'
# Install Sentinel1denoise package
RUN pip install https://github.com/nansencenter/sentinel1denoised/archive/v1.3.1.tar.gz

WORKDIR /src
# Add GDAL related env vars
ARG GDAL_ENABLE_DEPRECATED_DRIVER_DODS='YES'
