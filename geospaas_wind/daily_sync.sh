#! /bin/bash
export PROJ_LIB=$HOME/Miniconda3-4.6.14-Linux-x86_64/envs/py3openwind/share/proj
source $HOME/Miniconda3-4.6.14-Linux-x86_64/bin/activate py3openwind
cd $HOME/project

# Sync S1A
$HOME/Miniconda3-4.6.14-Linux-x86_64/envs/py3openwind/bin/python $HOME/project/manage.py ingest_thredds_crawl http://nbstds.met.no/thredds/catalog/NBS/S1A/`date -d "yesterday 00:00" '+%Y/%m/%d'`/catalog.html

# Sync S1B
$HOME/Miniconda3-4.6.14-Linux-x86_64/envs/py3openwind/bin/python $HOME/project/manage.py ingest_thredds_crawl http://nbstds.met.no/thredds/catalog/NBS/S1B/`date -d "yesterday 00:00" '+%Y/%m/%d'`/catalog.html

# Sync Arome Arctic
$HOME/Miniconda3-4.6.14-Linux-x86_64/envs/py3openwind/bin/python $HOME/project/manage.py ingest_thredds_crawl https://thredds.met.no/thredds/catalog/aromearcticarchive/`date -d "yesterday 00:00" '+%Y/%m/%d'`/catalog.html --filename arome_arctic_vtk*

# Process S1 wind
$HOME/Miniconda3-4.6.14-Linux-x86_64/envs/py3openwind/bin/python $HOME/project/manage.py process_sentinel1_wind --date `date -d "yesterday 00:00" '+%Y-%m-%d'`
