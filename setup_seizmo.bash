#!/bin/bash
#SETUP_SEIZMO.BASH    Download & Unpackage Extra SEIZMO Packages
#
# This was made to get around LANL's Matlab proxy issues by using wget to
# do the downloading.  After this is run you will still want to run
# INSTALL_SEIZMO to install SEIZMO components to the static path.  If you
# just want SEIZMO in your dynamic path then you can skip INSTALL_SEIZMO
# after this script and use STARTUP_SEIZMO & SHUTDOWN_SEIZMO as needed.
#
# REQUIRES: WGET, UNZIP, CD, MKDIR, MV

# get jars for taup
cd mattaup/lib
wget http://www.seis.sc.edu/software/maven2/edu/sc/seis/TauP/2.1.1/TauP-2.1.1.jar
wget http://www.seis.sc.edu/software/maven2/edu/sc/seis/seisFile/1.5.1/seisFile-1.5.1.jar
cd ../..

# get njtbx
mkdir -p njtbx
cd njtbx
wget http://epsc.wustl.edu/~ggeuler/codes/m/seizmo/njToolbox-2.0.zip
wget http://epsc.wustl.edu/~ggeuler/codes/m/seizmo/toolsUI-4.0.49.jar
wget http://epsc.wustl.edu/~ggeuler/codes/m/seizmo/njTools-2.0.12_jre1.5.jar
unzip njToolbox-2.0.zip
cd ..

# get m_map
wget http://www.eos.ubc.ca/%7Erich/m_map1.4.zip
unzip m_map1.4.zip

# get irisws
cd ws
wget http://www.iris.edu/files/IRIS-WS/2/2.0.6/IRIS-WS-2.0.6.jar
wget http://www.iris.edu/files/irisFetch.m/2-0-6/irisFetch.m
cd ..

# get gshhg
mkdir -p gshhg
cd gshhg
wget ftp://ftp.soest.hawaii.edu/gshhg/gshhg-bin-2.2.4.zip
unzip gshhg-bin-2.2.4.zip
mv gshhg-bin-2.2.4.zip ..
cd ..

# get extras
wget http://epsc.wustl.edu/~ggeuler/codes/m/seizmo/seizmo_iris_sacpzdb.zip
wget http://epsc.wustl.edu/~ggeuler/codes/m/seizmo/seizmo_3d_models.zip
wget http://epsc.wustl.edu/~ggeuler/codes/m/seizmo/seizmo_mapping_features.zip
unzip seizmo_iris_sacpzdb.zip
unzip seizmo_3d_models.zip
unzip seizmo_mapping_features.zip

# get export_fig
mkdir -p export_fig
cd export_fig
wget https://github.com/ojwoodford/export_fig/archive/master.zip
mv master export_fig.zip
unzip export_fig.zip
mv export_fig.zip ..
cd ..
