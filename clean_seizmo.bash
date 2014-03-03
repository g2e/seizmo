#!/bin/bash
#CLEAN_SEIZMO.BASH    Deletes Extra SEIZMO Packages
#
# This is meant to clean up a SEIZMO install to a "vanilla" state.  It is
# useful if you are having issues with downloaded packages.  Note that you
# will need to follow this with SETUP_SEIZMO.BASH and/or INSTALL_SEIZMO.
#
# REQUIRES: RM

# remove extras (both downloaded & expanded)
rm -r *.zip export_fig gshhg m_map njtbx
rm -r mattaup/lib/TauP-*.jar mattaup/lib/seisFile-*.jar
rm -r ws/IRIS-WS-*.jar ws/irisFetch.m event/globalcmt_full.mat
rm -r event/globalcmt_quick.mat response/sacpzdb.mat mapping/*.gz
rm -r mapping/*.mat models/block.desc models/CN*txt models/crust1*
rm -r models/*.mat models/*.bin

# now clean up messed up zips
rm sacpzdb.mat
rm CRUST*.mat CUB2.mat DZ04.mat HMSL06* MITP08.mat PRI05.mat
rm S20RTS.mat SAW24B16.mat SB4L18.* TX2006.mat TX2007.mat
rm CN*.txt crust1.* block.desc
