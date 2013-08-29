%% SEIZMO - A Matlab(R) & Octave Toolbox for Passive Seismology
%
%% CMB_CLUSTERING - Outlier Removal & Grouping
% Cluster analysis with <matlab:helpwin('cmb_clustering') CMB_CLUSTERING>
% uses the cross correlation info between the recordings of an earthquake
% from <matlab:helpwin('cmb_1st_pass') CMB_1ST_PASS> or
% <matlab:helpwin('cmb_2nd_pass') CMB_2ND_PASS> to help separate the
% waveforms into groups.  This also helps to remove strong signals that are
% unusual and are probably noise.  By separating the waveforms into groups
% with strong similarity, we can be assured that the cross-correlation
% based measurements (relative arrival times) are valid and low in error.
% To help understand how to use cluster analysis effectively, we go through
% the interface below.

%% When should I use CMB_CLUSTERING?
%
% Clustering the waveforms is generally helpful and should be done for all
% earthquakes.  For some earthquakes you will not remove any waveforms or
% split the waveforms into multiple groups, so CMB_CLUSTERING has no effect
% on the analysis of these earthquakes.  That is normal and okay.

%% The First Menu - A Number of Clustering Options
%
% Immediately upon running CMB_CLUSTERING you are presented with a menu and
% with a plot showing both the waveforms and a dendrogram colored according
% to the default dissimilarity cutoff (which is 0.2).  The menu has a
% number of options but you will mostly use 3 of the buttons.  A list of
% the buttons and if you will use them:
% (1) Yes
% (2) Maybe
% (3) Maybe
% (4) No
% (5) No
% (6) No
% (7) Yes
% (8) Yes
