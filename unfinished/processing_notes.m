
% Pdiff slow/decay pre processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. extract data from seed volume
% 2. fix header info
%    - need cmts
% 3. merge
% 4. detrend
% 5. taper
% 6. resample to 1sps
% 7. remove response to velocity
% 8. split into vert & horz datasets
% 9. rotate horz
% 10. add arrivals
% 11. timeshift to diff
%
% Pdiff initial processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. filter between 25s & 60s
% 2. get SNR estimate
% 3. align
% 4. scale
% 5. cluster
%
% Pdiff narrow band processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. loop over clusters
% 2. align
% 3. scale
%
% Pdiff corrections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. mantle (data only)
% 2. crust (data only)
% 3. ellipticity (data only)
% 4. geometrical spreading
% 5. Q-factor (synthetics only)