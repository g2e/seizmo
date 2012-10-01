function [st,nm]=lasa()
%LASA    Returns the long-period LASA station locations and names
%
%    Usage:    [st,kstnm]=lasa
%
%    Description:
%     [ST,KSTNM]=LASA returns the long period station locations and names
%     for the Large Aperture Seismic Array (LASA).  Locations are returned
%     in ST as a 21x2 matrix of [LAT LON] in degrees.  Names are the
%     station names in a 21x1 cellstr vector.  Note that the locations are
%     not exact (digitized from a figure in an article from 1965).
%
%    Notes:
%
%    Examples:
%     % The 0.1Hz array response of LASA:
%     plotarf(arf(lasa,50,201,0,0,1/10));
%     % and now compare that to a weighted cross-spectra version:
%     plotarf(arf(lasa,50,201,0,0,1/10,[],'coarray',pairweights(lasa)));
%
%    See also: YELLOWKNIFE, ARF

%     Version History:
%        Sep. 26, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 26, 2012 at 15:05 GMT

% locations are just digitized from a figure
st=[46.68251 -106.22628
    46.75065 -106.09788
    46.63272 -106.16863
    46.65631 -106.32063
    46.76113 -106.24987
    46.83975 -106.13456
    46.66679 -106.01664
    46.57245 -106.25511
    46.72969 -106.38352
    46.83451 -105.89085
    46.50169 -106.01664
    46.55410 -106.48835
    46.93409 -106.35732
    47.14898 -106.05332
    46.50955 -105.36935
    46.14791 -106.35469
    46.75589 -106.91288
    47.37435 -105.20687
    45.90157 -105.49252
    45.98019 -107.09895
    47.36911 -106.84999];
nm={'A0' 'B1' 'B2' 'B3' 'B4' ...
    'C1' 'C2' 'C3' 'C4' ...
    'D1' 'D2' 'D3' 'D4' ...
    'E1' 'E2' 'E3' 'E4' ...
    'F1' 'F2' 'F3' 'F4'}.';
end
