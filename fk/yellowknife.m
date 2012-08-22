function [loc,name]=yellowknife()
%YELLOWKNIFE    Returns locations and names for the Yellowknife array
%
%    Usage:    [location,name]=yellowknife
%
%    Description:
%     [LOCATION,NAME]=YELLOWKNIFE returns the [LAT LON ELEV] and station
%     names of the Yellowknife array.  Latitude and longitude are in
%     degrees, elevation is in meters.  
%
%    Notes:
%
%    Examples:
%     % Plot the array response function for the Yellowknife array:
%     loc=yellowknife;
%     fkarf(loc(:,1:2),20,201,1);
%
%    See also: ARRAYCENTER, ARRAYAPERTURE, FKARF

%     Version History:
%        May  25, 2011 - initial version
%        Aug.  6, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug.  6, 2012 at 23:55 GMT

% locations (lat, lon, depth)
loc=[62.61 -114.61 194.2
     62.40 -114.61 145.1
     62.42 -114.61 150.3
     62.45 -114.61 158.2
     62.47 -114.61 163.8
     62.52 -114.61 173.6
     62.54 -114.61 176.6
     62.56 -114.61 171
     62.58 -114.60 213.1
     62.49 -114.94 141.1
     62.49 -114.90 145
     62.49 -114.85 146.2
     62.49 -114.80 148.9
     62.49 -114.75 154.2
     62.49 -114.70 161
     62.49 -114.65 167.5
     62.49 -114.61 166.7
     62.49 -114.56 171.7
     62.48 -114.48 170.7
     62.43 -114.60 143.8
     62.56 -114.61 170.3
     62.49 -114.74 160.6];
name={'YKB0' 'YKB1' 'YKB2' 'YKB3' 'YKB4' 'YKB6' 'YKB7' 'YKB8' 'YKB9' ...
    'YKR1' 'YKR2' 'YKR3' 'YKR4' 'YKR5' 'YKR6' 'YKR7' 'YKR8' 'YKR9' ...
    'YKW1' 'YKW2' 'YKW3' 'YKW4'}';

end
