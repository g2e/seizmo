function [events]=readsodcsv(filename)
%READSODCSV    Read in a SOD .csv formatted file as an event structure

% time, latitude, longitude, depth, depthUnits, magnitude, magnitudeType, catalog, contributor, name, flinnEngdahlRegion, flinnEngdahlRegionType
%2002-03-03T12:08:03.200Z, 36.535, 70.645, 147.89999389648438, kiloMETER, 6.9, MB, ISCCD, ISC, 11968147, 718, 1
%2002-03-03T12:08:09.700Z, 35.515, 69.499, 225.60000610351562, kiloMETER, 6.9, MS, ISCCD, ISC, 11968151, 718, 1
%2002-03-03T07:16:15.690Z, -45.837, -76.119, 10.0, kiloMETER, 6.0, MW, WHDF, NEIC, 8653027, 143, 1

% check nargin
error(nargchk(1,1,nargin));

% check filename is char
if(~ischar(filename)); error('FILENAME must be a string!'); end

% open file
fid=fopen(filename);

% check if file exists
if(fid<0); error('%s not openable!',filename); end

% get the header
h=textscan(fid,'%s',1,'delimiter','\n','whitespace',''); h=h{:}{:};
h=textscan(h,'%s','delimiter',','); fields=h{:}{:};

% 

end