function [varargout]=geofkstructfix(fk)
%GEOFKSTRUCTFIX    Fix old geofk struct for the current geofk struct format
%
%    Usage:    geofk=geofkstructfix(geofk)
%
%    Description:
%     GEOFK=GEOFKSTRUCTFIX(GEOFK) upgrades a struct with a previous layout
%     to the current one.  This does the following currently:
%      .response ==> .beam
%                ==> .method & .center & .npairs
%
%    Notes:
%     - .npairs is probably not be accurate for volumes produced with
%       GEOFKXCVOLUME & GEOFKXCHORZVOLUME.
%
%    Examples:
%
%    See also: CHKGEOFKSTRUCT, FKSTRUCTFIX, CHKFKSTRUCT

%     Version History:
%        July  8, 2010 - initial version
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 13:05 GMT

% check nargin
error(nargchk(1,1,nargin));

% require a struct
if(~isstruct(fk))
    error('seizmo:geofkstructfix:badInput','GEOFK must be a struct!');
end

% require old fields
reqfields={'response' 'nsta' 'stla' 'stlo' 'stel' 'stdp' 'butc' 'eutc' ...
    'npts' 'delta' 'latlon' 'horzslow' 'freq' 'normdb' 'volume'};
fields=fieldnames(fk);
if(any(~ismember(reqfields,fields)))
    error('seizmo:geofkstructfix:badInput',...
        'Old GEOFK struct does not appear to be valid!');
end

% what is the difference?
% .response ==> .beam
%  ==> .method & .center & .npairs
for i=1:numel(fk)
    for j=1:numel(fields)
        switch fields{j}
            case 'response'
                varargout{1}(i).beam=fk(i).response;
            otherwise
                varargout{1}(i).(fields{j})=fk.(fields{j});
        end
    end
    varargout{1}(i).method='coarray';
    [clat,clon]=arraycenter(fk(i).stla,fk(i).stlo);
    varargout{1}(i).center=[clat clon];
    varargout{1}(i).npairs=fk(i).nsta*(fk(i).nsta-1)/2;
end

end
