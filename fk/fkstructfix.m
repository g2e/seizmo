function [varargout]=fkstructfix(fk)
%FKSTRUCTFIX    Migrates an old fk struct to the current fk struct format
%
%    Usage:    fk=fkstructfix(fk)
%
%    Description:
%     FK=FKSTRUCTFIX(FK) upgrades a struct with a previous layout to the
%     current one.  This does the following currently:
%      .response ==> .beam
%      .z        ==> .freq
%      .center   ==> .method & .center & .npairs
%
%    Notes:
%     - .npairs is probably not accurate for volumes produced with
%       FKXCVOLUME & FKXCHORZVOLUME.
%
%    Examples:
%
%    See also: CHKFKSTRUCT, GEOFKSTRUCTFIX, CHKGEOFKSTRUCT

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
    error('seizmo:fkstructfix:badInput','FK must be a struct!');
end

% require old fields
reqfields={'response' 'nsta' 'stla' 'stlo' 'stel' 'stdp' 'butc' 'eutc' ...
    'npts' 'delta' 'x' 'y' 'z' 'polar' 'center' 'normdb' 'volume'};
fields=fieldnames(fk);
if(any(~ismember(reqfields,fields)))
    error('seizmo:fkstructfix:badInput',...
        'Old FK struct does not appear to be valid!');
end

% what is the difference?
% .z ==> .freq
% .response ==> .beam
% .center ==> .method & .center & .npairs
for i=1:numel(fk)
    for j=1:numel(fields)
        switch fields{j}
            case 'response'
                varargout{1}(i).beam=fk(i).response;
            case 'z'
                varargout{1}(i).freq=fk(i).z;
            case 'center'
                % fix method/center/npairs
                if(ischar(fk(i).center))
                    varargout{1}(i).method=lower(fk(i).center);
                    switch fk(i).center
                        case 'coarray'
                            varargout{1}(i).npairs=fk(i).nsta ...
                                *(fk(i).nsta-1)/2;
                        case 'full'
                            varargout{1}(i).npairs=fk(i).nsta*fk(i).nsta;
                        case 'center'
                            varargout{1}(i).npairs=fk(i).nsta;
                    end
                    [clat,clon]=arraycenter(fk(i).stla,fk(i).stlo);
                    varargout{1}(i).center=[clat clon];
                else
                    varargout{1}(i).center=fk(i).center;
                    varargout{1}(i).method='user';
                    varargout{1}(i).npairs=fk(i).nsta;
                end
            otherwise
                varargout{1}(i).(fields{j})=fk.(fields{j});
        end
    end
end

end
