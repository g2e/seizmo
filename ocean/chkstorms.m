function [report]=chkstorms(storms)
%CHKSTORMS    Validate a storms struct
%
%    Usage:    msg=chkstorms(storms)
%
%    Description:
%     MSG=CHKSTORMS(STORMS) checks that STORMS is a struct as output by
%     READ_HURDAT or READ_GISS_STORMDB and can be plotted by MAPSTORMS.
%     See those functions for details on the struct layout.  MSG is an
%     error structure if a problem is found (otherwise it is empty).
%
%    Notes:
%
%    Examples:
%     % Validate the tropical cyclone dataset:
%     storms=load('tstorms.mat');
%     error(chkstorms(storms));
%
%    See also: ISSTORMS, READ_HURDAT, READ_GISS_STORMDB, MAPSTORMS

%     Version History:
%        Feb. 16, 2013 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 16, 2013 at 13:30 GMT

% todo:

% check nargin
report=[];
error(nargchk(1,1,nargin));

% require scalar struct
if(~isstruct(storms) || ~isscalar(storms))
    report.identifier='seizmo:chkstorm:dataNotStruct';
    report.message='STORMS must be a scalar struct!';
    return;
end

% check model has the required fields
reqfields={'time' 'lat' 'lon'};
if(~all(ismember(reqfields,fieldnames(storms))))
    report.identifier='seizmo:chkstorm:badFields';
    report.message=[...
        sprintf('STORMS must be a struct with fields:\n') ...
        sprintf('%s ',reqfields{:})];
    return;
end

% check each required field
if(~isequal(numel(storms.time),numel(storms.lat),numel(storms.lon)))
    report.identifier='seizmo:chkstorm:fieldsNotSynced';
    report.message='STORMS fields must be equal in size!';
    return;
elseif(size(storms.time,2)~=1 || ndims(storms.time)>2)
    report.identifier='seizmo:chkstorm:fieldsNotColumnVectors';
    report.message='STORMS fields must be column vectors!';
    return;
elseif(~iscell(storms.time) || ~iscell(storms.lat) || ~iscell(storms.lon))
    report.identifier='seizmo:chkstorm:fieldsNotCellArrays';
    report.message='STORMS fields must be cell arrays!';
    return;
elseif(any(~cellfun('isreal',storms.time)) ...
        || any(~cellfun('isreal',storms.lat)) ...
        || any(~cellfun('isreal',storms.lon)))
    report.identifier='seizmo:chkstorm:cellsNotReal';
    report.message='STORMS field cells must have real-valued vectors!';
    return;
elseif(any(cellfun('ndims',storms.time)~=2) ...
        || any(cellfun('ndims',storms.lat)~=2) ...
        || any(cellfun('ndims',storms.lon)~=2) ...
        || any(cellfun('size',storms.time,2)~=1) ...
        || any(cellfun('size',storms.lat,2)~=1) ...
        || any(cellfun('size',storms.lon,2)~=1))
    report.identifier='seizmo:chkstorm:realVectorsNotColumns';
    report.message='STORMS field contents not column vectors!';
    return;
elseif(~isequal(cellfun('size',storms.time,1), ...
        cellfun('size',storms.lat,1),cellfun('size',storms.lon,1)))
    report.identifier='seizmo:chkstorm:realVectorsNotSynced';
    report.message='STORMS field contents not synced!';
    return;
end

end
