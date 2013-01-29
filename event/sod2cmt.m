function [cmt,dist]=sod2cmt(sod,varargin)
%SOD2CMT    Returns a GlobalCMT struct from a SOD event csv struct
%
%    Usage:    [cmt,dist]=sod2cmt(sod)
%              [cmt,dist]=sod2cmt(sod,'opt1',val1,...,'optN',valN)
%
%    Description:
%     [CMT,DIST]=SOD2CMT(SOD) returns the most likely GlobalCMT centroid
%     moment tensor (cmt) for each of the specified events.  SOD is a
%     scalar structure formatted following the output from READSODEVENTCSV.
%     Each cmt is found using default FINDCMT settings.
%
%     [CMT,DIST]=SOD2CMT(SOD,'OPT1',VAL1,...,'OPTN',VALN) passes options on
%     to FINDCMT.  See FINDCMT for further details.
%
%    Notes:
%
%    Examples:
%     % Create a SOD event csv file for all available GlobalCMT CMTs:
%     writesodeventcsv('globalcmt.csv',cmt2sod(findcmts));
%
%    See also: CMT2SOD, FINDCMTS, FINDCMT,
%              READSODEVENTCSV, WRITESODEVENTCSV

%     Version History:
%        Feb. 29, 2012 - initial version
%        Jan. 18, 2013 - magnitude type bugfix, preallocation bugfix
%        Jan. 28, 2013 - force column vector output
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 28, 2013 at 17:25 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check input
if(~isstruct(sod) && ~isscalar(sod))
    error('seizmo:sod2cmt:badInput',...
        ['SOD must be a scalar struct formatted ' ...
        'like READSODEVENTCSV output!']);
end
req={'time' 'latitude' 'longitude' 'depth' 'magnitude' 'magnitudeType'};
fields=fieldnames(sod);
if(~all(ismember(req,fields)))
    error('seizmo:sod2cmt:missingFields',...
        ['SOD struct must have the following fields:\n' ...
        sprintf('%s ',req{:})]);
end

% loop over each event
nev=numel(sod.latitude);
dist=nan(nev,1);
for i=1:nev
    sod1=ssidx(sod,i);
    [cmt(i),dist(i)]=findcmt('time',sod1.time,...
        'location',[sod1.latitude sod1.longitude],...
        'magnitude',sod1.magnitude,'magtype',sod1.magnitudeType{:},...
        'depth',sod1.depth,varargin{:}); %#ok<AGROW>
    if(i==1); cmt(2:nev,1)=cmt; end
end

% concatenate struct elements into a scalar struct
cmt=sscat(cmt);

end
