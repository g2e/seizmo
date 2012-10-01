function [report]=chkfss(s)
%CHKFSS    Validate if a struct is as defined by FSS functions
%
%    Usage:    msg=chkfss(s)
%
%    Description:
%     MSG=CHKFSS(S) checks that S is a struct as output from FSS functions.
%     See those functions for details on the struct layout.  MSG is an
%     error structure if a problem is found (otherwise it is empty).
%
%    Notes:
%
%    Examples:
%     % Check that output from FSSXC is compatible:
%     s=fssxc(xcdata,50,101,[.01 .0125]);
%     error(chkfss(s));
%
%    See also: CHKGEOFSS, FSS, FSSXC, FSSHORZ, FSSHORZXC, ARF, ARFHORZ,
%              ISFSS

%     Version History:
%        June 22, 2010 - initial version
%        June 24, 2010 - bugfix
%        June 25, 2010 - more lenient to allow a variety of volume types
%        July  6, 2010 - major update for new struct
%        Apr.  4, 2012 - minor doc update
%        June  4, 2012 - adapted from chkgeofkstruct
%        June 10, 2012 - handle full method
%        June 13, 2012 - allow capon method
%        Sep. 12, 2012 - adapt from chkgeofss & chkfkstruct, drop capon
%        Sep. 27, 2012 - add capon, xc, caponxc, tdxc
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 27, 2012 at 18:50 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check map is proper struct
report=[];
fields={'spectra' 'nsta' 'st' 'butc' 'eutc' 'npts' 'delta' 'polar' ...
        'x' 'y' 'freq' 'npairs' 'method' 'center' 'whiten' 'weights'};
if(~isstruct(s) || ~all(ismember(fields,fieldnames(s))))
    report.identifier='seizmo:chkfss:dataNotStruct';
    report.message=['FSS data must be a struct with fields:' ...
        sprintf('\n') sprintf('%s ',fields{:})];
    return;
end

% valid method strings
valid.METHOD={'center' 'coarray' 'user' 'full' ...
    'capon' 'xc' 'caponxc' 'tdxc'};

% loop over each spectra
for i=1:numel(s)
    % check consistency
    ss=size(s(i).spectra);
    if(numel(ss)==2); ss(3)=1; end
    sx=numel(s(i).x);
    sy=numel(s(i).y);
    sf=numel(s(i).freq);
    if(~isreal(s(i).nsta) || ~isscalar(s(i).nsta) ...
            || s(i).nsta~=fix(s(i).nsta) || s(i).nsta<2 ...
            || size(s(i).st,1)~=s(i).nsta || size(s(i).st,2)<2 ...
            || ~isreal(s(i).st))
        report.identifier='seizmo:chkfss:fkStnInfoCorrupt';
        report.message='FSS station info appears corrupt!';
        return;
    elseif(~isequal(size(s(i).butc),[1 5]) || ~isreal(s(i).butc) ...
            || ~isequal(size(s(i).eutc),[1 5]) || ~isreal(s(i).eutc))
        report.identifier='seizmo:chkfss:fkTimingInfoCorrupt';
        report.message='FSS timing info appears corrupt!';
        return;
    elseif(~isreal(s(i).npts) || ~isscalar(s(i).npts) ...
            || s(i).npts~=fix(s(i).npts) || s(i).npts<0 ...
            || ~isreal(s(i).delta) || ~isscalar(s(i).delta) ...
            || s(i).delta<=0)
        report.identifier='seizmo:chkfss:fkWaveformInfoCorrupt';
        report.message='FSS npts/delta appears corrupt!';
        return;
    elseif(~islogical(s(i).polar))
        report.identifier='seizmo:chkfss:polarInvalid';
        report.message='FSS polar field must be TRUE/FALSE!';
        return;
    elseif(~isreal(s(i).npairs) ...
            || s(i).npairs~=fix(s(i).npairs) || s(i).npairs<0)
        report.identifier='seizmo:chkfss:npairsInvalid';
        report.message='FSS npairs field must be a positive integer!';
        return;
    elseif(~any(strcmpi(s(i).method,valid.METHOD)))
        report.identifier='seizmo:chkfss:methodInvalid';
        report.message='FSS method field contains an invalid method!';
        return;
    elseif(~isnumeric(s(i).center) || (~isreal(s(i).center) ...
            || ~numel(s(i).center)==2))
        report.identifier='seizmo:chkfss:centerInvalid';
        report.message='FSS center field must be [LAT LON]!';
        return;
    elseif(~isreal(s(i).x) || ~isreal(s(i).y) ...
            || ~isreal(s(i).freq) || any(s(i).freq<0))
        report.identifier='seizmo:chkfss:xyzCorrupt';
        report.message='FSS x/y/freq info is corrupt!';
        return;
    elseif(~isreal(s(i).spectra) ...
            || (ss(3)>1 && ~isequal(ss,[sy sx sf])) ...
            || (ss(3)==1 && ~isequal(ss,[sy sx 1])))
        report.identifier='seizmo:chkfss:spectraCorrupt';
        report.message='FSS spectra has wrong size or invalid data!';
        return;
    end
end

end
