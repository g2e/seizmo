function [s]=fsssub(s,frng,fidx)
%FSSSUB    Extracts a subspectra of a fss struct
%
%    Usage:    s=fsssub(s,frng)
%              s=fsssub(s,frng,fidx)
%
%    Description:
%     S=FSSSUB(S,FRNG) reduces the fss spectra in S to the frequency range
%     given by FRNG and returns S with all pertinent info updated.  S is a
%     struct (see FSS for details).  FRNG specifies the frequency range as
%     [FREQLOW FREQHIGH] in Hz.  By default FRNG is set to [] (full range).
%
%     S=FSSSUB(S,FRNG,FIDX) also reduces the fss spectra in S to the
%     frequencies with the indices corresponding to FIDX.  Note that only
%     those with within the range defined by FRNG are returned.  The
%     default is all frequencies.
%
%    Notes:
%
%    Examples:
%     % Split a spectra into 10s period windows:
%     s=fss(d,50,101,[1/50 1/20]);
%     s1=fsssub(s,[1/50 1/40]);
%     s2=fsssub(s,[1/40 1/30]);
%     s3=fsssub(s,[1/30 1/20]);
%
%    See also: FSS, FSSXC, FSSHORZ, FSSHORZXC, FSSAVG, PLOTFSS,
%              FSSFREQSLIDE, FSSFRAMESLIDE

%     Version History:
%        June 25, 2010 - initial version
%        July  6, 2010 - update for new struct
%        Apr.  4, 2012 - minor doc update
%        June  8, 2012 - adapted from geofksubvol
%        Sep. 12, 2012 - adapted from geofsssub
%        Sep. 22, 2012 - added fidx option
%        Sep. 29, 2012 - fix fidx bug
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 29, 2012 at 19:05 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% check struct
error(chkfss(s));

% don't allow averaged spectrums
for a=1:numel(s)
    if(size(s(a).spectra,3)~=numel(s(a).freq))
        error('seizmo:fsssub:badInput',...
            'S(%d) is already averaged: cannot extract a subspectrum!',a);
    end
end

% check frequency range
if(nargin<2); frng=[]; end
if(nargin<3); fidx=[]; end
sf=size(frng);
if(~isempty(frng) && (~isreal(frng) || numel(sf)~=2 || sf(1)~=1 ...
        || sf(2)~=2 || any(frng(:)<=0)))
    error('seizmo:fsssub:badInput',...
        'FRNG must be a 1x2 array of [FREQLOW FREQHIGH] in Hz!');
elseif(~islogical(fidx) && ~isnumeric(fidx))
    error('seizmo:fsssub:badInput',...
        'FIDX must be indices!');
end

% extract subspectra
for i=1:numel(s)
    % freq range
    if(isempty(frng))
        fmin=min(s(i).freq);
        fmax=max(s(i).freq);
    else
        fmin=frng(1);
        fmax=frng(2);
    end
    fidx0=s(i).freq>=fmin & s(i).freq<=fmax;
    
    % direct frequency indexing
    if(isempty(fidx))
        % do nothing
    elseif(islogical(fidx))
        if(~isequal(size(fidx0),size(fidx)))
            error('seizmo:fsssub:badInput',...
                'FIDX logical indexing does not match frequency size!');
        end
        fidx0=fidx0 & fidx;
    else
        if(any(fidx>numel(fidx0)) || any(fidx<0))
            error('seizmo:fsssub:badInput',...
                'FIDX indices outside indexing range!');
        end
        fidx0=intersect(find(fidx0),fidx);
    end
    
    % extract subspectra
    s(i).spectra=s(i).spectra(:,:,fidx0);
    
    % update info
    s(i).freq=s(i).freq(fidx0);
end

end
