function [s]=fsssub(s,frng)
%FSSSUB    Extracts a subspectra of a fss struct
%
%    Usage:    s=fsssub(s,frng)
%
%    Description:
%     S=FSSSUB(S,FRNG) reduces the fss spectra in S to the frequency range
%     given by FRNG and returns S with all pertinent info updated.  S is a
%     struct (see FSS for details).  FRNG specifies the frequency range as
%     [FREQLOW FREQHIGH] in Hz.
%
%     By default FRNG is set to [] (full range).
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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 12, 2012 at 19:05 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check struct
error(chkfss(s));

% don't allow averaged spectrums
for a=1:numel(s)
    if(size(s(a).spectra,3)~=numel(s(a).freq))
        error('seizmo:fsssub:badInput',...
            'S(%d) is already averaged: cannot extract a subspectrum!',...
            a);
    end
end

% check frequency range
if(nargin<2); frng=[]; end
sf=size(frng);
if(~isempty(frng) && (~isreal(frng) || numel(sf)~=2 || sf(1)~=1 ...
        || sf(2)~=2 || any(frng(:)<=0)))
    error('seizmo:fsssub:badInput',...
        'FRNG must be a 1x2 array of [FREQLOW FREQHIGH] in Hz!');
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
    fidx=s(i).freq>=fmin & s(i).freq<=fmax;
    
    % extract subspectra
    s(i).spectra=s(i).spectra(:,:,fidx);
    
    % update info
    s(i).freq=s(i).freq(fidx);
end

end
