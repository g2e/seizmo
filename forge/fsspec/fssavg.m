function [s]=fssavg(s,frng)
%FSSAVG    Average fss output over frequency
%
%    Usage:    s=fssavg(s,frng)
%
%    Description:
%     S=FSSAVG(S,FRNG) averages the frequency-slowness spectra in S over
%     the frequency range given by FRNG.  S must be a struct as output by
%     FSS.  FRNG is the frequency range to average as [FREQLOW FREQHIGH] in
%     Hz.
%
%    Notes:
%
%    Examples:
%     % Average over frequencies with about 25s period:
%     s=fss(d,50,101,1./[25.1 24.9]);
%     plotfss(fssavg(s));
%
%    See also: FSS, FSSXC, FSSHORZ, FSSHORZXC, FSSSUB, PLOTFSS

%     Version History:
%        June 24, 2010 - initial version
%        July  6, 2010 - update for new struct
%        Apr.  4, 2012 - minor doc update
%        June  7, 2012 - altered from geofkvol2map
%        Sep. 12, 2012 - altered from geofssavg
%        Sep. 14, 2012 - use nanmean instead of sum when averaging
%        Sep. 29, 2012 - error slightly earlier if no freq in frng
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 29, 2012 at 19:05 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check struct
error(chkfss(s));

% check frequency range
if(nargin<2); frng=[]; end
sf=size(frng);
if(~isempty(frng) && (~isreal(frng) || numel(sf)~=2 || sf(1)~=1 ...
        || sf(2)~=2 || any(frng(:)<0)))
    error('seizmo:fssavg:badInput',...
        'FRNG must be a 1x2 array of [FREQLOW FREQHIGH] in Hz!');
end

% average spectrum over frequencies/slownesses
for i=1:numel(s)
    % handle empty frng
    if(isempty(frng))
        fmin=min(s(i).freq);
        fmax=max(s(i).freq);
    else
        fmin=frng(1);
        fmax=frng(2);
    end
    
    % average (if possible)
    if(numel(s(i).freq)==size(s(i).spectra,3))
        fidx=s(i).freq>=fmin & s(i).freq<=fmax;
        if(~any(fidx))
            error('seizmo:fssavg:badFRNG',...
                'FRNG outside S(%d) frequency range!',i);
        end
        s(i).spectra=nanmean(s(i).spectra(:,:,fidx),3);
        s(i).freq=s(i).freq(fidx);
    end
end

end
