function [vol]=fksubvol(vol,frng)
%FKSUBVOL    Extracts a frequency-based subvolume of a fk volume
%
%    Usage:    vol=fksubvol(vol,frng)
%
%    Description:
%     VOL=FKSUBVOL(VOL,FRNG) extracts the fk beam data in VOL for the
%     frequency range given by FRNG, returning the result.  VOL is a fk
%     struct produced either by FKVOLUME or FK4D.  FRNG gives the frequency
%     range to extract as [FREQLOW FREQHIGH] in Hz.
%
%    Notes:
%
%    Examples:
%     % Split a volume into 10s period windows:
%     svol=fkvolume(data,50,201,[1/50 1/20]);
%     svol1=fksubvol(svol,[1/50 1/40]);
%     svol2=fksubvol(svol,[1/40 1/30]);
%     svol3=fksubvol(svol,[1/30 1/20]);
%     svol4=fksubvol(svol,[1/20 1/10]);
%
%    See also: FKFREQSLIDE, PLOTFKMAP, FKVOLUME, FK4D, FKMAP, FKVOL2MAP

%     Version History:
%        June 11, 2010 - initial version
%        July  6, 2010 - update for new struct
%        July 14, 2010 - bugfix related to new struct
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 12:55 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check fk struct
error(chkfkstruct(vol));

% don't allow map
if(any(~[vol.volume]))
    error('seizmo:fksubvol:badInput',...
        'VOL must be a fk beam volume!');
end

% check frequency range
if(nargin<2); frng=[]; end
sf=size(frng);
if(~isempty(frng) && (~isreal(frng) || numel(sf)~=2 || sf(1)~=1 ...
        || sf(2)~=2 || any(frng(:)<=0)))
    error('seizmo:fksubvol:badInput',...
        'FRNG must be a 1x2 array of [FREQLOW FREQHIGH] in Hz!');
end

% extract frequencies in frng
for i=1:numel(vol)
    if(isempty(frng))
        fmin=min(vol(i).freq);
        fmax=max(vol(i).freq);
    else
        fmin=frng(1);
        fmax=frng(2);
    end
    fidx=vol(i).freq>=fmin & vol(i).freq<=fmax;
    vol(i).beam=vol(i).beam(:,:,fidx);
    maxdb=max(vol(i).beam(:));
    vol(i).beam=vol(i).beam-maxdb;
    vol(i).normdb=vol(i).normdb+maxdb;
    vol(i).freq=vol(i).freq(fidx);
end

end
