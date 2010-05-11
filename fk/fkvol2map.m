function [vol]=fkvol2map(vol,frng)
%FKVOL2MAP    Converts a fk volume to a fk map
%
%    Usage:    map=fkvol2map(vol,frng)
%
%    Description: MAP=FKVOL2MAP(VOL,FRNG) averages the fk response data in
%     VOL across the frequency range given by FRNG, returning the result in
%     MAP.  MAP may be then plotted using PLOTFKMAP.  VOL is a fk struct
%     produced either by FKVOLUME or FK4D.  FRNG gives the frequency range
%     to average over as [FREQLOW FREQHIGH] in Hz.
%
%    Notes:
%
%    Examples:
%     Show frequency averaged energies for a dataset at 20-50s periods:
%      svol=fkvolume(data,50,201,[1/50 1/20]);
%      smap=fkvol2map(svol,[1/50 1/20]);
%      plotfkmap(smap);
%
%    See also: FKFREQSLIDE, PLOTFKMAP, FKVOLUME, FK4D, FKMAP

%     Version History:
%        May  11, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  11, 2010 at 14:50 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check fk struct
error(chkfkstruct(vol));

% don't allow array/map
if(any(~[vol.volume]))
    error('seizmo:fkvol2map:badInput',...
        'VOL must be a fk volume response!');
end

% check frequency range
if(nargin<2); frng=[]; end
sf=size(frng);
if(~isempty(frng) && (~isreal(frng) || numel(sf)~=2 || sf(1)~=1 ...
        || sf(2)~=2 || any(frng(:)<=0)))
    error('seizmo:fkvol2map:badInput',...
        'FRNG must be a 1x2 array of [FREQLOW FREQHIGH] in Hz!');
end

% average frequencies in frng
for i=1:numel(vol)
    if(isempty(frng))
        fmin=min(vol(i).z);
        fmax=max(vol(i).z);
    else
        fmin=frng(1);
        fmax=frng(2);
    end
    vol(i).volume=false;
    fidx=vol(i).z>=fmin & vol(i).z<=fmax;
    vol(i).response=sum(vol(i).response(:,:,fidx),3)/sum(fidx);
    vol(i).z=vol(i).z(fidx);
end

end
