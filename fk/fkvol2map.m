function [vol]=fkvol2map(vol,frng)
%FKVOL2MAP    Converts a fk volume to a fk map
%
%    Usage:    map=fkvol2map(vol,frng)
%
%    Description:
%     MAP=FKVOL2MAP(VOL,FRNG) averages the fk beam data in VOL across the
%     frequency range given by FRNG, returning the result in MAP.  MAP may
%     be then plotted using PLOTFKMAP.  VOL is a fk struct produced either
%     by FKVOLUME or FK4D.  FRNG gives the frequency range to average over
%     as [FREQLOW FREQHIGH] in Hz.
%
%    Notes:
%     - Averages across dB values rather than the underlying values.  This
%       provides slightly better results.  This doesn't matter much because
%       the underlying values are not quite right anyways for 'full' or
%       'coarray': the real values were already absoluted.  (I should just
%       save the true real values and force other functions to deal with
%       this...meh)
%
%    Examples:
%     % Show frequency averaged energies for a dataset at 20-50s periods:
%     svol=fkvolume(data,50,201,[1/50 1/20]);
%     smap=fkvol2map(svol,[1/50 1/20]);
%     plotfkmap(smap);
%
%    See also: FKFREQSLIDE, PLOTFKMAP, FKVOLUME, FK4D, FKMAP

%     Version History:
%        May  11, 2010 - initial version
%        May  23, 2010 - renormalize (set normdb appropriately)
%        July  6, 2010 - update for new struct
%        Nov. 18, 2010 - include code for averaging across true values
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 18, 2010 at 16:05 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check fk struct
error(chkfkstruct(vol));

% don't allow map
if(any(~[vol.volume]))
    error('seizmo:fkvol2map:badInput',...
        'VOL must be a fk beam volume!');
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
        fmin=min(vol(i).freq);
        fmax=max(vol(i).freq);
    else
        fmin=frng(1);
        fmax=frng(2);
    end
    vol(i).volume=false;
    fidx=vol(i).freq>=fmin & vol(i).freq<=fmax;
    vol(i).freq=vol(i).freq(fidx);
    
    % average across true values
    %vol(i).beam=vol(i).beam+vol(i).normdb;
    %vol(i).beam=10.^(vol(i).beam/10);
    %vol(i).beam=sum(vol(i).beam(:,:,fidx),3)/sum(fidx);
    %vol(i).beam=10*log10(vol(i).beam);
    %vol(i).normdb=max(vol(i).beam(:));
    %vol(i).beam=vol(i).beam-vol(i).normdb;
    
    % average across dB
    vol(i).beam=sum(vol(i).beam(:,:,fidx),3)/sum(fidx);
    maxdb=max(vol(i).beam(:));
    vol(i).beam=vol(i).beam-maxdb;
    vol(i).normdb=vol(i).normdb+maxdb;
end

end
