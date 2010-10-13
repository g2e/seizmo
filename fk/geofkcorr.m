function [r]=geofkcorr(vol1,vol2)
%GEOFKCORR    Zero lag correlation coefficient between geofk spectra
%
%    Usage:    r=geofkcorr(vol1,vol2)
%
%    Description:
%     R=GEOFKCORR(VOL1,VOL2) computes the correlation coeffient between the
%     two geofk structs VOL1 & VOL2.  VOL1 & VOL2 must have the same input
%     dimensions (number of slowness points, latlon points, and frequency
%     points).  VOL1 and VOL2 may also be geofk maps or geofk ARFs.
%
%    Notes:
%
%    Examples:
%     % Compare geofk results for 2 different frequency bands:
%     vol1=geofkxcvolume(data,50,201,1./[20 15]);
%     vol2=geofkarf(data,50,201,1./[15 10]);
%     r=fkcorr(vol1,vol2)
%
%    See also: FKCORR

%     Version History:
%        Sep. 23, 2010 - initial version
%        Oct. 10, 2010 - add docs, better checking
%        Oct. 12, 2010 - proper output
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 12, 2010 at 20:00 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check inputs
if(~isgeofkstruct(vol1) && ~isgeofkarfstruct(vol1))
    error(chkgeofkstruct(vol1));
elseif(~isgeofkstruct(vol1) && ~isgeofkarfstruct(vol1))
    error(chkgeofkstruct(vol2));
end

% require volumes are equal in size
if(~isequal(size(vol1.beam),size(vol2.beam)))
    error('seizmo:geofkcorr:badInput',...
        'GEOFK volumes must have equal dimensions!');
end

% require volumes are equal in size
if(~isequal(size(vol1.beam),size(vol2.beam)))
    error('seizmo:geofkcorr:badInput',...
        'GEOFK volumes must have equal dimensions!');
end

% to correlate we need to convert to linear space
% - y=10*log10(x) => x=10.^(y/10);
r=corrcoef([10.^(vol1.beam(:)/10) 10.^(vol2.beam(:)/10)]);
r=r(2); % only want the correlation coeff for VOL1 vs VOL2

end
