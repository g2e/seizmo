function [r]=fkcorr(vol1,vol2)
%FKCORR    Zero lag correlation coefficient between fk spectra
%
%    Usage:    r=fkcorr(vol1,vol2)
%
%    Description:
%     R=FKCORR(VOL1,VOL2) computes the correlation coefficient between the
%     two fk structs VOL1 & VOL2.  VOL1 & VOL2 must have the same slowness
%     and frequency dimensions (the .x, .y, & .freq fields must have the
%     same size).  VOL1 or VOL2 may also be fk maps or fk ARFs.
%
%    Notes:
%
%    Examples:
%     % Compare fk results for 2 different frequency bands:
%     vol1=fkmap(data,50,201,1./[20 15]);
%     vol2=fkmap(data,50,201,1./[15 10]);
%     r=fkcorr(vol1,vol2)
%
%    See also: GEOFKCORR, FKVOLUME, FKMAP, FKARF

%     Version History:
%        Sep. 23, 2010 - initial version
%        Oct. 10, 2010 - add docs, better checking, allow ARFs
%        Oct. 12, 2010 - proper output
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 12, 2010 at 20:00 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check inputs
if(~isfkstruct(vol1) && ~isfkarfstruct(vol1))
    error(chkfkstruct(vol1));
elseif(~isfkstruct(vol1) && ~isfkarfstruct(vol1))
    error(chkfkstruct(vol2));
end

% require scalar
if(~isscalar(vol1) || ~isscalar(vol2))
    error('seizmo:fkcorr:badInput',...
        'VOL1 and VOL2 must be scalar fk structs!');
end

% require volumes are equal in size
if(~isequal(size(vol1.beam),size(vol2.beam)))
    error('seizmo:fkcorr:badInput',...
        'FK volumes must have equal dimensions!');
end

% to correlate we need to convert to linear space
% - y=10*log10(x) => x=10.^(y/10);
r=corrcoef([10.^(vol1.beam(:)/10) 10.^(vol2.beam(:)/10)]);
r=r(2); % only want the correlation coeff for VOL1 vs VOL2

end
