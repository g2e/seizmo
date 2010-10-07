function [r]=fkcorr(vol1,vol2)
%FKCORR    Zero lag correlation coefficient between fk spectra
%
%    Usage:    r=fkcorr(vol1,vol2)
%
%    Description:
%     R=FKCORR(VOL1,VOL2)
%
%    Notes:
%
%    Examples:
%
%    See also: GEOFKCORR

%     Version History:
%        Sep. 23, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 23, 2010 at 20:00 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check inputs
error(chkfkstruct(vol1));
error(chkfkstruct(vol2));

% require volumes are equal in size
if(~isequal(size(vol1.beam),size(vol2.beam)))
    error('seizmo:fkcorr:badInput',...
        'FK volumes must have equal dimensions!');
end

% to correlate we need to convert to linear space
% - y=10*log10(x) => x=10.^(y/10);
r=corrcoef([10.^(vol1.beam(:)/10) 10.^(vol2.beam(:)/10)]);

end
