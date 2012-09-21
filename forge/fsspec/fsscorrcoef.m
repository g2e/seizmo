function [r,p,rlo,rhi]=fsscorrcoef(s1,s2,varargin)
%FSSCORRCOEF    Zero lag correlation coefficient between fss spectra
%
%    Usage:    r=fsscorrcoef(s1,s2)
%              [r,p,rlo,rhi]=fsscorrcoef(s1,s2,...)
%
%    Description:
%     R=FSSCORRCOEF(S1,S2) computes the correlation coefficient between
%     the frequency-slowness spectra in S1 & S2.  S1 & S2 must have spectra
%     with the same dimensions.  Also S1 & S2 must be scalar.
%
%     [R,P,RLO,RHI]=FSSCORRCOEF(S1,S2,...) provides access to additional
%     inputs and outputs of the function CORRCOEF.  See that function for
%     more details.
%
%    Notes:
%
%    Examples:
%     % Compare 2 different frequency bands:
%     s1=fssavg(fss(data,50,201,1./[20 15]));
%     s2=fssavg(fss(data,50,201,1./[15 10]));
%     r=fsscorrcoef(s1,s2)
%
%    See also: FSSCORRCOEF, ARF, FSS, FSSXC, FSSHORZ, FSSHORZXC, CORRCOEF

%     Version History:
%        Sep. 23, 2010 - initial version
%        Oct. 10, 2010 - add docs, better checking, allow ARFs
%        Oct. 12, 2010 - proper output
%        Sep. 14, 2012 - adapted from fkcorr, use corrcoef fully
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 14, 2012 at 20:00 GMT

% todo:

% check nargin
error(nargchk(2,inf,nargin));

% check inputs
error(chkfss(s1));
error(chkfss(s2));

% require scalar spectra
if(~isscalar(s1) || ~isscalar(s2))
    error('seizmo:fsscorrcoef:badInput',...
        'S1 and S2 are not scalar FSS structs!');
end

% require spectra to be equal size
if(~isequal(size(s1.spectra),size(s2.spectra)))
    error('seizmo:fsscorrcoef:badInput',...
        'S1 and S2 spectra are not equally sized!');
end

% pass info to corrcoef
[r,p,rlo,rhi]=corrcoef(s1.spectra(:),s2.sprectra(:),varargin{:});
r=r(2); p=p(2); rlo=rlo(2); rhi=rhi(2);

end
