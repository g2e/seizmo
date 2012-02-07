function [mt,mo]=mt_norm(mt)
%MT_NORM    Returns moment tensors normalized by their scalar moment
%
%    Usage:    [mt,mo]=mt_norm(mt)
%
%    Description:
%     [MT,MO]=MT_NORM(MT) normalizes the moment tensors given by the 3x3xN
%     or Nx6 array MT by its scalar moment.  The output MT gives the
%     normalized moment tensors and MO contains their scalar moments.  This
%     function is useful for comparison and decomposition.
%
%    Notes:
%     - Warns for tensors with Mo=0.
%
%    Examples:
%     % Can you explain this histogram?
%     [mt,mo]=mt_norm(mt_s2g(findcmts));
%     hist(mo,100)
%
%    See also: MT_DIAG, MT_DECOMP, SCALARMOMENT

%     Version History:
%        June  7, 2011 - initial version
%        Feb.  7, 2012 - add warning about Mo=0 tensors
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  7, 2012 at 13:50 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check tensor format
mtsz=size(mt);
if(~isnumeric(mt) || ~isreal(mt))
    error('seizmo:mt_norm:badInput',...
        'MT must be a real-valued numeric array!');
elseif(isequal(mtsz(1:2),[3 3]) && any(numel(mtsz)==[2 3]))
    if(numel(mtsz)>2); n=mtsz(3); else n=1; end
    convert=false;
elseif(mtsz(2)==6 && numel(mtsz)==2)
    mt=mt_v2g(mt); % convert from Nx6 to 3x3xN
    n=mtsz(1);
    convert=true;
else
    error('seizmo:mt_norm:badInput',...
        'MT must be a harvard moment tensor array as 3x3xN or Nx6!');
end

% have to work one at a time
mo=nan(n,1);
for i=1:n
    % get scalar moment
    mo(i,1)=sqrt(trace(mt(:,:,i)*mt(:,:,i))/2);
    
    % normalize
    % - skip mo==0 (single couple?)
    if(mo(i,1)); mt(:,:,i)=mt(:,:,i)/mo(i,1); end
end

% warn about Mo==0
if(any(~mo))
    warning('seizmo:mt_norm:noMoment',...
        ['Mo=0 for some tensors so they are unnormalized:\n' ...
        sprintf('%d ',find(~mo))]);
end

% convert back if input was Nx6
if(convert); mt=mt_g2v(mt); end

end
