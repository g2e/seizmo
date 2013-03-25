function [varargout]=mt_norm(varargin)
%MT_NORM    Returns moment tensors normalized by their scalar moment
%
%    Usage:    [mt,mo]=mt_norm(mt)
%              [m1,m2,m3,m4,m5,m6,mo]=mt_norm(m1,m2,m3,m4,m5,m6)
%
%    Description:
%     [MT,MO]=MT_NORM(MT) normalizes the moment tensors given by the 3x3xN
%     or Nx6 array MT by its scalar moment.  MT may also be a scalar struct
%     as output by FINDCMT or FINDCMTS.  The output MT gives the normalized
%     moment tensors and MO contains their scalar moments.  This function
%     is useful for comparison and decomposition.  The output is the same
%     format as the input.
%
%     [M1,M2,M3,M4,M5,M6,MO]=MT_NORM(M1,M2,M3,M4,M5,M6) allows for moment
%     tensor component input/output.  
%
%    Notes:
%     - Warns for tensors with MO=0.
%     - Struct input only modifies the moment tensor component fields not
%       the scalarmoment or exponent fields.
%
%    Examples:
%     % Can you explain this histogram?
%     [mt,mo]=mt_norm(findcmts);
%     hist(log10(mo),100)
%
%    See also: MT_DIAG, MT_DECOMP, SCALARMOMENT, RADPAT, FINDCMT, FINDCMTS

%     Version History:
%        June  7, 2011 - initial version
%        Feb.  7, 2012 - add warning about MO=0 tensors
%        Mar. 23, 2013 - minor doc update
%        Mar. 25, 2013 - update for mt_check/mt_change
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 25, 2013 at 13:50 GMT

% todo:

% check nargin
error(nargchk(1,6,nargin));

% check tensor format
error(mt_check(varargin{:}));

% handle struct input separately (b/c we have to alter in place)
if(isstruct(varargin{1}))
    mt=varargin{1};
    mo=mt.scalarmoment;
    if(any(mo==0))
        warning('seizmo:mt_norm:noMoment',...
            ['MO=0 for some tensors so they are unnormalized:\n' ...
            sprintf('%d ',find(~mo))]);
    end
    bad=mo==0;
    mo(bad)=1;
    mt.mrr=mt.mrr./mo;
    mt.mtt=mt.mtt./mo;
    mt.mpp=mt.mpp./mo;
    mt.mrt=mt.mrt./mo;
    mt.mrp=mt.mrp./mo;
    mt.mtp=mt.mtp./mo;
    mo(bad)=0;
    varargout={mt mo};
else
    [mt,from]=mt_change('g',varargin{:});
    n=size(mt,3);
    
    % have to work one at a time
    mo=nan(n,1);
    for i=1:n
        % get scalar moment
        mo(i,1)=sqrt(trace(mt(:,:,i)*mt(:,:,i))/2);
        
        % normalize
        % - skip mo==0 (single couple / linear vector dipole)
        if(mo(i,1)); mt(:,:,i)=mt(:,:,i)/mo(i,1); end
    end
    
    % warn about mo==0
    if(any(~mo))
        warning('seizmo:mt_norm:noMoment',...
            ['MO=0 for some tensors so they are unnormalized:\n' ...
            sprintf('%d ',find(~mo))]);
    end
    
    % convert back
    [varargin{:}]=mt_change(from,mt);
    varargout=[varargin {mo}];
end

end
