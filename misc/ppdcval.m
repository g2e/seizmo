function [y]=ppdcval(pp,x,dcpos)
%PPDCVAL    Evaluate piecewise polynomial with discontinuity support
%
%    Usage:    y=ppdcval(pp,x)
%              y=ppdcval(pp,x,dcpos)
%
%    Description:
%     Y=PPDCVAL(PP,X) or Y=PPDCVAL(X,PP) works exactly like PPVAL.  In
%     terms of handling discontinuities in the piecewise polynomial, that
%     means that sites occurring on a break point (potential discontinuity)
%     in the piecewise polynomial are evaluated on the polynomial to the
%     positive side of the break point (discontinuity).
%
%     Y=PPDCVAL(PP,X,DCPOS) allows changing to which side a site on a
%     discontinuity is evaluated.  The default DCPOS is TRUE and works
%     exactly like PPVAL or PPDCVAL without the third input.  Setting DCPOS
%     to FALSE will return values evaluated using the polynomial on the
%     negative side of the break point (discontinuity).
%
%    Notes:
%
%    Examples:
%     % Simple step function piecewise polynomial that is -1 below 0 and
%     % 1 above.  PPVAL assumes 0 occurs on the positive side while PPDCVAL
%     % allows you to choose on which side the break point evaluates:
%     pp=mkpp([-inf 0 inf],[-1 1]);
%     y=ppval(pp,0)
%     y=ppdcval(pp,0,true)
%     y=ppdcval(pp,0,false)
%
%    See also: INTERPDC1, INTERP1, SPLINE, PCHIP, PPVAL, MKPP, UNMKPP

%     Version History:
%        June  1, 2010 - initial version
%        Aug.  8, 2010 - minor doc updates
%        Jan. 11, 2011 - code refactor & doc rewrite, fixes bug in PPVAL
%                        for n-dimension poly & scalar input
%        Jan. 24, 2011 - fix variable name bug
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 24, 2011 at 13:00 GMT

% todo:

% check number of inputs
error(nargchk(2,3,nargin));

% for compatibility with Matlab, allow first two inputs to be switched
if(isstruct(x)); [pp,x]=deal(x,pp); end

% default discontinuity flag
if(nargin<3 || isempty(dcpos)); dcpos=true; end

% check inputs
if(~isscalar(dcpos) || (~isreal(dcpos) && ~islogical(dcpos)))
    error('seizmo:ppdcval:badInput',...
        'DCPOS must be TRUE or FALSE!');
end

% x => row vector
nx=numel(x);
sx=size(x);
xv=x(:).';

% hide 1st dim of x if was a row vector
% - this suppresses the first dimension (1) for row vector input
%   when evaluating polynomials that return multiple values per site
%   thus making it equivalent to if it had been input as a column vector
if(sx(1)==1 && numel(sx)==2); sx(1)=[]; end

% extract info about polynomial
% break points, coefficients, number of pieces, order, dimension
[breaks,coef,l,k,d]=unmkpp(pp);

% figure out which piece each site is in
% - action depends on discontinuity flag choice
% - zero evaluation points handling
if(dcpos) % positive side
    % empty case check
    if(nx)
        % histc returns the positive side so no adjustment necessary
        [p,p]=histc(xv,[-inf breaks(2:l) inf]);
    else % empty (set to 1 so no indexing error)
        p=ones(1,nx);
    end
else % negative side
    % empty case check
    if(nx)
        % adjust sites on break points down one
        [p,p]=histc(xv,[-inf breaks(2:l) inf]);
        onbp=ismember(xv,breaks(2:l));
        p(onbp)=p(onbp)-1;
    else % empty (set to 1 so no indexing error)
        p=ones(1,nx);
    end
end

% go ahead and pre-evaluate NaN/Inf values as NaN
% - -inf is in correct piece
%   +inf needs correction (down one)
%    nan is in none
p(xv==inf)=l;
nans=find(p==0);
p(nans)=1;

% shift site to be relative to minimum site of the piece
xv=xv-breaks(p);

% tile sites and pieces for polynomials that produce non-scalar output
nd=prod(d);
if(nd>1) % vector output
    % replicate downward
    xv=xv(ones(nd,1),:);
    
    % correct piece indexing for site replication
    p=nd*p;
    poff=(nd:-1:1).'-1;
    p=p(ones(nd,1),:)-poff(:,ones(1,nx));
else % scalar output
    % check for row vector case
    if(numel(sx)==1) % row vector
        d=1; % makes output match input
    else
        d=[];
    end
end

% make column vectors
xv=xv(:);
p=p(:);

% evaluate (((a0)*v+a1)*v+a2)*...
y=coef(p,1);
for a=2:k; y=xv.*y+coef(p,a); end

% if the piecewise polynomial is just a single constant
% then c=ppdcval(nan).  If not, then nan=ppdcval(nan).
if(k==1 && l>1); y=reshape(y,[nd nx]); y(:,nans)=nan; end

% shape output with polynomial dimensions followed by input dimensions
y=reshape(y,[d sx]);

% correct orientation if interp was used
% & both polynomial/input are non-scalar
if(nd>1 && nx>1 && isfield(pp,'orient') && strcmpi(pp.orient,'first'))
    % check for vector
    ny=ndims(y);
    if(isvector(x))
        py=[ny 1:ny-1];
    else
        ndx=ndims(x);
        py=[(ny-ndx+1):ny 1:(ny-ndx)];
    end
    y=permute(y,py);
end

end
