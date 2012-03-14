function [n,d]=rrat(x,tol)
%RRAT    Relative rational approximation
%
%    Usage:    [n,d]=rrat(x,tol)
%              s=rrat(x,tol)
%
%    Description:
%     [N,D]=RRAT(X,TOL) returns the numerator and denominator matrices
%     that express the elements in X as a fraction of two small integers.
%     The integer fraction will be within TOL*X of X.  TOL is optional and
%     is by default 1e-6*norm(X(:),1).  This differs from RAT in that TOL
%     actually is relative to X (note that the RAT documentation indicates
%     that it is but in actuality it is not) rather than an absolute
%     tolerance value.
%
%     RRAT(X,TOL) or S=RRAT(X,TOL) returns the continued fraction expansion
%     as a string.
%
%    Notes:
%
%    Examples:
%     % this shows the benefit of relative vs absolute tolerance
%     a=1/300+rand/1e8;
%     rat(a)
%     rrat(a)
%     rat(a,1e-2)
%     rrat(a,1e-2)
%
%    See also: RAT, FORMAT, RATS

%     Version History:
%        Aug. 16, 2010 - gplv3 version derived from GNU Octave's RAT
%        Mar. 13, 2012 - drop license comment, reformat code, many bugfixes
%
%     Written by Paul Kienzle (RAT in GNU Octave)
%                Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 13, 2012 at 11:00 GMT

% todo:

% check io
error(nargchk(1,2,nargin));
error(nargchk(0,2,nargout));

% get size, make vector and find infinite elements
sx=size(x);
nx=numel(x);
x=x(:);

% Replace Inf with nan while calculating ratios.
isinfx=isinf(x);
infx=x(isinfx);
x(isinfx)=0;

% detect zeros
iszerox=x==0;

% default norm
if(nargin<2); tol=1e-6*norm(x,1); end

% First step in the approximation is the integer portion

% First element in the continued fraction.
n=round(x);
d=ones(nx,1);
frac=x-n;
lastn=ones(nx,1);
lastd=zeros(nx,1);
steps=zeros(nx,0);

% nan the zeros
x(iszerox)=nan;

% Grab new factors until all continued fractions converge.
while(1) % goes until break is reached
    % Determine which fractions have not yet converged.
    % edit by Garrett Euler - makes code match described behavior
    % (tol is in relative terms, not absolute terms)
    idx=find(abs((x-n./d)./x)>=tol);
    if(isempty(idx))
        if(isempty(steps)); steps=NaN.*ones(nx,1); end
        break;
    end
    
    % Grab the next step in the continued fraction.
    flip=1./frac(idx);
    % Next element in the continued fraction.
    step=round(flip);
    
    if(nargout<2)
        tsteps=NaN.*ones(nx,1);
        tsteps(idx)=step;
        steps=[steps tsteps];
    end
    
    frac(idx)=flip-step;
    
    % Update the numerator/denominator.
    nextn=n;
    nextd=d;
    n(idx)=n(idx).*step+lastn(idx);
    d(idx)=d(idx).*step+lastd(idx);
    lastn=nextn;
    lastd=nextd;
end

if(nargout==2)
    % Move the minus sign to the top.
    n=n.*sign(d);
    d=abs(d);
    
    % Return the same shape as you receive.
    n=reshape(n,sx);
    d=reshape(d,sx);
    
    % Use 1/0 for Inf.
    n(isinfx)=sign(infx);
    d(isinfx)=0;
    
    % Reshape the output.
    n=reshape(n,sx);
    d=reshape(d,sx);
else
    n='';
    nsteps=size(steps,2);
    x(iszerox)=0;
    x(isinfx)=infx;
    for i=1:nx
        s=int2str(x(i));
        j=1;
        
        while(1) % goes until break is reached
            k=j; j=j+1;
            if(isnan(steps(i,k))); break; end
            if(j>nsteps || isnan(steps(i,k)))
                if(steps(i,k)<0)
                    s=cat(2,s,[' + 1/(' int2str(steps(i,k)) ')']);
                else
                    s=cat(2,s,[' + 1/' int2str(steps(i,k))]);
                end
                break;
            else
                s=cat(2,s,[' + 1/(' int2str(steps(i,k)) ')']);
            end
        end
        s=[s repmat(')',1,j-2)];
        if(~isempty(n))
            n_nc=size(n,2);
            s_nc=size(s,2);
            if(n_nc>s_nc)
                s(:,s_nc+1:n_nc)=' ';
            elseif(s_nc>n_nc)
                n(:,n_nc+1:s_nc)=' ';
            end
        end
        n=cat(1,n,s);
    end
end

end
