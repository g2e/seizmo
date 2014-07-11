function [varargout]=fisher(varargin)
%FISHER    Converts correlation coefficients to the Z statistic
%
%    Usage:    z=fisher(r)
%              z=fisher(c)
%              [z1,z2,...,zN]=fisher(r1,r2,...,rN)
%
%    Description:
%     Z=FISHER(R) uses the Fisher transform to convert correlation
%     coefficients R to the z statistic Z.  Correlation coefficients should
%     be bounded by -1 and 1.  Any correlation coefficient outside the
%     range is set to NaN.
%
%     Z=FISHER(C) performs the Fisher transform on the complex input C.
%     This is done by taking the transform on the complex modulus
%     multiplied by the complex unit vector of the elements of C.  Note
%     that a complex modulus larger than 1 will return a complex NaN.
%
%     [Z1,Z2,...,ZN]=FISHER(R1,R2,...,RN) performs a multidimensional
%     Fisher transform assuming the N-dimensional data R=<R1,R2,...,RN> are
%     contained by the N-dimensional unit hypersphere.  If the vector
%     length of any element of R is larger than 1 the result for those
%     elements will be NaNs.
%
%    Notes:
%     - Multidimensional Fisher transform inspired by (6) of:
%        Prieto et al 2009, JGR, doi:10.1029/2008JB006067
%
%    Examples:
%     % Convert correlations to z statistic get the mean and
%     % std dev and convert them back to correlation values:
%     z=fisher(r);
%     z_mean=mean(z);
%     z_std=std(z);
%     r_lower=ifisher(z_mean-z_std);
%     r_mean=ifisher(z_mean);
%     r_upper=ifisher(z_mean+z_std);
%
%    See also: IFISHER

%     Version History:
%        Sep.  9, 2009 - minor doc update
%        Mar. 11, 2010 - fixed example
%        Apr.  4, 2012 - minor doc update
%        May  27, 2014 - simplify algorithm for speed/simplicity
%        July 10, 2014 - multidimensional tranform added
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 10, 2014 at 10:35 GMT

% todo:

% SPLIT CODE BY NUMBER OF INPUTS
if(nargin==1)
    % REQUIRE NUMERIC
    if(~isnumeric(varargin{1}))
        error('seizmo:fisher:badInput','FISHER requires numeric input!');
    end
    
    % INPUT IS REAL
    if(isreal(varargin{1}))
        % CONVERT CORRELATIONS OUT OF RANGE TO NaN
        varargin{1}(abs(varargin{1})>1)=nan;
        
        % FISHER TRANSFORM -> CORRELATION SPACE TO Z STATISTIC SPACE
        %Z=0.5*log((1+R)./(1-R)); % OLD WAY
        varargout{1}=atanh(varargin{1});
    else % INPUT IS COMPLEX
        % COMPLEX MODULUS
        cm=abs(varargin{1});
        
        % CONVERT CORRELATIONS OUT OF RANGE TO NaN
        cm(cm>1)=nan;
        
        % FISHER TRANSFORM -> CORRELATION SPACE TO Z STATISTIC SPACE
        varargout{1}=varargin{1}.*atanh(cm)./cm;
    end
else % MULTIDIMENSIONAL
    % REQUIRE NUMERIC
    if(~all(cellfun('isreal',varargin)))
        error('seizmo:fisher:badInput','FISHER requires numeric input!');
    end
    
    % EXPAND AND COMBINE
    [varargin{:}]=expandscalars(varargin{:});
    sumdim=ndims(varargin{1})+1;
    varargin{1}=cat(sumdim,varargin{:});
    varargin(2:end)=[];
    
    % VECTOR LENGTH
    vl=sqrt(sum(varargin{1}.^2,sumdim));
    
    % CONVERT CORRELATIONS OUT OF RANGE TO NaN
    vl(vl>1)=nan;
    
    % FISHER TRANSFORM -> CORRELATION SPACE TO Z STATISTIC SPACE
    varargin{1}=varargin{1}.*submat(atanh(vl)./vl,sumdim,ones(1,nargin));
    
    % SPLIT
    for i=1:nargin; varargout{i}=submat(varargin{1},sumdim,i); end
end

end
