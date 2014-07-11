function [varargout]=ifisher(varargin)
%IFISHER    Converts Z statistics to correlation coefficients
%
%    Usage:    r=ifisher(z)
%              [r1,r2,...,rN]=ifisher(z1,z2,...,zN)
%
%    Description:
%     R=IFISHER(Z) uses the inverse Fisher transform to convert z
%     statistics Z to correlation coefficients R.  Correlation coefficients
%     will be bounded by -1 and 1.  NaNs will be preserved as NaNs.  If Z
%     is complex then the output is also complex (see FISHER for details on
%     the complex Fisher transform.
%
%     [R1,R2,...,RN]=IFISHER(Z1,Z2,...,ZN) performs the multidimensional
%     inverse Fisher transform such that R=<R1,R2,...,RN> are contained by
%     the N-dimensional unit hypersphere.
%
%    Notes:
%     - Multidimensional inverse Fisher transform inspired by (6) of:
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
%    See also: FISHER

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
        error('seizmo:ifisher:badInput','IFISHER requires numeric input!');
    end
    
    % INPUT IS REAL
    if(isreal(varargin{1}))
        % INVERSE FISHER TRANSFORM -> Z STATISTIC TO CORRELATION SPACE
        %R=(exp(2*Z)-1)./(exp(2*Z)+1);
        %R(Z==inf)=1; % otherwise inf returns nan
        varargout{1}=tanh(varargin{1});
    else % INPUT IS COMPLEX
        % COMPLEX MODULUS
        cm=abs(varargin{1});
        
        % INVERSE FISHER TRANSFORM -> Z STATISTIC TO CORRELATION SPACE
        varargout{1}=varargin{1}.*tanh(cm)./cm;
    end
else % MULTIDIMENSIONAL
    % REQUIRE NUMERIC
    if(~all(cellfun('isreal',varargin)))
        error('seizmo:ifisher:badInput','IFISHER requires numeric input!');
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
    varargin{1}=varargin{1}.*submat(tanh(vl)./vl,sumdim,ones(1,nargin));
    
    % SPLIT
    for i=1:nargin; varargout{i}=submat(varargin{1},sumdim,i); end
end
