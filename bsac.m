function [data]=bsac(varargin)
%BSAC   Arrange timeseries into SAClab structure format
%
%    Description: Takes vectors of x-values and y-values and formats them 
%     to be compatible with SAClab routines such as rsac, wsac, lh, gh, ch.
%
%    Usage:    saclab_struct=bsac(x_vec1,y_vec1,x_vec2,y_vec2...)
%
%    Examples:
%     To create a square root function in Matlab and then convert the array
%     information into a SAClab compatible structure and ultimately write
%     to a SAC formatted binary file:
%
%     xarray=linspace(0,30,1000);
%     yarray=sqrt(xarray);
%     data=bsac(xarray,yarray);
%     data.name='myfile';
%     wsac(data);
%
%    by Michael Thorne (5/2004)   mthorne@asu.edu
%       Garrett Euler  (2/2008)   ggeuler@wustl.edu
%
%    See also:  wsac, lh, ch, gh, rsac, sachp 

% check number of inputs
if (mod(nargin,2))
    error('Unpaired x/y vectors!')
end

% preferred SAC/SAClab version
pref=6;

% get preferred header layout
h=sachi(pref);

% undefined numeric header
undef=zeros(h.size,1,h.store);
for i=1:length(h.ntype)
    for j=1:length(h.(h.ntype{i}))
        undef(h.(h.ntype{i})(j).minpos:h.(h.ntype{i})(j).maxpos)=h.(h.ntype{i})(j).undef;
    end
end

% undefined char header
for i=1:length(h.stype)
    for j=1:length(h.(h.stype{i}))
        sfields=fieldnames(h.(h.stype{i})(j).pos);
        for k=1:length(sfields)
            m=h.(h.stype{i})(j).pos.(sfields{k});
            n=m(2)-m(1)+1; o=length(h.(h.stype{i})(j).undef);
            undef(m(1):m(2))=[h.(h.stype{i})(j).undef repmat(32,1,n-o)];
        end
    end
end

% create structure
data(1:nargin/2,1)=struct('version',pref,'endian','ieee-le','name',[],...
    'head',undef,'x',[]);

% loop for each pair
for i=1:2:nargin
    % output index
    j=round((i+1)/2);
    
    % check vector lengths
    if(~isvector(varargin{i}) || length(varargin{i})<2)
        error('xarray is not a vector')
    end
    if(~isvector(varargin{i+1}) || length(varargin{i+1})<2)
        error('yarray is not a vector')
    end
    npts=length(varargin{i});
    if(npts~=length(varargin{i+1}))
        error('x and y series are not the same length')
    end
    
    % ensure column vectors
    varargin{i}=varargin{i}(:);
    varargin{i+1}=varargin{i+1}(:);
    
    % fill in dependent variable
    data(j).x(:,1)=varargin{i+1};
    
    % fill in knowns/presets
    delta=(varargin{i}(end)-varargin{i}(1))/(npts-1);
    data(j)=ch(data(j),...
        'delta',delta,'b',varargin{i}(1),'e',varargin{i}(end),...
        'npts',npts,'depmin',-norm(min(data(j).x)),...
        'depmax',norm(max(data(j).x)),'depmen',norm(mean(data(j).x)),...
        'iftype','itime','leven',h.true,'lcalda',h.true,...
        'lovrok',h.true,'lpspol',h.false,'nvhdr',pref,...
        'knetwk','SAClab');
    
    % if timeseries is unevenly spaced add proper info
    if(abs(delta-(varargin{i}(2)-varargin{i}(1)))>eps)
        data(j).t(:,1)=varargin{i};
        data(j)=ch(data(j),'leven',h.false,...
            'odelta',varargin{i}(2)-varargin{i}(1));
    end
end

end
