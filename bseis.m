function [data]=bseis(varargin)
%BSEIS   Arranges xy data into a SAClab data structure
%
%    Description: BSEIS(IND1,DEP1,IND2,DEP2,...) takes arrays of
%     independent and dependent components and arranges them into a SAClab
%     data structure (one record per IND/DEP pair) to be compatible with
%     SAClab functions.  Independent data must be a vector as multiple
%     independent components are not currently supported.  If there are
%     multiple dependent components for a inddependent datum, they should
%     be arranged such that DEP contains each component in separate
%     columns.
%
%    Notes:
%     - outputs records as equivalent to SAC version 6
%     - the byte-order is set to match the current architecture
%     - the filetype is set as 'General X vs Y file'
%     - automatically figures out if data is evenly sampled
%
%    System requirements: Matlab 7
%
%    Input/Output requirements: numeric arrays
%
%    Header changes: 
%     CREATES HEADER INFO: 
%      DELTA, B, E, NPTS, DEPMEN, DEPMIN, DEPMAX, IFTYPE, LEVEN, LOVROK,
%      NVHDR, KNETWK, and for unevenly spaced data ODELTA
%
%    Usage:    data=bseis(IND1,DEP1,IND2,DEP2...)
%
%    Examples:
%     To create a square root function in Matlab and then convert the array
%     information into a SAClab compatible structure and ultimately write
%     to a formatted binary file (SAC version 6 equivalent):
%      times=linspace(0,30,1000);
%      amps=sqrt(times);
%      data=bseis(times,amps);
%      data.name='myfile';
%      wseis(data);
%
%    See also:  wseis, rseis

%     Version History:
%        Oct. 29, 2007 - initial version
%        Nov.  7, 2007 - doc update
%        Jan. 28, 2008 - complete rewrite, SACHP support
%        Feb. 28, 2008 - renamed from BSAC to BSEIS
%        Feb. 29, 2008 - minor doc update
%        Mar.  2, 2008 - use GENUMDESC
%        Mar.  3, 2008 - fixed unevenly sampled behavior
%        Mar.  4, 2008 - use platform native byte-order
%        May  12, 2008 - DEP* fix
%        June 12, 2008 - output XY by default and doc update
%        June 28, 2008 - fixed default header settings, added dataless
%                        support, records now have names by default,
%                        .dep and .ind rather than .x and .t
%        June 30, 2008 - history fix
%        Sep. 24, 2008 - multi-component support, fixed some behavior bugs,
%                        global options support (alt. header version)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 24, 2008 at 19:40 GMT

% todo:
% - make sure dataless & 1point have LEVEN==undef

% check number of inputs
if (mod(nargin,2)) 
    error('SAClab:bseis:badNargs','Unpaired IND/DEP data!')
end

% preferred SAClab header version
pref=6;

% get preferred header layout
global SACLAB
if(isfield(SACLAB,'PREFERREDHEADERVERSION') ...
        && ~isempty(SACLAB.PREFERREDHEADERVERSION) ...
        && isscalar(SACLAB.PREFERREDHEADERVERSION) ...
        && isnumeric(SACLAB.PREFERREDHEADERVERSION) ...
        && any(vvseis==SACLAB.PREFERREDHEADERVERSION))
    pref=SACLAB.PREFERREDHEADERVERSION;
end
h=seisdef(pref);

% undefine numeric header
undef=zeros(h.size,1,h.store);
for i=1:length(h.ntype)
    for j=1:length(h.(h.ntype{i}))
        undef(h.(h.ntype{i})(j).minpos:h.(h.ntype{i})(j).maxpos)=h.undef.ntype;
    end
end

% undefine char header
for i=1:length(h.stype)
    for j=1:length(h.(h.stype{i}))
        sfields=fieldnames(h.(h.stype{i})(j).pos);
        for k=1:length(sfields)
            m=h.(h.stype{i})(j).pos.(sfields{k});
            n=m(2)-m(1)+1; o=length(h.undef.stype);
            undef(m(1):m(2))=[h.undef.stype repmat(32,1,n-o)];
        end
    end
end

% get this platform's native byte-order
[platform,maxint,endian]=computer;
clear platform maxint
if(strcmpi(endian,'L')); endian='ieee-le';
else endian='ieee-be'; end

% create structure
data(1:nargin/2,1)=struct('version',pref,'endian',endian,...
    'name','','head',undef,'dep',[]);

% loop for each pair
for i=1:2:nargin
    % output index
    j=round((i+1)/2);
    
    % check type
    if(~isnumeric(varargin{i}))
        error('SAClab:bseis:badInput',...
            'Independent data must be numeric: pair %d !',j);
    elseif(~isnumeric(varargin{i+1}))
        error('SAClab:bseis:badInput',...
            'Dependent data must be numeric: pair %d !',j);
    end
    
    % check size
    if(~isvector(varargin{i}))
        error('SAClab:bseis:badInput',...
            'Independent data must be a vector: pair %d !',j);
    end
    npts=length(varargin{i});
    if(isvector(varargin{i+1}))
        % vectors can be row or column vectors
        if(npts~=length(varargin{i+1}))
            error('SAClab:bseis:badInput',...
                ['Dependent data does not match independent data length'...
                 ': pair %d !'],j);
        end
        
        % fill in dependent variable
        data(j).dep=varargin{i+1}(:);
    else
        % arrays must have components oriented down columns
        if(npts~=size(varargin{i+1},1))
            error('SAClab:bseis:badInput',...
                ['Dependent data does not match independent data length'...
                 ': pair %d !'],j);
        end
        
        % fill in dependent variable
        data(j).dep=varargin{i+1};
    end
    
    % edit name
    data(j).name=['SAClab.' num2str(j) '.sac'];
    
    % handle dataless separately
    if(npts==0)
        data(j)=ch(data(j),'npts',0,'iftype','General X vs Y file',...
            'lovrok','true','nvhdr',pref,'knetwk','SAClab');
        continue;
    end
    
    % fill in knowns/presets
    delta=(varargin{i}(end)-varargin{i}(1))/(npts-1);
    data(j)=ch(data(j),...
        'delta',delta,'b',varargin{i}(1),'e',varargin{i}(end),...
        'npts',npts,'depmin',min(data(j).dep(:)),...
        'depmax',max(data(j).dep(:)),'depmen',mean(data(j).dep(:)),...
        'iftype','General X vs Y file','leven','true',...
        'lovrok','true','nvhdr',pref,'knetwk','SAClab');
    
    % handle single point data
    if(npts==1)
        data(j)=ch(data(j),'leven',nan);
        continue;
    end
    
    % add proper info to unevenly spaced data
    if(abs(delta-diff(varargin{i}))>10*eps)
        data(j).ind=varargin{i}(:);
        data(j)=ch(data(j),'leven','false',...
            'odelta',data(j).ind(min([2 end]))-data(j).ind(1));
    end
end

end
