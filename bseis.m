function [data]=bseis(varargin)
%BSEIS   Arranges xy data into a SAClab data structure
%
%    Description: BSEIS(IND1,DEP1,IND2,DEP2,...) takes arrays of
%     independent and dependent components and arranges them into a SAClab
%     data structure (one record per IND/DEP pair) to be compatible with
%     SAClab functions.  Independent data must be a vector as multiple
%     independent components are not currently supported.  If there are
%     multiple dependent components for an independent datum, they should
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
%        Oct.  2, 2008 - SACLAB global options cleaned up, fixed uneven
%                        datafile detection
%        Oct. 15, 2008 - hasdata field support, possible bugfix for struct
%                        setup
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 15, 2008 at 07:15 GMT

% todo:

% check number of inputs
if (mod(nargin,2)) 
    error('SAClab:bseis:badNargs','Unpaired IND/DEP data!')
end

% defaults
option.PREFHDRVER=6;
option.ENDIAN=[];

% get options from SACLAB global
global SACLAB
if(isfield(SACLAB,'BSEIS'))
    if(isfield(SACLAB.BSEIS,'ENDIAN'))
        if(strcmpi(SACLAB.BSEIS.ENDIAN,{'ieee-le' 'ieee-be'}))
            option.ENDIAN=lower(SACLAB.BSEIS.ENDIAN);
        else
            warning('SAClab:bseis:badOption','ENDIAN incomprehensible!');
        end
    end
    if(isfield(SACLAB.BSEIS,'PREFHDRVER') ...
        && isscalar(SACLAB.BSEIS.PREFHDRVER) ...
        && isnumeric(SACLAB.BSEIS.PREFHDRVER) ...
        && any(vvseis==SACLAB.BSEIS.PREFHDRVER))
        option.PREFHDRVER=SACLAB.BSEIS.PREFHDRVER;
    else
        warning('SAClab:bseis:badOption','PREFHDRVER incomprehensible!');
    end
end

% get the version definition
h=seisdef(option.PREFHDRVER);

% get the native byteorder if needed
if(isempty(option.ENDIAN))
    option.ENDIAN=nativeendian;
end

% undefine numeric header
undef=nan(h.size,1,h.store);
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

% create structure
data(1:nargin/2,1)=struct('name',[],'endian',[],...
    'version',[],'hasdata',[],'head',[],'dep',[]);

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
    
    % fill in other fields
    data(j).endian=option.ENDIAN;
    data(j).version=option.PREFHDRVER;
    data(j).hasdata=true;
    data(j).head=undef;
    
    % handle dataless separately
    if(npts==0)
        data(j)=ch(data(j),'npts',0,'iftype','General X vs Y file',...
            'lovrok','true','nvhdr',option.PREFHDRVER,'knetwk','SAClab');
        continue;
    end
    
    % fill in knowns/presets
    delta=(varargin{i}(end)-varargin{i}(1))/(npts-1);
    data(j)=ch(data(j),...
        'delta',delta,'b',varargin{i}(1),'e',varargin{i}(end),...
        'npts',npts,'depmin',min(data(j).dep(:)),...
        'depmax',max(data(j).dep(:)),'depmen',mean(data(j).dep(:)),...
        'iftype','General X vs Y file','leven','true',...
        'lovrok','true','nvhdr',option.PREFHDRVER,'knetwk','SAClab');
    
    % handle 1pt data
    if(npts==1)
        data(j)=ch(data(j),'leven',nan);
        continue;
    end
    
    % add proper info to unevenly spaced data
    if(any(abs(delta-diff(varargin{i}))>10*eps))
        data(j).ind=varargin{i}(:);
        data(j)=ch(data(j),'leven','false',...
            'odelta',data(j).ind(min([2 end]))-data(j).ind(1));
    end
end

end
