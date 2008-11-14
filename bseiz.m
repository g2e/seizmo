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
%    Tested on: Matlab r2007b
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
%        Oct. 27, 2008 - update for struct changes, better SACLAB global
%                        handling, made consistent with new requirements,
%                        single CH call
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 27, 2008 at 06:55 GMT

% todo:

% turn off checking
oldseischkstate=get_seischk_state;
set_seischk_state(true);

% check number of inputs
if (mod(nargin,2)) 
    error('SAClab:bseis:badNargs','Unpaired IND/DEP data!')
end

% defaults
option.FILETYPE='SAClab Binary';
option.VERSION=6;
option.ENDIAN=[];

% get options from SACLAB global
global SACLAB
if(isfield(SACLAB,'BSEIS'))
    if(isfield(SACLAB.BSEIS,'ENDIAN'))
        if(strcmpi(SACLAB.BSEIS.ENDIAN,{'ieee-le' 'ieee-be'}))
            option.ENDIAN=lower(SACLAB.BSEIS.ENDIAN);
        else
            warning('SAClab:bseis:badOption',...
                'ENDIAN option from SACLAB global unknown => skipping!');
        end
    end
    if(isfield(SACLAB.BSEIS,'FILETYPE'))
        if(isfield(SACLAB.BSEIS,'VERSION')...
            && isscalar(SACLAB.BSEIS.VERSION)...
            && isnumeric(SACLAB.BSEIS.PREFHDRVER)...
            && any(SACLAB.BSEIS.VERSION==vvseis(SACLAB.BSEIS.FILETYPE)))
            option.FILETYPE=SACLAB.BSEIS.FILETYPE;
            option.VERSION=SACLAB.BSEIS.VERSION;
        else
            warning('SAClab:bseis:badOption',...
                'FILETYPE option from SACLAB global unknown => skipping!');
        end
    end
    if(isfield(SACLAB.BSEIS,'VERSION'))
        if(isscalar(SACLAB.BSEIS.VERSION)...
            && isnumeric(SACLAB.BSEIS.VERSION)...
            && any(SACLAB.BSEIS.VERSION==vvseis(option.FILETYPE)))
            option.VERSION=SACLAB.BSEIS.VERSION;
        else
            warning('SAClab:bseis:badOption',...
                'VERSION option from SACLAB global unknown => skipping!');
        end
    end
end

% get the version definition
h=seisdef(option.FILETYPE,option.VERSION);

% get the native byteorder if needed
if(isempty(option.ENDIAN))
    option.ENDIAN=nativeendian;
end

% undefine numeric header
undef=nan(h.size,1,h.store);
for i=1:length(h.ntype)
    for j=1:length(h.(h.ntype{i}))
        undef(h.(h.ntype{i})(j).minpos:h.(h.ntype{i})(j).maxpos)=...
            h.undef.ntype;
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
nrecs=nargin/2;
data(1:nrecs,1)=struct('location',[],'name',[],'filetype',[],...
    'version',[],'endian',[],'hasdata',[],'head',[],'dep',[]);

% loop for each pair
leven=true(nrecs,1); delta=ones(nrecs,1);
[b,e,npts,ncmp,depmen,depmin,depmax]=swap(nan(nrecs,1));
for i=1:2:nargin
    % output index
    j=ceil(i/2);
    
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
    
    % get npts
    npts(j)=length(varargin{i});
    
    % cross check
    if(isvector(varargin{i+1}))
        % vectors can be row or column vectors
        if(npts(j)~=length(varargin{i+1}))
            error('SAClab:bseis:badInput',...
                ['Dependent data does not match independent data length'...
                 ': pair %d !'],j);
        end
        
        % fill in dependent variable (as column vector)
        data(j).dep=varargin{i+1}(:);
        ncmp(j)=1;
    else
        % arrays must have components oriented down columns
        [ndpts,ncmp(j)]=size(varargin{i+1});
        if(npts(j)~=ndpts)
            error('SAClab:bseis:badInput',...
                ['Dependent data does not match independent data length'...
                 ': pair %d !'],j);
        end
        
        % fill in dependent variable (assure 2D output)
        data(j).dep=varargin{i+1}(:,:);
    end
    
    % edit name
    data(j).name=['SAClab.' num2str(j) '.sac'];
    
    % fill in other fields
    data(j).location='.';
    data(j).endian=option.ENDIAN;
    data(j).filetype=option.FILETYPE;
    data(j).version=option.VERSION;
    data(j).hasdata=true;
    data(j).head=undef;
    
    % handle 0pt
    if(npts(j)==0); continue; end
    
    % get b,e,dep*
    b(j)=varargin{i}(1);
    e(j)=varargin{i}(end);
    depmen(j)=mean(data(j).dep(:));
    depmin(j)=min(data(j).dep(:));
    depmax(j)=max(data(j).dep(:));
    
    % handle 1pt
    if(npts(j)==1); continue; end
    
    % get delta and handle uneven
    delta(j)=diff(varargin{i}([1 end]))/(npts(j)-1);
    if(any(abs(delta(j)-diff(varargin{i}))>10*eps))
        data(j).ind=varargin{i}(:); leven(j)=false;
    end
end

% write header changes
data=ch(data,'b',b,'e',e,'delta',delta,'npts',npts,'ncmp',ncmp,...
    'depmen',depmen,'depmin',depmin,'depmax',depmax,...
    'iftype','General X vs Y file','lovrok','true',...
    'nvhdr',option.VERSION,'knetwk','SAClab','leven',leven);


% toggle checking back
set_seischk_state(oldseischkstate);

end
