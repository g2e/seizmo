function [data,x,t]=deconvolve_source_timefunction(data,varargin)
%DECONVOLVE_SOURCE_TIMEFUNCTION   Deconvolve source from SEIZMO records
%
%    Usage:    data=deconvolve_source_timefunction(data,hwidth)
%              data=deconvolve_source_timefunction(data,hwidth,type)
%              data=deconvolve_source_timefunction(data,hwidth,type,h2o)
%              data=deconvolve_source_timefunction(...
%                       data,hwidth,type,h2o,frange)
%              data=deconvolve_source_timefunction(...
%                       data,hwidth,type,h2o,frange,zi)
%              data=deconvolve_source_timefunction(...
%                       data,hwidth,type,h2o,frange,zi,zf)
%              [data,x,t]=deconvolve_source_timefunction(...)
%
%    Description: DATA=DECONVOLVE_SOURCE_TIMEFUNCTION(DATA,HWIDTH)
%     deconvolves a gaussian source function with halfwidth HWIDTH from
%     records in SEIZMO struct DATA.  HWIDTH is in seconds and must be a
%     scalar or an array of size equal to the number of records in DATA.
%     The gaussian source function has unit area so that the operation
%     preserves the energy in the records (see the last Usage format to get
%     the actual source function).
%
%     DATA=DECONVOLVE_SOURCE_TIMEFUNCTION(DATA,HWIDTH,TYPE) specifies the
%     type of source function used in the deconvolution.  See function
%     MAKE_SOURCE_TIMEFUNCTION for valid values.  TYPE should be a string
%     like 'gaussian' or 'triangle' or a cellstr array with one string per
%     record in DATA.
%
%     See function DECONVOLVE for details on the following options:
%     DATA=DECONVOLVE_SOURCE_TIMEFUNCTION(DATA,HWIDTH,TYPE,H2O)
%     DATA=DECONVOLVE_SOURCE_TIMEFUNCTION(DATA,HWIDTH,TYPE,H2O,FRANGE)
%     DATA=DECONVOLVE_SOURCE_TIMEFUNCTION(DATA,HWIDTH,TYPE,H2O,FRANGE,ZI)
%     DATA=DECONVOLVE_SOURCE_TIMEFUNCTION(DATA,HWIDTH,TYPE,H2O,FRANGE,ZI,ZF)
%
%     [DATA,X,T]=DECONVOLVE_SOURCE_TIMEFUNCTION(...) also returns the
%     source time functions used in the deconvolution in cell arrays X and
%     T.  See MAKE_SOURCE_TIMEFUNCTION for more info.
%
%    Notes:
%     - gaussian-type functions extend from about -1.5*HWIDTH to 1.5*HWIDTH
%     - triangle-type functions extend from about -HWIDTH to HWIDTH
%
%    Header changes: DEPMIN, DEPMEN, DEPMAX
%
%    Examples:
%     Deconvolve a gaussian with a halfwidth of 10s from some data:
%      data1=deconvolve_source_timefunction(data,10,'gaussian');
%
%    See also: DECONVOLVE, CONVOLVE_SOURCE_TIMEFUNCTION, CONVOLVE,
%              MAKE_SOURCE_TIMEFUNCTION, TRIANGLETF, GAUSSIANTF

%     Version History:
%        Oct. 29, 2009 - initial version
%        Jan. 30, 2010 - fix checking state functions, better messages
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 30, 2010 at 20:15 GMT

% todo:

% check nargin
msg=nargchk(2,7,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check header
    data=checkheader(data);
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

% attempt deconvolution
try
    % number of records
    nrecs=numel(data);
    
    % get header info
    [b,e,delta]=getheader(data,'b','e','delta');
    leven=getlgc(data,'leven');
    iftype=getenumid(data,'iftype');
    
    % cannot do spectral/xyz records
    if(any(~strcmpi(iftype,'itime') & ~strcmpi(iftype,'ixy')))
        error('seizmo:deconvolve_source_timefunction:badIFTYPE',...
            ['Record(s):\n' sprintf('%d ',...
            find(~strcmpi(iftype,'itime') & ~strcmpi(iftype,'ixy'))) ...
            '\nDatatype of record(s) in DATA must be Timeseries or XY!']);
    end
    
    % cannot do unevenly sampled records
    if(any(strcmpi(leven,'false')))
        error('seizmo:deconvolve_source_timefunction:badLEVEN',...
            ['Record(s):\n' sprintf('%d ',find(strcmpi(leven,'false'))) ...
            '\nInvalid operation on unevenly sampled record(s)!']);
    end
    
    % pass to make_source_timefunction
    [x,t]=make_source_timefunction(delta,varargin{1:min(nargin-1,2)});
    
    % get delay (in samples!)
    delay=nan(nrecs,1);
    for i=1:nrecs; delay(i)=t{i}(1); end
    delay=round(delay./delta);
    
    % deconvolve records
    data=deconvolve(data,x,delay,varargin{3:min(nargin-1,6)});
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror)
end

end
