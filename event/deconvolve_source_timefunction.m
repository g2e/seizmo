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
%    Description:
%     DATA=DECONVOLVE_SOURCE_TIMEFUNCTION(DATA,HWIDTH) deconvolves a
%     gaussian source function with halfwidth HWIDTH from records in SEIZMO
%     struct DATA.  HWIDTH is in seconds and must be a scalar or an array
%     of size equal to the number of records in DATA.  The gaussian source
%     function has unit area so that the operation preserves the energy in
%     the records (see the last Usage format to get the actual source
%     function).
%
%     DATA=DECONVOLVE_SOURCE_TIMEFUNCTION(DATA,HWIDTH,TYPE) specifies the
%     type of source function used in the deconvolution.  See function
%     MAKE_SOURCE_TIMEFUNCTION for valid values.  TYPE should be a string
%     like 'gaussian' or 'triangle' or a cellstr array with one string per
%     record in DATA.
%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     SEE FUNCTION DECONVOLVE FOR DETAILS ON THE REMAINING INPUTS!
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%     % Deconvolve a gaussian with a halfwidth of 10s from some data:
%     data1=deconvolve_source_timefunction(data,10,'gaussian');
%
%    See also: DECONVOLVE, CONVOLVE_SOURCE_TIMEFUNCTION, CONVOLVE,
%              MAKE_SOURCE_TIMEFUNCTION, TRIANGLETF, GAUSSIANTF

%     Version History:
%        Oct. 29, 2009 - initial version
%        Jan. 30, 2010 - fix checking state functions, better messages
%        Feb. 11, 2011 - mass nargchk fix, mass seizmocheck fix, use
%                        checkheader more effectively
%        Mar. 24, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 24, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(2,7,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check header
    data=checkheader(data,...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR');
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% attempt deconvolution
try
    % number of records
    nrecs=numel(data);
    
    % get header info
    delta=getheader(data,'delta');
    
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
    error(lasterror);
end

end
