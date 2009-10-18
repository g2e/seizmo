function [data]=convolve_source_timefunction(data,varargin)
%CONVOLVE_SOURCE_TIMEFUNCTION   Convolve source function on SEIZMO records
%
%    Usage:    data=convolve_source_timefunction(data,hwidth)
%              data=convolve_source_timefunction(data,hwidth,type)
%
%    Description: DATA=CONVOLVE_SOURCE_TIMEFUNCTION(DATA,HWIDTH) convolves
%     a gaussian source function onto records in SEIZMO struct DATA.
%     HWIDTH defines the half width of the source function and must be a
%     scalar or an array of size equal to the number of records in DATA.
%     Note that the gaussian is centered on each point, so the convolution
%     is an acausal operation.  The gaussian has unit area so that the
%     energy is preserved (use MAKE_SOURCE_TIMEFUNCTION to get the actual
%     source function).  The returned records include extra points before
%     and after the time limits of the original records.  These points
%     are included because energy has been given to those points through
%     the convolution operation.
%
%     DATA=CONVOLVE_SOURCE_TIMEFUNCTION(DATA,HWIDTH,TYPE) specifies the
%     type of source function to be used in the convolution.  See function
%     MAKE_SOURCE_TIMEFUNCTION for valid values.  TYPE should be a string
%     like 'gaussian' or 'triangle' or a cellstr array with one string per
%     record in DATA.
%
%    Notes:
%     - gaussian-type functions extend from -1.5*HWIDTH to 1.5*HWIDTH
%     - triangle-type functions extend from -HWIDTH to HWIDTH
%
%    Header changes: DEPMIN, DEPMEN, DEPMAX, NPTS, B, E
%
%    Examples:
%     Convolve a 10 second triangle source function onto some synthetic
%     data read into a SEIZMO dataset:
%      data=convolve_source_timefunction(data,10,'triangle');
%
%    See also: MAKE_SOURCE_TIMEFUNCTION, TRIANGLETF, GAUSSIANTF

%     Version History:
%        Oct. 17, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 18, 2009 at 00:25 GMT

% todo:

% check nargin
msg=nargchk(2,3,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% attempt header check
try
    % check header
    data=checkheader(data);
    
    % turn off header checking
    oldcheckheaderstate=get_checkheader_state;
    set_checkheader_state(false);
catch
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    rethrow(lasterror)
end

% attempt rest
try
    % number of records
    nrecs=numel(data);
    
    % get header info
    delta=getheader(data,'delta');
    leven=getlgc(data,'leven');
    iftype=getenumid(data,'iftype');
    
    % cannot do spectral/xyz records
    if(any(~strcmpi(iftype,'itime') & ~strcmpi(iftype,'ixy')))
        error('seizmo:convolve_source_timefunction:badIFTYPE',...
            'Datatype of records in DATA must be Timeseries or XY!');
    end
    
    % cannot do unevenly sampled records
    if(any(strcmpi(leven,'false')))
        error('seizmo:convolve_source_timefunction:badLEVEN',...
            'Invalid operation on unevenly sampled records!');
    end
    
    % pass to make_source_timefunction
    [x,t]=make_source_timefunction(delta,varargin{:});
    
    % convolve with records
    [data,zf]=convolve(data,x);
    
    % attach final conditions
    data=attach(data,'ending',zf);
    
    % get/apply time shift
    t1=nan(nrecs,1);
    for i=1:nrecs; t1(i)=t{i}(1); end
    data=timeshift(data,t1);
    
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    set_checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    set_checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    rethrow(lasterror)
end

end
