function [data]=clip(data,thresh,value)
%CLIP    Clip SEIZMO data above a threshhold
%
%    Usage:    data=clip(data,thresh)
%              data=clip(data,thresh,value)
%
%    Description:
%     DATA=CLIP(DATA,THRESH) clips values in the records in SEIZMO struct
%     DATA above the value THRESH.  The clipped data are replaced with
%     THRESH.
%
%     DATA=CLIP(DATA,THRESH,VALUE) uses VALUE to replace the clipped data.
%
%    Notes:
%
%    Examples:
%     % Clip anything above 3*RMS:
%     rms=getvaluefun(data,@(x)sqrt(mean(x.^2)));
%     data=clip(data,3*rms);
%
%     % Replace clipped sections with NaNs and then interpolate them
%     % (Note this requires the function INPAINT_NANS from Matlab FEX):
%     data=solofun(clip(data,thresh,nan),@inpaint_nans);
%
%    See also: SOLOFUN

%     Version History:
%        Jan. 12, 2012 - initial version
%        June 11, 2012 - fix rms formula in example, skip double conversion
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 11, 2012 at 11:15 GMT

% check nargin
error(nargchk(2,3,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% number of records
nrecs=numel(data);

% default value to thresh if none
if(nargin<3 || isempty(value)); value=thresh; end

% check thresh & value
if(~isnumeric(thresh) || ~isreal(thresh) || ~any(numel(thresh)==[1 nrecs]))
    error('seizmo:clip:badInput',...
        'THRESH must be a real-valued scalar or vector w/ 1 value/record');
end
thresh=abs(thresh); % force positive
if(isscalar(thresh)); thresh(1:nrecs,1)=thresh; end
if(~isnumeric(value) || ~isreal(value) || ~any(numel(value)==[1 nrecs]))
    error('seizmo:clip:badInput',...
        'VALUE must be a real-valued scalar or vector w/ 1 value/record');
end
value=abs(value); % force positive
if(isscalar(value)); value(1:nrecs,1)=value; end

% verbosity
verbose=seizmoverbose;

% detail message
if(verbose)
    disp('Adding Constant to Record(s)');
    print_time_left(0,nrecs);
end

% loop over records and clip
[depmin,depmen,depmax]=deal(nan(nrecs,1));
for i=1:nrecs
    % skip dataless
    if(isempty(data(i).dep))
        % detail message
        if(verbose); print_time_left(i,nrecs); end
        continue;
    end
    
    % clip
    clipped=abs(data(i).dep)>thresh(i);
    data(i).dep(clipped)=value(i)*sign(data(i).dep(clipped));
    
    % dep*
    depmen(i)=nanmean(data(i).dep(:)); 
    depmin(i)=min(data(i).dep(:)); 
    depmax(i)=max(data(i).dep(:));
    
    % detail message
    if(verbose); print_time_left(i,nrecs); end
end

% update header
data=changeheader(data,'depmin',depmin,'depmen',depmen,'depmax',depmax);

end
