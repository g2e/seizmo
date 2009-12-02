% how?
%
% need to decide on inputs
% - best if we can have the first input be all the data
% - rest of the inputs would be options
%
% need to decide on options
% - 'to' (sac compatible)
%   - only horizontals + orthogonal
%   - allow 'gcp' or any numeric header field
%       - who's header?  ...only 'gcp' until this is decided
%   - allow any real number (reduce to 0-360)
% - 'through'
%   - require orthogonal
%   - this could get confusing...how do we decide which 2 channels to pick?
%       - option should denote this ('pair')
%           - we need a system here...
%
%              1
%              ^   2
%              |  7
%              | /
%              |/
%              +------> 3
%             
%             so the 3 options are:
%             12, 13, 23 (default)
%
%             but what defines 1,2,3?
%               easy case:
%               - if exist a vertical it is 1
%               - if 2 horizontals - 2 is 90deg ccw from 3
%               general case:
%               - 3 channels:
%                   - 1 is most vertical: max(abs(90-abs(cmpinc)))
%                   - failing that -- the most northward
%                   - 2 and 3 follow from finding 1
%               - 2 channels:
%                   - 1 vs 2 vs 3?
%                       - only handle 2 horizontal case
%                           - might be better to keep dead channels around!
% - '3d'
%   - require 3 linearly independent components
%       - throw warning if not and move on to next
%   - just specify the output axis orientations (allow 1+ channel output)
%       - [inc az] -> [0 0 90 0 90 90 ...]
%       - how could we handle 'gcp' here?
%           - eval?

function [data]=rotate(data,varargin)
%ROTATE    Rotates SEIZMO records that are orthogonal channel pairs
%
%    Usage:    data=rotate(data)
%              data=rotate(...,'to',azimuth,...)
%              data=rotate(...,'include',cmplist,...)
%              data=rotate(...,'exclude',cmplist,...)
%              data=rotate(...,'reverse',true|false,...)
%              data=rotate(...,'adjust',method,...)
%              data=rotate(...,'shiftmax',value,...)
%              data=rotate(...,'shiftunits',units,...)
%              data=rotate(...,'interpolate',method,...)
%              data=rotate(...,'useabsolutetiming',logical,...)
%              data=rotate(...,'timing',standard,...)
%              data=rotate(...,'requiredcharfields',fields,...)
%              data=rotate(...,'requiredrealfields',fields,...)
%              data=rotate(...,'verbose',true|false,...)
%
%    Description:
%
%    Notes:
%
%    Header changes: 
%
%    Examples:
%
%    See also: CHECKORIENTATION, ISORTHOGONAL, ISPARALLEL

%     Version History:
%        Nov.  3, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov.  3, 2009 at 19:50 GMT

% todo:

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% attempt header check
try
    % check header of first
    data=checkheader(data);
    
    % 1 or 2 datasets
    nrecs=numel(data);
    nvarargin=numel(varargin);
    %nrecs1=numel(data1); onedata=true;
    %if(nargin>1 && isseizmo(varargin{1},'dep'))
    %    data2=varargin{1};
    %    varargin(1)=[];
    %    nvarargin=nvarargin-1;
    %    nrecs2=numel(data2);
    %    onedata=false;
    %    data2=checkheader(data2);
    %end
    
    % turn off header checking
    oldcheckheaderstate=get_checkheader_state;
    set_checkheader_state(false);
catch
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

% get global
global SEIZMO

% attempt rest
try
    % check number of variable arguments
    if(mod(nvarargin,2))
        error('seizmo:rotate:unpairedOption','Unpaired Option/Value!');
    end
    
    % valid values for string options
    valid.INTERPOLATE={'spline' 'pchip' 'linear' 'nearest'};
    valid.ADJUST={'longer' 'shorter' 'first' 'last' 'one' 'two'};
    valid.SHIFTUNITS={'seconds' 'intervals'};
    valid.TIMING={'utc' 'tai'};
    
    % defaults
    option.TO='gcp'; % rotate to backazimuth
    option.INCLUDE='*'; % include all components
    option.EXCLUDE=''; % exclude no components
    option.REVERSE=false; % do not reverse orientation of 2nd cmp
    option.INTERPOLATE='spline'; % spline/pchip/linear/nearest
    option.ADJUST='shorter'; % longer/shorter/first/last
    option.SHIFTMAX=0.01; % interval: 0-0.5 , seconds: 0+
    option.SHIFTUNITS='intervals'; % seconds/intervals
    option.TIMING='utc'; % utc/tai
    option.USEABSOLUTETIMING=true; % true/false
    option.REQUIREDCHARFIELDS={'knetwk' 'kstnm' 'khole' 'kcmpnm'};
    option.REQUIREDREALFIELDS={'delta' 'cmpinc' 'cmpaz'};
    option.VERBOSE=seizmoverbose; % default to seizmoverbose state
    
    % get options from SEIZMO global
    ME=upper(mfilename);
    try
        fields=fieldnames(SEIZMO.(ME));
        for i=1:numel(fields)
            if(~isempty(SEIZMO.(ME).(fields{i})))
                option.(fields{i})=SEIZMO.(ME).(fields{i});
            end
        end
    catch
    end
    
    % get options from command line
    for i=1:2:nvarargin
        if(~ischar(varargin{i}))
            error('seizmo:rotate:badInput',...
                'Options must be specified as a string!');
        end
        if(~isempty(varargin{i+1}))
            option.(upper(varargin{i}))=varargin{i+1};
        end
    end
    
    % check options
    fields=fieldnames(option);
    for i=1:numel(fields)
        % get value of field
        value=option.(fields{i});
        
        % specific checks
        switch lower(fields{i})
            case 'shiftmax'
                if(~isnumeric(value) || ~isscalar(value))
                    error('seizmo:rotate:badInput',...
                        '%s must be a scalar real number!',fields{i});
                end
            case {'interpolate' 'adjust' 'shiftunits' 'timing'}
                if(~ischar(value) || size(value,1)~=1 ...
                        || ~any(strcmpi(value,valid.(fields{i}))))
                    error('seizmo:rotate:badInput',...
                        ['%s option must be one of the following:\n'...
                        sprintf('%s ',valid.(fields{i}){:})],fields{i});
                end
            case{'include' 'exclude'}
                % fix char arrays
                if(ischar(value))
                    value=cellstr(value);
                    option.(fields{i})=value;
                end
                if(~iscellstr(value))
                    error('seizmo:rotate:badInput',...
                        '%s must be a cellstr of component names!',...
                        fields{i});
                end
            case {'requiredcharfields' 'requiredrealfields'}
                % fix char arrays
                if(ischar(value))
                    value=cellstr(value);
                    option.(fields{i})=value;
                end
                if(~iscellstr(value))
                    error('seizmo:rotate:badInput',...
                        '%s option must be a cellstr of header fields!',...
                        fields{i});
                end
            case {'reverse' 'useabsolutetiming' 'verbose'}
                if(~islogical(value) || ~isscalar(value))
                    error('seizmo:rotate:badInput',...
                        '%s option must be a logical!',fields{i});
                end
            otherwise
                error('seizmo:rotate:badInput',...
                    'Unknown option: %s !',fields{i});
        end
    end
    
    % get header fields
    if(option.USEABSOLUTETIMING)
        [b,e,delta,npts,ncmp,depmin,depmax,depmen,nz]=getheader(data,...
            'b','e','delta','npts','ncmp','depmin','depmax','depmen','nz');
    else
        [b,e,delta,npts,ncmp,depmin,depmax,depmen]=getheader(data,...
            'b','e','delta','npts','ncmp','depmin','depmax','depmen');
        nz=nan(nrecs,6);
    end
    szreal=size(option.REQUIREDREALFIELDS); reqreal=cell(szreal);
    szchar=size(option.REQUIREDCHARFIELDS); reqchar=cell(szchar);
    if(prod(szreal)~=0)
        [reqreal{:}]=getheader(data,option.REQUIREDREALFIELDS{:});
    end
    if(prod(szchar)~=0)
        [reqchar{:}]=getheader(data,option.REQUIREDCHARFIELDS{:});
    end
    iftype=getenumid(data,'iftype');
    leven=strcmpi(getlgc(data,'leven'),'true');

    % require timeseries and general x vs y
    if(any(~strcmpi(iftype,'itime') & ~strcmpi(iftype,'ixy')))
        error('seizmo:rotate:badRecordType',...
            'Records must be Time Series or General X vs Y !');
    end

    % get start and end of records in absolute time
    if(option.USEABSOLUTETIMING)
        if(strcmp(option.TIMING,'utc'))
            ab=gregorian2modserial(utc2tai(...
                [nz(:,1:4) nz(:,5)+nz(:,6)/1000+b]));
            ae=gregorian2modserial(utc2tai(...
                [nz(:,1:4) nz(:,5)+nz(:,6)/1000+e]));
        else
            ab=gregorian2modserial([nz(:,1:4) nz(:,5)+nz(:,6)/1000+b]);
            ae=gregorian2modserial([nz(:,1:4) nz(:,5)+nz(:,6)/1000+e]);
        end
    else
        ab=[zeros(nrecs,1) b];
        ae=[zeros(nrecs,1) e];
    end
    
    % change real to char
    for i=1:prod(szreal)
        reqreal{i}=num2str(reqreal{i},'%16.16e');
    end

    % make groups (require at least leven and ncmp to be the same)
    [f,h,h]=unique(char(strcat(strcat('',reqchar{:}),'_',...
        strcat('',reqreal{:}),'_',num2str(leven),'_',num2str(ncmp))),...
        'rows');
    
    
    
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    set_checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    set_checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror)
end

end



