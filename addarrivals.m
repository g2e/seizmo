function [data]=addarrivals(data,varargin)
%ADDARRIVALS    Adds the indicated phases to the SEIZMO data header
%
%    Usage:    data=addarrivals(data)
%              data=addarrivals(...,'m|mod|model',modelname,...)
%              data=addarrivals(...,'p|ph|phase|phases',phaselist,...)
%              data=addarrivals(...,'f|field|fields',indices,...)
%
%    Description: ADDARRIVALS(DATA) calls TAUPTIME to insert phase arrival
%     times into the headers of records in SEIZMO struct DATA.  The model
%     used to find phase arrival times is IASP91. The default phase list is
%     'ttall', which grabs nearly all phases of interest (more than can be
%     put into the SAC data header).  The first 10 arrivals from this list
%     are inserted into the t, kt, and user fields (t gets arrival times,
%     kt gets the phase names, and user gets the ray parameters).  Event
%     depth and distance are determined from the header fields evdp and
%     gcarc or dist or stla+stlo+evla+evlo.
%
%     ADDARRIVALS(...,'M|MOD|MODEL',MODELNAME,...) sets the model for the
%     arrival time calculation to MODELNAME.  MODELNAME must be a cellstr
%     or char array of one model per record or a single model for all
%     records.  Accepts lots of different 1D models (see TauP for how to
%     add more models).  By default MODELNAME is 'iasp91'.
%
%     ADDARRIVALS(...,'P|PHASE|PHASES',PHASELIST,...) sets the phases to be
%     added.  PHASELIST must be a cellstr or char array with one comma-
%     separated phase list per record or a single comma-separated list for
%     all records.  Accepts lots of different phases (see TauP for a list
%     and conventions).  By default PHASELIST is 'ttbasic', which returns
%     a list similar to Brian Kennett's ttimes program.
%
%     ADDARRIVALS(...,'F|FIELD|FIELDS',INDICES,...) indicates the indices
%     of the t,kt,user header fields to put arrivals in.  The default is
%     0:9.  The indices are for all records and can not be individually
%     set.
%
%    Notes:
%     - TauP was written by:
%       H. Philip Crotwell, Thomas J. Owens, Jeroen Ritsema
%     - matTaup was written by:
%       Qin Li
%
%    Header changes: T, KT, USER
%
%    Examples:
%     Fill t, kt, user fields 6:9 with several S phases:
%      data=addarrivals(data,'phases','tts+','fields',6:9);
%
%    See also: GETARRIVAL, TAUPTIME, CHANGEHEADER

%     Version History:
%        June 29, 2009 - initial version
%        Aug. 25, 2009 - description update (forgot fields option)
%        Sep.  2, 2009 - now uses tauptime
%        Dec.  9, 2009 - works with newer tauptime
%        Jan. 26, 2010 - seizmoverbose support, properly SEIZMO handling
%        Feb.  3, 2010 - versioninfo caching
%        Mar.  8, 2010 - versioninfo caching dropped
%        Apr. 20, 2010 - doc update, flexible option fields, drop global
%        Aug. 21, 2010 - update undef checks
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 21, 2010 at 23:45 GMT

% todo:
% - really need to catch taup messages

% check nargin
if(mod(nargin-1,2))
    error('seizmo:addarrivals:badNumInputs',...
        'Bad number of arguments!');
end

% check data structure & headers
data=checkheader(data);

% toggle off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt adding arrivals
try
    % verbosity
    verbose=seizmoverbose;
    
    % check all options preceeded by char field
    if(~iscellstr(varargin(1:2:end)))
        error('seizmo:addarrivals:badInput',...
            'Options must be specified as a strings!');
    end
    
    % throw defaults in front
    varargin=[{'m' 'iasp91' 'p' 'ttbasic' 'f' 0:9} varargin];

    % check options
    nrecs=numel(data);
    for i=1:2:numel(varargin)
        % specific checks
        switch lower(varargin{i})
            case {'m' 'mod' 'model'}
                if(iscellstr(varargin{i+1}))
                    varargin{i+1}=char(varargin{i+1});
                end
                if(~ischar(varargin{i+1}) ...
                        || ~any(size(varargin{i+1},1)==[1 nrecs]))
                    error('seizmo:addarrivals:badInput',...
                        ['MODEL must be a cellstr/char array with one\n'...
                        'model per record or a single model for all!']);
                end
                if(size(varargin{i+1},1)==1)
                    varargin{i+1}=varargin{i+1}(ones(nrecs,1),:);
                end
                model=cellstr(varargin{i+1});
            case {'p' 'ph' 'phase' 'phases'}
                if(iscellstr(varargin{i+1}))
                    varargin{i+1}=char(varargin{i+1});
                end
                if(isempty(varargin{i+1}) || ...
                        ~ischar(varargin{i+1}) ...
                        || ~any(size(varargin{i+1},1)==[1 nrecs]))
                    error('seizmo:addarrivals:badInput',...
                        ['PHASES must be a cellstr/char array w/ one\n'...
                        'comma-separated phase list per record or a\n'...
                        'single comma-separated phase list for all!']);
                end
                if(size(varargin{i+1},1)==1)
                    varargin{i+1}=varargin{i+1}(ones(nrecs,1),:);
                end
                phases=cellstr(varargin{i+1});
            case {'f' 'field' 'fields'}
                if(~isempty(varargin{i+1}) && (numel(varargin{i+1})>10 ...
                        || any(fix(varargin{i+1})~=varargin{i+1}) ...
                        || any(varargin{i+1}<0 | varargin{i+1}>9)))
                    error('seizmo:addarrivals:badInput',...
                        ['FIELDS must be a index array w/ values from\n'...
                        '0 to 9 indicating the kt,t,user fields to be\n'...
                        'overwrote.  List is common to all records!']);
                end
                fields=varargin{i+1};
            otherwise
                error('seizmo:addarrivals:badOption',...
                    'Unknown Option: %s',varargin{i});
        end
    end

    % get relevant header info
    [evdp,gcarc,dist,stla,stlo,evla,evlo,o,t,kt,user]=getheader(data,...
        'evdp','gcarc','dist','stla','stlo',...
        'evla','evlo','o','t','kt','user');
    
    % detail message
    if(verbose)
        disp('Adding Arrival Info to Record(s)');
        print_time_left(0,nrecs);
    end
    
    % loop over records adding info to header
    for i=1:nrecs
        % set progress bar to overwrite
        redraw=false;
        
        % check header info
        if(isnan(evdp(i)))
            %redraw=true;
            warning('seizmo:addarrivals:badEVDP',...
                'Record: %d\nEVDP undefined! Treating as zero.',i);
            evdp(i)=0;
        end
        if(isnan(o(i)))
            %redraw=true;
            warning('seizmo:addarrivals:badO',...
                'Record: %d\nO field undefined! Treating as zero.',i);
            o(i)=0;
        end
        if(~isnan(gcarc(i)))
            location={'deg' gcarc(i)};
        elseif(~isnan(dist(i)))
            location={'km' dist(i)};
        elseif(all(~isnan([stla(i) stlo(i) evla(i) evlo(i)])))
            location={'sta' [stla(i) stlo(i)] 'evt' [evla(i) evlo(i)]};
        else
            error('seizmo:addarrivals:badLocation',...
                ['Record: %d\nGCARC, DIST, or '...
                'STLA+STLO+EVLA+EVLO must be set to get arrivals!'],i)
        end

        % get arrivals
        arrivals=tauptime('mod',model{i},'h',evdp(i)/1000,...
            'ph',phases{i},location{:});

        % add arrivals
        for j=1:min(numel(fields),numel(arrivals))
            t(i,fields(j)+1)=arrivals(j).time+o(i);
            kt{i,fields(j)+1}=arrivals(j).phase;
            user(i,fields(j)+1)=arrivals(j).rayparameter;
        end
        
        % detail message
        if(verbose); print_time_left(i,nrecs,redraw); end
    end

    % update header
    data=changeheader(data,'t',t,'kt',kt,'user',user);
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

end
