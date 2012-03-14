function [win,peak,bad]=autowindow(data,thresh,reach,env,warn)
%AUTOWINDOW    Automatic windowing of SEIZMO records
%
%    Usage:    [win,peak]=autowindow(data,thresh,reach)
%              [win,peak]=autowindow(data,thresh,reach,env)
%              [win,peak,bad]=autowindow(data,thresh,reach,env,warn)
%
%    Description:
%     [WIN,PEAK]=AUTOWINDOW(DATA,THRESH,REACH) determines a window around
%     the strongest peak in each record in SEIZMO struct DATA.  The window
%     is REACH times the distance to the first points from the peak that
%     are below THRESH times the peak value.  The records are transformed
%     to the envelope of their analytic signal before the window
%     determination to aid in windowing oscillatory signals.  If the peak
%     is on the edge of the record or the window extends past the record
%     limits, a warning is issued and the window is reset to extend only to
%     the record time limits.  WIN is in relative time and is a NRECSx2
%     matrix of [START STOP].  PEAK gives the time of the peak in the
%     windows. Spectral records are not supported.
%
%     [WIN,PEAK]=AUTOWINDOW(DATA,THRESH,REACH,ENV) sets if the records are
%     enveloped before the window determination.  This is useful if the
%     records are already envelopes (ie for stacking).  The default is
%     TRUE, which means it will calculate the envelope.
%
%     [WIN,PEAK,BAD]=AUTOWINDOW(DATA,THRESH,REACH,ENV,WARN) alters the
%     behavior in response to edge conditions (peak at 1st or last point
%     in record, window extending past 1st or last point).  WARN=TRUE is
%     the default.  BAD is a logical array indicating the records that had
%     edge issues.
%
%    Notes:
%
%    Examples:
%     % Get an automatic window and implement it:
%     win=autowindow(data,0.3,1.5);
%     windata=cut(data,win(:,1),win(:,2));
%     plot0(windata);
%
%    See also: USERWINDOW, ENVELOPE, CUT

%     Version History:
%        May  20, 2010 - initial version
%        Nov.  2, 2011 - doc update
%        Mar. 13, 2012 - seizmocheck fix, use getheader improvements
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 13, 2012 at 18:50 GMT

% todo:

% check nargin
error(nargchk(3,5,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check headers
    data=checkheader(data,...
        'MULCMP_DEP','ERROR',...
        'NONTIME_IFTYPE','ERROR');
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
end

% do auto window
try
    % verbosity
    verbose=seizmoverbose;
    
    % number of records
    nrecs=numel(data);

    % default options
    if(nargin<4 || isempty(env)); env=true; end
    if(nargin<5 || isempty(warn)); warn=true; end

    % check options
    if(~isreal(thresh) || (~isscalar(thresh) && numel(thresh)~=nrecs))
        error('seizmo:autowindow:badInput',...
            ['THRESH must be a real-valued scalar ' ...
            'or array with NRECS elements in the range 0 to 1!']);
    elseif(~isreal(reach) || (~isscalar(reach) && numel(reach)~=nrecs))
        error('seizmo:autowindow:badInput',...
            ['REACH must be a real-valued scalar ' ...
            'or array with NRECS elements!']);
    elseif(~islogical(env) || (~isscalar(env) && numel(env)~=nrecs))
        error('seizmo:autowindow:badInput',...
            ['ENV must be a logical scalar ' ...
            'or array with NRECS elements!']);
    elseif(~islogical(warn) || ~isscalar(warn))
        error('seizmo:autowindow:badInput',...
            'WARN must be a logical scalar!');
    end
    
    % expand scalars
    if(isscalar(thresh)); thresh(1:nrecs,1)=thresh; end
    if(isscalar(reach)); reach(1:nrecs,1)=reach; end
    if(isscalar(env)); env(1:nrecs,1)=env; end
    
    % header info
    [b,e,delta,npts,leven]=getheader(data,...
        'b','e','delta','npts','leven lgc');
    leven=~strcmpi(leven,'false');
    
    % envelope those desired
    if(any(env)); data(env)=envelope(data(env)); end
    
    % detail message
    if(verbose)
        disp('Determining Automatic Window Limits for Record(s)');
        print_time_left(0,nrecs);
    end
    
    % loop over records
    win=nan(nrecs,2); bad=false(nrecs,1); peak=nan(nrecs,1);
    for i=1:nrecs
        % handle dataless
        if(npts(i)==0);
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            bad(i)=true;
            continue;
        end
        
        % get peak
        [pval,pidx]=max(data(i).dep);
        
        % get peak time
        if(leven(i))
            peak(i)=b(i)+(pidx-1)*delta(i);
        else % uneven
            peak(i)=data(i).idx(pidx);
        end
        
        % is peak on edge?
        skip1=false; skip2=false;
        if(pidx==1)
            skip1=true;
            win(i,1)=b(i);
            bad(i)=true;
        end
        if(pidx==npts(i))
            skip2=true;
            win(i,2)=e(i);
            bad(i)=true;
        end
        
        % get window limits
        if(~skip1)
            tidx1=find(data(i).dep(1:pidx-1)<=pval*thresh(i),1,'last');
            if(isempty(tidx1))
                if(leven(i))
                    win(i,1)=b(i);
                else % uneven
                    win(i,1)=data(i).ind(1);
                end
                bad(i)=true;
            end
            if(leven(i))
                win(i,1)=peak(i)...
                    -reach(i)*(peak(i)-(b(i)+(tidx1-1)*delta(i)));
                if(win(i,1)<b(i))
                    win(i,1)=b(i);
                    bad(i)=true;
                end
            else % uneven
                win(i,1)=peak(i)-reach(i)*(peak(i)-data(i).idx(tidx1));
                if(win(i,1)<data(i).ind(1))
                    win(i,1)=data(i).ind(1);
                    bad(i)=true;
                end
            end
        end
        if(~skip2)
            tidx2=find(data(i).dep(pidx+1:end)<=pval*thresh(i),1,'first');
            if(isempty(tidx2))
                if(leven(i))
                    win(i,2)=e(i);
                else % uneven
                    win(i,2)=data(i).ind(end);
                end
                bad(i)=true;
            else
                tidx2=tidx2+pidx;
                if(leven(i))
                    win(i,2)=peak(i)...
                        +reach(i)*((b(i)+(tidx2-1)*delta(i))-peak(i));
                    if(win(i,2)>e(i))
                        win(i,2)=e(i);
                        bad(i)=true;
                    end
                else % uneven
                    win(i,2)=peak(i)+reach(i)*(data(i).idx(tidx2)-peak(i));
                    if(win(i,2)<data(i).ind(end))
                        win(i,2)=data(i).ind(end);
                        bad(i)=true;
                    end
                end
            end
        end
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end
    
    % warn about problems
    if(warn && any(bad))
        warning('seizmo:autowindow:edgeIssues',...
            ['Edge issues encountered for record(s):\n' ...
            sprintf('%d ',find(bad))]);
    end
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror);
end

end
