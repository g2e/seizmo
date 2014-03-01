function []=globalcmt_update()
%GLOBALCMT_UPDATE    Updates GlobalCMT catalogs (requires internet)
%
%    Usage:    globalcmt_update
%
%    Description:
%     GLOBALCMT_UPDATE will search the GlobalCMT Project's website for
%     additions to their catalogs and will add any new CMTs found to
%     SEIZMO's local catalogs.  GLOBALCMT_UPDATE does not check for changes
%     in catalogs previously downloaded from the GlobalCMT Project (the
%     exception is the quick CMT catalog which is redownloaded every time
%     GLOBALCMT_UPDATE is run).
%
%    Notes:
%     - Try not to use GLOBALCMT_UPDATE too often as it downloads & updates
%       the entire quick CMT catalog every run.
%     - Needs write permission to SEIZMO directories.
%     - Also updates the cached catalogs under SEIZMO.GLOBALCMT
%
%    Examples:
%     % Update your catalog, then find CMTs from the last week:
%     globalcmt_update
%     findcmts('st',datevec(now-7),'nd',7)
%
%    See also: READNDK, FINDCMTS, FINDCMT, SSIDX, SETEVENT

%     Version History:
%        Aug.  3, 2010 - initial version
%        Jan.  5, 2011 - improved docs, fixed download bug
%        Nov.  1, 2011 - condensed code to remove some redundancies
%        Feb. 15, 2012 - say 'new cmts' to avoid confusion
%        Mar.  1, 2012 - minor doc update, Octave save workaround
%        Mar. 20, 2013 - no error if cannot write (just warn and move on)
%        Jan. 14, 2014 - catches urlread errors and gives warning
%        Jan. 25, 2014 - indented verbose messages for visual improvement
%        Feb. 27, 2014 - catch error if no catalog, doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 27, 2014 at 21:30 GMT

% todo:

% verbosity
verbose=seizmoverbose;

% SEIZMO global access
global SEIZMO

% load full globalcmt catalog
try
    full=SEIZMO.GLOBALCMT.FULL;
catch
    try
        full=load('globalcmt_full.mat');
        SEIZMO.GLOBALCMT.FULL=full;
    catch
        warning('seizmo:globalcmt_update:noCatalog',...
            'GLOBALCMT catalog does not exist!  Run GLOBALCMT_CREATE!');
        return;
    end
end
fields=fieldnames(full);
nf=numel(full.name);

% get latest month in full
lastyr=max(full.year);
lastmon=max(full.month(full.year==lastyr));

% save for quick cmts
final=[lastyr lastmon];

% get month before current
time=datevec(now);
time(2)=time(2)-1;
time(3)=1;
time=fixdates(time(1:3));
maxyr=time(1);
maxmon=time(2);

% month strings
month={'jan' 'feb' 'mar' 'apr' 'may' 'jun' ...
    'jul' 'aug' 'sep' 'oct' 'nov' 'dec'};

% update full
skip=false; updated=false;
url='http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_MONTHLY/';
for i=lastyr:maxyr
    % early exit after last monthly catalog found
    if(skip); break; end
    
    % get months to go through
    if(lastyr==maxyr)
        mon=lastmon:maxmon;
    elseif(i==lastyr)
        mon=lastmon:12;
    elseif(i==maxyr)
        mon=1:maxmon;
    else
        mon=1:12;
    end
    
    % loop over months of this year
    for j=mon
        % get catalog
        try
            [ndk,ok]=urlread([url '/' num2str(i) '/' ...
                month{j} num2str(i-2000,'%02d') '.ndk']);
        catch % no internet connection?
            le=lasterror;
            warning(le.identifier,le.message);
            skip=true;
            break;
        end
        
        % check that file exists and has entries
        % NOTE: this fails if the monthly catalog was actually empty
        if(~ok || isempty(ndk))
            skip=true;
            break;
        else
            % detail message
            if(verbose)
                disp([' Retrieved ' month{j} ...
                    num2str(i-2000,'%02d') '.ndk']);
            end
            
            % updated flag
            updated=true;
            
            % convert to struct
            cmt=readndk(ndk,true);
            
            % combine with full
            for k=1:numel(fields)
                full.(fields{k})=[full.(fields{k}); cmt.(fields{k})];
            end
            
            % update final
            final=[i j];
        end
    end
end

% only update if new cmts found
if(updated)
    % delete dupes by name (also sorts by name)
    [keep,keep]=unique(full.name);
    full=ssidx(full,sort(keep));
    
    % detail message
    if(verbose)
        nf2=numel(full.name);
        disp([' Found ' num2str(nf2-nf) ' New CMTs']);
    end
    
    % save full
    path=fileparts(mfilename('fullpath'));
    SEIZMO.GLOBALCMT.FULL=full;
    try
        if(isoctave)
            save([path filesep 'globalcmt_full.mat'],...
                '-struct','-7','full');
        else % matlab
            save([path filesep 'globalcmt_full.mat'],'-struct','full');
        end
    catch
        le=lasterror;
        warning(le.identifier,le.message);
    end
else
    % detail message
    if(verbose); disp(' Found 0 CMTs'); end
end

% get quick catalog
qcmt='http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_QUICK/';
try
    [qndk,ok]=urlread([qcmt 'qcmt.ndk']);
catch % no internet connection?
    le=lasterror;
    warning(le.identifier,le.message);
    ok=false;
    qndk=[];
end

% skip if could not read
if(ok && ~isempty(qndk))
    % fix final
    final=fixdates([final(1) final(2)+1 1]);
    
    % convert to struct
    quick=readndk(qndk,true);
    
    % remove dupes
    keep=(quick.year==final(1) & quick.month>=final(2)) ...
        | quick.year>final(1);
    quick=ssidx(quick,keep);
    
    % detail message
    nq=numel(quick.name);
    if(verbose); disp([' Found ' num2str(nq) ' New Quick CMTs']); end
    
    % save
    path=fileparts(mfilename('fullpath'));
    SEIZMO.GLOBALCMT.QUICK=quick;
    try
        if(isoctave)
            save([path filesep 'globalcmt_quick.mat'],...
                '-struct','-7','quick');
        else % matlab
            save([path filesep 'globalcmt_quick.mat'],'-struct','quick');
        end
    catch
        le=lasterror;
        warning(le.identifier,le.message);
    end
end

end
