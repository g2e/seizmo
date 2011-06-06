function []=globalcmt_update()
%GLOBALCMT_UPDATE    Updates GlobalCMT catalogs (requires internet)
%
%    Usage:    globalcmt_update
%
%    Description:
%     GLOBALCMT_UPDATE will search the GlobalCMT Project's website for
%     updates to their catalogs and will add any new CMTs found to SEIZMO's
%     catalogs.  GLOBALCMT_UPDATE does not check for changes in the old
%     catalogs from the GlobalCMT Project (old being those that are already
%     a part of the SEIZMO catalogs).  Try not to use GLOBALCMT_UPDATE too
%     often as it downloads & updates the quick CMT catalog every run.
%
%    Notes:
%     - needs write permission to SEIZMO directories
%     - also updates the cached catalogs under SEIZMO.GLOBALCMT
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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan.  5, 2011 at 21:30 GMT

% todo:

% verbosity
verbose=seizmoverbose;

% SEIZMO global access
global SEIZMO

% load full globalcmt catalog
try
    full=SEIZMO.GLOBALCMT.FULL;
catch
    full=load('globalcmt_full.mat');
    SEIZMO.GLOBALCMT.FULL=full;
end
fields=fieldnames(full);
nf=numel(full.name);

% get latest month in full
lastyr=max(full.year);
lastmon=max(full.month(full.year==lastyr));
lastyr=lastyr-2000;

% save for quick cmts
final=[lastyr lastmon];

% get month before current
time=datevec(now);
time(2)=time(2)-1;
time(3)=1;
time=fixdates(time(1:3));
maxyr=time(1)-2000;
maxmon=time(2);

% month strings
month={'jan' 'feb' 'mar' 'apr' 'may' 'jun' ...
    'jul' 'aug' 'sep' 'oct' 'nov' 'dec'};

% update full
skip=false; updated=false;
url='http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_MONTHLY/';
if(lastyr==maxyr)
    for j=lastmon:maxmon
        [ndk,ok]=urlread([url '/' num2str(maxyr+2000) '/' ...
            month{j} num2str(maxyr,'%02d') '.ndk']);
        if(~ok || isempty(ndk))
            break;
        else
            % detail message
            if(verbose)
                disp(['Retrieved ' month{j} num2str(maxyr,'%02d') '.ndk']);
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
            final=[maxyr j];
        end
    end
else
    for i=lastyr:maxyr
        if(skip); break; end
        if(i==lastyr)
            for j=lastmon:12
                [ndk,ok]=urlread([url '/' num2str(i+2000) '/' ...
                    month{j} num2str(i,'%02d') '.ndk']);
                if(~ok || isempty(ndk))
                    skip=true;
                    break;
                else
                    % detail message
                    if(verbose)
                        disp(['Retrieved ' month{j} ...
                            num2str(i,'%02d') '.ndk']);
                    end
                    
                    % updated flag
                    updated=true;
                    
                    % convert to struct
                    cmt=readndk(ndk,true);
                    
                    % combine with full
                    for k=1:numel(fields)
                        full.(fields{k})=[full.(fields{k}); ...
                            cmt.(fields{k})];
                    end
                    
                    % update final
                    final=[i j];
                end
            end
        elseif(i==maxyr)
            for j=1:maxmon
                [ndk,ok]=urlread([url '/' num2str(i+2000) '/' ...
                    month{j} num2str(i,'%02d') '.ndk']);
                if(~ok || isempty(ndk))
                    skip=true;
                    break;
                else
                    % detail message
                    if(verbose)
                        disp(['Retrieved ' month{j} ...
                            num2str(i,'%02d') '.ndk']);
                    end
                    
                    % updated flag
                    updated=true;
                    
                    % convert to struct
                    cmt=readndk(ndk,true);
                    
                    % combine with full
                    for k=1:numel(fields)
                        full.(fields{k})=[full.(fields{k}); ...
                            cmt.(fields{k})];
                    end
                    
                    % update final
                    final=[i j];
                end
            end
        else
            for j=1:12
                [ndk,ok]=urlread([url '/' num2str(i+2000) '/' ...
                    month{j} num2str(i,'%02d') '.ndk']);
                if(~ok || isempty(ndk))
                    skip=true;
                    break;
                else
                    % detail message
                    if(verbose)
                        disp(['Retrieved ' month{j} ...
                            num2str(i,'%02d') '.ndk']);
                    end
                    
                    % updated flag
                    updated=true;
                    
                    % convert to struct
                    cmt=readndk(ndk,true);
                    
                    % combine with full
                    for k=1:numel(fields)
                        full.(fields{k})=[full.(fields{k}); ...
                            cmt.(fields{k})];
                    end
                    
                    % update final
                    final=[i j];
                end
            end
        end
    end
end

% only update if new cmts found
if(updated)
    % delete dupes
    [keep,keep]=unique(full.name);
    full=ssidx(full,sort(keep));
    
    % detail message
    if(verbose)
        nf2=numel(full.name);
        disp(['Found ' num2str(nf2-nf) ' CMTs']);
    end
    
    % save full
    path=fileparts(mfilename('fullpath'));
    SEIZMO.GLOBALCMT.FULL=full;
    save([path filesep 'globalcmt_full.mat'],'-struct','full');
else
    % detail message
    if(verbose)
        nf2=numel(full.name);
        disp(['Found ' num2str(nf2-nf) ' CMTs']);
    end
end

% get quick catalog
qcmt='http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_QUICK/';
[qndk,ok]=urlread([qcmt 'qcmt.ndk']);

% skip if could not read
if(ok && ~isempty(qndk))
    % fix final
    final(1)=final(1)+2000;
    final=fixdates([final(1) final(2)+1 1]);
    
    % convert to struct
    quick=readndk(qndk,true);
    
    % remove dupes
    keep=(quick.year==final(1) & quick.month>=final(2)) ...
        | quick.year>final(1);
    quick=ssidx(quick,keep);
    
    % detail message
    nq=numel(quick.name);
    if(verbose); disp(['Found ' num2str(nq) ' Quick CMTs']); end
    
    % save
    path=fileparts(mfilename('fullpath'));
    SEIZMO.GLOBALCMT.QUICK=quick;
    save([path filesep 'globalcmt_quick.mat'],'-struct','quick');
end

end
