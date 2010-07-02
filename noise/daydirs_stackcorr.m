function []=daydirs_stackcorr(indir,outdir,stns,cmp1,cmp2,o)
%DAYDIRS_STACKCORR    Stacks correlograms in day directories
%
%    Usage:    daydirs_stackcorr(indir,outdir,stns,cmp1,cmp2)
%              daydirs_stackcorr(indir,outdir,stns,cmp1,cmp2,overwrite)
%
%    Description: DAYDIRS_STACKCORR(INDIR,OUTDIR,STNS,CMP1,CMP2) stacks
%     correlations between stations with names given in char/cellstr array
%     STNS and component names given by CMP1 & CMP2.  CMP1 & CMP2 must be a
%     3 character component code like 'LHZ'.  There are several stack sets
%     computed: ones for the full time, individual months and year-
%     independent monthly stacks (for example Jan 2006 & 2007 are stacked
%     together).  In addition, the stacks are returned as the raw
%     correlograms or as empirical Green's functions (which is further
%     broken down into the positive component, negative component,
%     symmetric component and normal two-way Green's function).  As with
%     all daydirs functions this operates on records within a day directory
%     layout under INDIR, returning the result in OUTDIR.  Note that the
%     directory layout in OUTDIR is quite different:
%      OUTDIR
%       |
%       +-> XX
%            |
%            +-> CORR_FULLSTACK
%                CORR_1MONSTACK
%                CORR_MONSTACK
%                GREEN_FULLSTACK
%                GREEN_FULLSTACK_NEGCMP
%                GREEN_FULLSTACK_POSCMP
%                GREEN_FULLSTACK_SYMCMP
%                GREEN_1MONSTACK
%                GREEN_1MONSTACK_NEGCMP
%                GREEN_1MONSTACK_POSCMP
%                GREEN_1MONSTACK_SYMCMP
%                GREEN_MONSTACK
%                GREEN_MONSTACK_NEGCMP
%                GREEN_MONSTACK_POSCMP
%                GREEN_MONSTACK_SYMCMP
%     
%     where XX is the last char of CMP1 & CMP2 combined (so 'LHZ' 'LHZ'
%     would produce 'ZZ').
%
%     DAYDIRS_STACKCORR(INDIR,OUTDIR,STNS,CMP1,CMP2,OVERWRITE) quietly
%     overwrites pre-existing records in OUTDIR when OVERWRITE is set to
%     TRUE.  By default OVERWRITE is FALSE.
%
%    Notes:
%
%    Header changes: SCALE (number of records in stack), DEP*
%
%    Examples:
%
%    See also: DAYDIRS_MERGECUT_25HRS, DAYDIRS_RESAMPLE, DAYDIRS_NORMALIZE,
%              DAYDIRS_CORRELATE, DAYDIRS_ROTCORR, DAYDIRS_RINST,
%              DAYDIRS_MAKE

%     Version History:
%        June 20, 2010 - added to seizmo, fixed bug that would replace
%                        correlograms rather than combine with reversed
%                        ones, added docs
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 20, 2010 at 12:55 GMT

% todo:

% check nargin
error(nargchk(5,6,nargin));

% defaults
if(nargin<6 || isempty(o)); o=false; end
if(~isscalar(o) || ~islogical(o))
    error('seizmo:daydirs_rotcorr:badInput',...
        'OVERWRITE flag must be a scalar logical!');
end

% check directories
if(~ischar(indir) || ~isvector(indir))
    error('seizmo:daydirs_stackcorr:fileNotString',...
        'INDIR must be a string!');
end
if(~exist(indir,'dir'))
    error('seizmo:daydirs_stackcorr:dirConflict',...
        ['Input Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(~ischar(outdir) || ~isvector(outdir))
    error('seizmo:daydirs_stackcorr:fileNotString',...
        'OUTDIR must be a string!');
end
if(exist(outdir,'file'))
    if(~exist(outdir,'dir'))
        error('seizmo:daydirs_stackcorr:dirConflict',...
            'Output Directory: %s\nIs a file!',outdir);
    end
    if(~o)
        fprintf('Output Directory: %s\nDirectory Exists!\n',outdir);
        reply=input('Overwrite? Y/N [N]: ','s');
        if(isempty(reply) || ~strncmpi(reply,'y',1))
            disp('Not overwriting!');
            return;
        end
        disp('Overwriting!');
    end
end

% check stns & components
if(ischar(stns)); stns=cellstr(stns); end
if(~iscellstr(stns) || any(cellfun('ndims',stns)~=2) ...
        || any(cellfun('size',stns,1)~=1))
    error('seizmo:daydirs_stackcorr:badInput',...
        'STNS must be a char/cellstr array of station names!');
end
if(~ischar(cmp1) || ~isequal(size(cmp1),[1 3]))
    error('seizmo:daydirs_stackcorr:badInput',...
        'CMP1 must be a 3 character component code!');
end
if(~ischar(cmp2) || ~isequal(size(cmp2),[1 3]))
    error('seizmo:daydirs_stackcorr:badInput',...
        'CMP2 must be a 3 character component code!');
end

% directory separator
fs=filesep;

% parallel processing setup (8 instances)
matlabpool(8);

% number of stations
nstn=numel(stns);

% get orientation codes for output
cname=upper([cmp1(3) cmp2(3)]);

% alter verbosity
verbose=seizmoverbose(false);
if(verbose); disp('Stacking Correlograms'); end

for m=1:nstn
    parfor s=1:nstn
        % skip if equal
        if(m==s); continue; end
        
        % where are we?
        disp([stns{m} '.' cmp1 ' <-> ' stns{s} '.' cmp2])
        
        % read correlograms handling if no records found
        skip=false;
        try
            data=readseizmo(...
                [indir fs '*' fs '*' fs 'CORR_-_*' stns{m} '*' cmp1 '*' ...
                stns{s} '*' cmp2]);
        catch
            skip=true;
            try
                data=reverse_correlations(readseizmo(...
                    [indir fs '*' fs '*' fs 'CORR_-_*' stns{s} '*' ...
                    cmp2 '*' stns{m} '*' cmp1]));
            catch
                continue;
            end
        end
        if(~skip)
            try
                data=[data; ...
                    reverse_correlations(readseizmo(...
                    [indir fs '*' fs '*' fs 'CORR_-_*' stns{s} ...
                    '*' cmp2 '*' stns{m} '*' cmp1]))];
            catch
                % no reverse
            end
        end
        
        % skip if none
        nrecs=numel(data);
        if(verbose); disp(['FOUND ' num2str(nrecs) ' RECORDS']); end
        if(~nrecs); continue; end
        
        % get year and month
        [nzyear,nzmonth,nzcday]=getheader(data,...
            'nzyear','nzmonth','nzcday');
        yrs=unique(nzyear); yrs=yrs(:).';
        mons=unique(nzmonth); mons=mons(:).';
        dates=unique([nzyear nzmonth nzcday],'rows');

        % fix user0/1
        data=changeheader(data,'user0',m,'user1',s);
        
        % stack
        if(verbose); disp('STACKING'); end
        binoperr('ref','ignore');
        stack_full=divide(addrecords(data),nrecs);
        stack_full.name=['STACK_FULL_' stns{m} '_' stns{s} '_' cname];
        stack_full=changeheader(stack_full,'z',{[dates(1,:) 0 0 0]},...
            'scale',nrecs);
        count1=0;
        for yr=yrs
            for mon=mons
                try
                    tmp=addrecords(data(nzyear==yr & nzmonth==mon));
                    tmp=divide(tmp,sum(nzyear==yr & nzmonth==mon));
                catch
                    continue;
                end
                count1=count1+1;
                if(count1==1)
                    stack_1mon=tmp;
                else
                    stack_1mon(count1,1)=tmp;
                end
                stack_1mon(count1,1).name=...
                    ['STACK_' num2str(yr) '.' sprintf('%02d',mon) ...
                    '_' stns{m} '_' stns{s} '_' cname];
                stack_1mon(count1,1)=changeheader(stack_1mon(count1,1),...
                    'z',{[yr mon 1 0 0 0]},...
                    'scale',sum(nzyear==yr & nzmonth==mon));
            end
        end
        count2=0;
        for mon=mons
            try
                tmp=addrecords(data(nzmonth==mon));
                tmp=divide(tmp,sum(nzmonth==mon));
            catch
                continue;
            end
            count2=count2+1;
            if(count2==1)
                stack_mon=tmp;
            else
                stack_mon(count2,1)=tmp;
            end
            stack_mon(count2,1).name=['STACK_' sprintf('%02d',mon) ...
                '_' stns{m} '_' stns{s} '_' cname];
            stack_mon(count2,1)=changeheader(stack_mon(count2,1),...
                'z',{[0 mon 1 0 0 0]},'scale',sum(nzmonth==mon));
        end
        
        % detail message
        if(verbose); disp('GETTING COMPONENTS'); end
        
        % negative time derivative gives empirical green's function
        dif_stack_full=seizmofun(stack_full,@(x)(-gradient(x)));
        dif_stack_1mon=seizmofun(stack_1mon,@(x)(-gradient(x)));
        dif_stack_mon=seizmofun(stack_mon,@(x)(-gradient(x)));
        
        % get positive component
        dif_pos_stack_full=cut(dif_stack_full,0);
        dif_pos_stack_1mon=cut(dif_stack_1mon,0);
        dif_pos_stack_mon=cut(dif_stack_mon,0);
        
        % get negative component
        dif_neg_stack_full=cut(reverse(dif_stack_full),0);
        dif_neg_stack_1mon=cut(reverse(dif_stack_1mon),0);
        dif_neg_stack_mon=cut(reverse(dif_stack_mon),0);
        
        % get symmetric component
        dif_sym_stack_full=seizmofun(ch(dif_stack_full,'b',0),...
            @(x)[x(ceil(end/2),:); ...
                 x(floor(end/2):-1:1,:)+x(ceil(end/2)+1:end,:)]);
        dif_sym_stack_1mon=seizmofun(ch(dif_stack_1mon,'b',0),...
            @(x)[x(ceil(end/2),:); ...
                 x(floor(end/2):-1:1,:)+x(ceil(end/2)+1:end,:)]);
        dif_sym_stack_mon=seizmofun(ch(dif_stack_mon,'b',0),...
            @(x)[x(ceil(end/2),:); ...
                 x(floor(end/2):-1:1,:)+x(ceil(end/2)+1:end,:)]);
        
        % write it all
        if(verbose); disp('WRITING'); end
        w(stack_full,'path',[outdir fs cname fs 'CORR_FULLSTACK' fs]);
        w(stack_1mon,'path',[outdir fs cname fs 'CORR_1MONSTACK' fs]);
        w(stack_mon,'path',[outdir fs cname fs 'CORR_MONSTACK' fs]);
        w(dif_stack_full,'path',...
            [outdir fs cname fs 'GREEN_FULLSTACK' fs],'prepend','GREEN_');
        w(dif_stack_1mon,'path',...
            [outdir fs cname fs 'GREEN_1MONSTACK' fs],'prepend','GREEN_');
        w(dif_stack_mon,'path',...
            [outdir fs cname fs 'GREEN_MONSTACK' fs],'prepend','GREEN_');
        w(dif_pos_stack_full,...
            'path',[outdir fs cname fs 'GREEN_FULLSTACK_POSCMP' fs],...
            'prepend','POSCMP_GREEN_');
        w(dif_pos_stack_1mon,...
            'path',[outdir fs cname fs 'GREEN_1MONSTACK_POSCMP' fs],...
            'prepend','POSCMP_GREEN_');
        w(dif_pos_stack_mon,...
            'path',[outdir fs cname fs 'GREEN_MONSTACK_POSCMP' fs],...
            'prepend','POSCMP_GREEN_');
        w(dif_neg_stack_full,...
            'path',[outdir fs cname fs 'GREEN_FULLSTACK_NEGCMP' fs],...
            'prepend','NEGCMP_GREEN_');
        w(dif_neg_stack_1mon,...
            'path',[outdir fs cname fs 'GREEN_1MONSTACK_NEGCMP' fs],...
            'prepend','NEGCMP_GREEN_');
        w(dif_neg_stack_mon,...
            'path',[outdir fs cname fs 'GREEN_MONSTACK_NEGCMP' fs],...
            'prepend','NEGCMP_GREEN_');
        w(dif_sym_stack_full,...
            'path',[outdir fs cname fs 'GREEN_FULLSTACK_SYMCMP' fs],...
            'prepend','SYMCMP_GREEN_');
        w(dif_sym_stack_1mon,...
            'path',[outdir fs cname fs 'GREEN_1MONSTACK_SYMCMP' fs],...
            'prepend','SYMCMP_GREEN_');
        w(dif_sym_stack_mon,...
            'path',[outdir fs cname fs 'GREEN_MONSTACK_SYMCMP' fs],...
            'prepend','SYMCMP_GREEN_');
    end
end

% parallel processing takedown
matlabpool close;
seizmoverbose(verbose);

end
