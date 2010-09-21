function []=daydirs_mergecut_25hrs(indir,outdir,sec,per,tap,o)
%DAYDIRS_MERGECUT_25HRS    Creates 25 hour records from a day directories
%
%    Usage:    daydirs_mergecut_25hrs(indir,outdir)
%              daydirs_mergecut_25hrs(indir,outdir,secs)
%              daydirs_mergecut_25hrs(indir,outdir,secs,min%)
%              daydirs_mergecut_25hrs(indir,outdir,secs,min%,taper%)
%              daydirs_mergecut_25hrs(indir,outdir,secs,min%,taper%,overwr)
%
%    Description: DAYDIRS_MERGECUT_25HRS(INDIR,OUTDIR) cuts and merges
%     2-day sets to create 25-hour records.  This 1 hour overlap is
%     intended to assure 100% usage of the data in the correlation stage.
%     INDIR is the input directory under which is a day directory layout of
%     records.  The processed records are output in an equivalent layout
%     under OUTDIR.
%
%     DAYDIRS_MERGECUT_25HRS(INDIR,OUTDIR,SECS) adjusts the maximum length
%     of records to SECS number of seconds.  The default is 90000 (25
%     hours).
%
%     DAYDIRS_MERGECUT_25HRS(INDIR,OUTDIR,SECS,MIN%) changes the minimum
%     percent of the maximum length a record must be to not be removed.
%     The default MIN% is 70.
%
%     DAYDIRS_MERGECUT_25HRS(INDIR,OUTDIR,SECS,MIN%,TAPER%) sets what
%     percent of the record to taper after merge, cut & trend removal.  The
%     default is 1 (1% of 25 hours is 15 minutes).
%
%     DAYDIRS_MERGECUT_25HRS(INDIR,OUTDIR,SECS,MIN%,TAPER%,OVERWR) quietly
%     overwrites pre-existing records in OUTDIR when OVERWRITE is set to
%     TRUE.  By default OVERWRITE is FALSE.
%
%    Notes:
%
%    Header changes: NPTS, B, E, DEP*
%
%    Examples:
%
%    See also: DAYDIRS_MAKE, DAYDIRS_RINST, DAYDIRS_NORMALIZE,
%              DAYDIRS_CORRELATE, DAYDIRS_ROTCORR, DAYDIRS_STACKCORR,
%              DAYDIRS_RESAMPLE

%     Version History:
%        June 18, 2010 - initial version
%        June 30, 2010 - bugfixes
%        Sep. 21, 2010 - commented out parallel processing lines
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 21, 2010 at 11:15 GMT

% todo:

% check nargin
error(nargchk(2,6,nargin));

% defaults
if(nargin<3 || isempty(sec)); sec=90000; end
if(nargin<4 || isempty(per)); per=70; end
if(nargin<5 || isempty(tap)); tap=1; end
if(nargin<6 || isempty(o)); o=false; end
if(~isscalar(sec) || ~isreal(sec) || sec<=0)
    error('seizmo:daydirs_mergecut_25hrs:badInput',...
        'SECONDS must be a positive real-valued scalar!');
end
if(~isscalar(per) || ~isreal(per) || per<0 || per>100)
    error('seizmo:daydirs_mergecut_25hrs:badInput',...
        'MIN%% must be a real-valued scalar in the range 0 to 100!');
end
if(~isscalar(tap) || ~isreal(tap) || tap<0 || tap>50)
    error('seizmo:daydirs_mergecut_25hrs:badInput',...
        'TAPER%% must be a real-valued scalar in the range 0 to 50!');
end
if(~isscalar(o) || ~islogical(o))
    error('seizmo:daydirs_mergecut_25hrs:badInput',...
        'OVERWR flag must be a scalar logical!');
end

% percent to fraction
per=per/100;
tap=tap/100;

% check directories
if(~ischar(indir) || ~isvector(indir))
    error('seizmo:daydirs_mergecut_25hrs:fileNotString',...
        'INDIR must be a string!');
end
if(~exist(indir,'dir'))
    error('seizmo:daydirs_mergecut_25hrs:dirConflict',...
        ['Input Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(~ischar(outdir) || ~isvector(outdir))
    error('seizmo:daydirs_mergecut_25hrs:fileNotString',...
        'OUTDIR must be a string!');
end
if(exist(outdir,'file'))
    if(~exist(outdir,'dir'))
        error('seizmo:daydirs_mergecut_25hrs:dirConflict',...
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

% minimum number of seconds
minsec=per*sec;

% directory separator
fs=filesep;

% parallel processing setup (8 instances)
%matlabpool(8);

% get year directories and day directories
dirs=xdir([indir fs]);
dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dirs
years=str2double({dirs.name});
nyears=numel(years);
if(any(isnan(years)))
    error('seizmo:daydirs_mergecut_25hrs:badLayout',...
        'Improper directory layout!');
end
jdays=cell(size(years));
for i=1:nyears
    % get day directories
    dirs=xdir([indir fs num2str(years(i))]);
    dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dir
    jdays{i}=str2double({dirs.name});
    if(any(isnan(jdays{i})))
        error('seizmo:daydirs_mergecut_25hrs:badLayout',...
            'Improper directory layout!');
    end
end

% verbosity (turn it off for the loop)
verbose=seizmoverbose(false);
if(verbose); disp('Merging and cutting record(s)'); end

% loop over years
for i=1:nyears
    % working year
    yr=years(i);
    syr=num2str(yr);
    
    % loop over days
    for j=1:numel(jdays{i})
    %parfor j=1:numel(jdays{i})
        % working julian day
        jday=jdays{i}(j);
        sjday=num2str(jday,'%03d');
        
        % get next day (using fixdates)
        oneday=false;
        nxt=fixdates([yr jday+1]);
        syr2=num2str(nxt(1));
        sjday2=num2str(nxt(2),'%03d');
        if(~any(nxt(1)==years) || ~any(nxt(2)==jdays{nxt(1)==years}))
            % nothing the next day so set special no merge flag
            oneday=true;
        end
        
        % detail message
        if(verbose); disp(['PROCESSING DAY ' syr '.' sjday]); end
        
        % attempt merge & cut
        try
            % read in data
            if(verbose); disp(['DAY 1: ' indir fs syr fs sjday fs]); end
            data1=readseizmo([indir fs syr fs sjday fs]);
            if(~oneday)
                if(verbose); disp(['DAY 2: ' indir fs syr2 fs sjday2 fs]); end
                data2=readseizmo([indir fs syr2 fs sjday2 fs]);
                
                % adjust ref time to match 1st day
                data2=synchronize(data2,[yr jday 0 0 0],[],'iday');
            end
            
            % remove dead records
            data1=removedeadrecords(data1);
            if(~oneday); data2=removedeadrecords(data2); end
            
            % skip if none left
            if(oneday)
                if(isempty(data1)); continue; end
            else
                if(isempty(data1) && isempty(data2)); continue; end
                
                % convert to oneday if one is empty
                if(isempty(data1))
                    oneday=true;
                    data1=data2;
                    data2=[];
                elseif(isempty(data2))
                    oneday=true;
                end
            end
            
            % merge datasets
            % - allow for 1 second gaps/overlaps
            % - always adjust later record to that of first
            if(oneday)
                data=merge(data1,'tolerance',1.1,'adjust','last');
            else
                data=merge([data1; data2],...
                    'tolerance',1.1,'adjust','last');
            end
            
            % window records to some number of seconds (default is 25hrs)
            data=cut(data,0,sec);
            
            % skip if none left
            if(isempty(data)); continue; end
            
            % remove records with less than 70% of 25 hours (17.5)
            [b,e]=getheader(data,'b','e');
            data(e-b<minsec)=[];
            
            % skip if none left
            if(isempty(data)); continue; end
            
            % remove trend and taper (1% of 25hrs == 15min)
            data=taper(removetrend(data),tap);
            
            % write
            writeseizmo(data,'pathchange',{indir outdir});
        catch
            % close pool & fix verbosity
            disp([syr '.' sjday])
            tmp=lasterror;
            warning(tmp.message);
            %matlabpool close;
            seizmoverbose(verbose);
            
            % ???
            error(lasterror);
        end
    end
end

% parallel processing takedown & fix verbosity
%matlabpool close;
seizmoverbose(verbose);

end
