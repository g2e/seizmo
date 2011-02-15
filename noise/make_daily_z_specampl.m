function [ampl]=make_daily_z_specampl(indir,outdir,startdate,enddate)
%MAKE_DAILY_Z_SPECAMPL    Return daily array spectral amplitude
%
%    Usage:    ampl=make_daily_z_specampl(corr_dir)
%              ampl=make_daily_z_specampl(corr_dir,odir)
%              ampl=make_daily_z_specampl(corr_dir,odir,startdate,enddate)
%
%    Description:
%     AMPL=MAKE_DAILY_Z_SPECAMPL(CORR_DIR) returns spectral amplitudes
%     for the entire array for daily time spans using the correlograms
%     given in CORR_DIR.  CORR_DIR must have been created by
%     DAYDIRS_CORRELATE.
%
%     AMPL=MAKE_DAILY_Z_SPECAMPL(CORR_DIR,ODIR) will also save the records
%     to the directory given by ODIR.
%
%     AMPL=MAKE_DAILY_Z_SPECAMPL(CORR_DIR,ODIR,STARTDATE,ENDDATE)
%     explicitly sets the day range allowed.  The dates must be appropriate
%     for DATENUM translation or of the form [yr jday] or YEAR.DAY where
%     DAY is the thousandths term (ie 2006.365 & 2006.001).
%
%    Notes:
%     - AMPL is a general xy filetype.
%
%    Examples:
%     % Get daily array spectral amplitudes and plot them up as an image:
%     dailyampl=make_daily_z_specampl(my_corr_dir);
%     days=datenum(cell2mat(getheader(dailyampl,'z6')));
%     % need to fill in the day gaps with blanks
%     % create image matrix
%     % plot it using imagesc
%     % use dateticks to make it look good
%
%    See also: MAKE_FULL_Z_SPECAMPL, MAKE_MONTHLY_Z_SPECAMPL,
%              DAYDIRS_CORRELATE, MAKE_YRMO_Z_SPECAMPL

%     Version History:
%        Feb. 14, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 14, 2011 at 15:55 GMT

% todo:

% check nargin
error(nargchk(1,4,nargin));

% directory separator
fs=filesep;

% check stack dir
if(~ischar(indir) || ~isvector(indir))
    error('seizmo:make_daily_z_specampl:fileNotString',...
        'CORR_DIR must be a string!');
end
if(~exist(indir,'dir'))
    error('seizmo:make_daily_z_specampl:dirConflict',...
        ['CORR_DIR Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(nargin>1 && ~isempty(outdir))
    wo=true;
    if(~ischar(outdir) || ~isvector(outdir))
        error('seizmo:make_daily_z_specampl:fileNotString',...
            'OUTDIR must be a string!');
    end
    o=false;
    if(exist(outdir,'file'))
        if(~exist(outdir,'dir'))
            error('seizmo:make_daily_z_specampl:dirConflict',...
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
else
    wo=false;
end

% get year directories and day directories
dirs=xdir([indir fs]);
dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dirs
years=str2double({dirs.name});
nyears=numel(years);
if(any(isnan(years)))
    error('seizmo:make_daily_z_specampl:badLayout',...
        'Improper directory layout!');
end
jdays=cell(size(years));
for i=1:nyears
    % get day directories
    dirs=xdir([indir fs num2str(years(i))]);
    dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dir
    jdays{i}=str2double({dirs.name});
    if(any(isnan(jdays{i})))
        error('seizmo:make_daily_z_specampl:badLayout',...
            'Improper directory layout!');
    end
end

% default date limits
if(nargin<3 || isempty(startdate))
    [ymin,ymini]=min(years); jdmin=min(jdays{ymini});
    startdate=datenum(doy2cal([ymin jdmin]));
end
if(nargin<4 || isempty(enddate))
    [ymax,ymaxi]=max(years); jdmax=max(jdays{ymaxi});
    enddate=datenum(doy2cal([ymax jdmax]));
end

% check date limits
if(isscalar(startdate) && startdate<3000)
    % year.jday
    startdate=datenum(doy2cal(...
        [round(startdate) round(1e3*(startdate-round(startdate)))]));
elseif(numel(startdate)==2)
    % [year jday]
    startdate=datenum(doy2cal(startdate));
else
    % something else
    startdate=datenum(startdate);
end
if(~isscalar(startdate) || ~isfinite(startdate))
    error('seizmo:make_daily_z_specampl:badInput',...
        'STARTDATE is not formatted correctly!');
end
if(isscalar(enddate) && enddate<3000)
    % year.jday
    startdate=datenum(doy2cal(...
        [round(enddate) round(1e3*(enddate-round(enddate)))]));
elseif(numel(enddate)==2)
    % [year jday]
    enddate=datenum(doy2cal(enddate));
else
    % something else
    enddate=datenum(enddate);
end
if(~isscalar(enddate) || ~isfinite(enddate))
    error('seizmo:make_daily_z_specampl:badInput',...
        'ENDDATE is not formatted correctly!');
end

% verbosity (turn it off for the loop)
verbose=seizmoverbose(false);
if(verbose); disp('Computing daily array amplitude spectra'); end

% loop over years
cnt=1;
for i=1:nyears
    % working year
    yr=years(i);
    syr=num2str(yr);
    
    % loop over days
    for j=1:numel(jdays{i})
        % working julian day
        jday=jdays{i}(j);
        sjday=num2str(jday,'%03d');
        
        % skip if not in range
        date=datenum(doy2cal([yr jday]));
        if(date<startdate || date>enddate)
            continue;
        end
        
        % detail message
        if(verbose); disp(['PROCESSING DAY ' syr '.' sjday]); end
        
        % read in data
        try
            data=readseizmo(...
                [indir fs syr fs sjday fs 'CORR_-_MASTER_-_*.BHZ_-_SLAVE_-_*.BHZ'],...
                [indir fs syr fs sjday fs 'CORR_-_MASTER_-_*.BHZ_-_SLAVE_-_*.LHZ'],...
                [indir fs syr fs sjday fs 'CORR_-_MASTER_-_*.LHZ_-_SLAVE_-_*.BHZ'],...
                [indir fs syr fs sjday fs 'CORR_-_MASTER_-_*.LHZ_-_SLAVE_-_*.LHZ']);
        catch
            continue;
        end
        
        % weighted average
        ampl(cnt,1)=divide(stack(keepam(dft(data))),numel(data));
        
        % clean headers somewhat
        ampl(cnt,1)=changeheader(ampl(i,1),'st',nan,'ev',nan,'cmp',nan,...
            'kname',nan,'scale',numel(data),'kstnm',[syr '.' sjday]);
        
        % change name
        ampl(cnt,1)=changename(ampl(i,1),'name',...
            ['ARRAY_SPECAMPL_' syr '.' sjday '_ZZ']);
        if(wo)
            ampl(cnt,1)=changepath(ampl(i,1),'path',outdir);
        else
            ampl(cnt,1)=changepath(ampl(i,1),'path','.');
        end
        
        % write if wanted
        if(~nargout && wo); writeseizmo(ampl); end
        
        % increment
        if(nargout); cnt=cnt+1; end
    end
end

% write if wanted
if(nargout && wo); writeseizmo(ampl); end

% return verbosity to normal
seizmoverbose(verbose);

end
