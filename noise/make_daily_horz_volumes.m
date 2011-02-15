function []=make_daily_horz_volumes(indir,startdate,enddate)
%MAKE_DAILY_HORZ_VOLUMES    Computes daily fk volumes for horizontals
%
%    Usage:    make_daily_horz_volumes(rotcorr_dir)
%              make_daily_horz_volumes(rotcorr_dir,startdate,enddate)
%
%    Description:
%     MAKE_DAILY_HORZ_VOLUMES(ROTCORR_DIR) creates fk-based slowness
%     response volumes for an array at daily time spans using the
%     correlograms given in ROTCORR_DIR.  ROTCORR_DIR must be created
%     by DAYDIRS_ROTCORR.  The period range is 4 to 100s, the maximum
%     slowness is 50sec/deg and the slowness resolution is 1/3 sec/deg.
%
%     MAKE_DAILY_HORZ_VOLUMES(ROTCORR_DIR,STARTDATE,ENDDATE) explicitly
%     sets the day range allowed.  The dates should be appropriate for
%     DATENUM translation.
%
%    Notes:
%
%    Examples:
%     % 
%
%    See also: MAKE_FULL_HORZ_VOLUMES, MAKE_MONTHLY_HORZ_VOLUMES,
%              DAYDIRS_ROTCORR, MAKE_YRMO_HORZ_VOLUMES, FKXCHORZVOLUME,
%              MAKE_DAILY_Z_VOLUMES

%     Version History:
%        Feb. 14, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 14, 2011 at 15:55 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% directory separator
fs=filesep;

% check stack dir
if(~ischar(indir) || ~isvector(indir))
    error('seizmo:make_daily_horz_volumes:fileNotString',...
        'ROTCORR_DIR must be a string!');
end
if(~exist(indir,'dir'))
    error('seizmo:make_daily_horz_volumes:dirConflict',...
        ['ROTCORR_DIR Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end

% get year directories and day directories
dirs=xdir([indir fs]);
dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dirs
years=str2double({dirs.name});
nyears=numel(years);
if(any(isnan(years)))
    error('seizmo:make_daily_horz_volumes:badLayout',...
        'Improper directory layout!');
end
jdays=cell(size(years));
for i=1:nyears
    % get day directories
    dirs=xdir([indir fs num2str(years(i))]);
    dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dir
    jdays{i}=str2double({dirs.name});
    if(any(isnan(jdays{i})))
        error('seizmo:make_daily_horz_volumes:badLayout',...
            'Improper directory layout!');
    end
end

% default date limits
if(nargin<2 || isempty(startdate))
    [ymin,ymini]=min(years); jdmin=min(jdays{ymini});
    startdate=datenum(doy2cal([ymin jdmin]));
end
if(nargin<3 || isempty(enddate))
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
    error('seizmo:make_daily_horz_volumes:badInput',...
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
    error('seizmo:make_daily_horz_volumes:badInput',...
        'ENDDATE is not formatted correctly!');
end

% verbosity (turn it off for the loop)
verbose=seizmoverbose;
if(verbose); disp('Computing daily horizontal fk volumes'); end

% loop over years
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
            rr=readseizmo(...
                [indir fs syr fs sjday fs 'CORR_-_MASTER_-_*.BHR_-_SLAVE_-_*.BHR'],...
                [indir fs syr fs sjday fs 'CORR_-_MASTER_-_*.BHR_-_SLAVE_-_*.LHR'],...
                [indir fs syr fs sjday fs 'CORR_-_MASTER_-_*.LHR_-_SLAVE_-_*.BHR'],...
                [indir fs syr fs sjday fs 'CORR_-_MASTER_-_*.LHR_-_SLAVE_-_*.LHR']);
            rt=readseizmo(...
                [indir fs syr fs sjday fs 'CORR_-_MASTER_-_*.BHR_-_SLAVE_-_*.BHT'],...
                [indir fs syr fs sjday fs 'CORR_-_MASTER_-_*.BHR_-_SLAVE_-_*.LHT'],...
                [indir fs syr fs sjday fs 'CORR_-_MASTER_-_*.LHR_-_SLAVE_-_*.BHT'],...
                [indir fs syr fs sjday fs 'CORR_-_MASTER_-_*.LHR_-_SLAVE_-_*.LHT']);
            tr=readseizmo(...
                [indir fs syr fs sjday fs 'CORR_-_MASTER_-_*.BHT_-_SLAVE_-_*.BHR'],...
                [indir fs syr fs sjday fs 'CORR_-_MASTER_-_*.BHT_-_SLAVE_-_*.LHR'],...
                [indir fs syr fs sjday fs 'CORR_-_MASTER_-_*.LHT_-_SLAVE_-_*.BHR'],...
                [indir fs syr fs sjday fs 'CORR_-_MASTER_-_*.LHT_-_SLAVE_-_*.LHR']);
            tt=readseizmo(...
                [indir fs syr fs sjday fs 'CORR_-_MASTER_-_*.BHT_-_SLAVE_-_*.BHT'],...
                [indir fs syr fs sjday fs 'CORR_-_MASTER_-_*.BHT_-_SLAVE_-_*.LHT'],...
                [indir fs syr fs sjday fs 'CORR_-_MASTER_-_*.LHT_-_SLAVE_-_*.BHT'],...
                [indir fs syr fs sjday fs 'CORR_-_MASTER_-_*.LHT_-_SLAVE_-_*.LHT']);
        catch
            continue;
        end
        
        % get fk volume
        [rvol,tvol]=fkxchorzvolume(rr,rt,tr,tt,50,301,[1/100 1/4]);
        save(['fkvol.r.' syr '.' sjday '.mat'],'-struct','rvol');
        save(['fkvol.t.' syr '.' sjday '.mat'],'-struct','tvol');
    end
end

end
