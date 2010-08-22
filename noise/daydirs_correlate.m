function []=daydirs_correlate(indir,outdir,zpad,maxlag,o)
%DAYDIRS_CORRELATE    Correlates records in day directories
%
%    Usage:    daydirs_correlate(indir,outdir)
%              daydirs_correlate(indir,outdir,zpad,maxlag)
%              daydirs_correlate(indir,outdir,zpad,maxlag,overwrite)
%
%    Description: DAYDIRS_CORRELATE(INDIR,OUTDIR) correlates records under
%     INDIR, outputing the correlograms under OUTDIR.  INDIR should have a
%     day directory layout (ie INDIR/YEAR/JDAY/).  OUTDIR will have an
%     equivalent directory layout.  Records are zero padded from 0 to 90000
%     seconds prior to correlation and the correlograms are truncated to
%     +/- 4000 seconds lag time.
%
%     DAYDIRS_CORRELATE(INDIR,OUTDIR,ZPAD,MAXLAG) alters the default zero
%     padding limits and maximum lag time.  The defaults are [0 90000] for
%     ZPAD and 4000 for MAXLAG.
%
%     DAYDIRS_CORRELATE(INDIR,OUTDIR,ZPAD,MAXLAG,OVERWRITE) quietly
%     overwrites pre-existing records in OUTDIR when OVERWRITE is set to
%     TRUE.  By default OVERWRITE is FALSE.
%
%    Notes:
%
%    Header changes: See CORRELATE
%
%    Examples:
%
%    See also: DAYDIRS_MERGECUT_25HRS, DAYDIRS_RESAMPLE, DAYDIRS_NORMALIZE,
%              DAYDIRS_RINST, DAYDIRS_ROTCORR, DAYDIRS_STACKCORR,
%              DAYDIRS_MAKE

%     Version History:
%        June 20, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 20, 2010 at 12:55 GMT

% todo:

% check nargin
error(nargchk(2,5,nargin));

% defaults
if(nargin<3 || isempty(zpad)); zpad=[0 90000]; end
if(nargin<4 || isempty(maxlag)); maxlag=4000; end
if(nargin<5 || isempty(o)); o=false; end
if(numel(zpad)~=2 || ~isreal(zpad))
    error('seizmo:daydirs_correlate:badInput',...
        'ZPAD must be a 2 element real-valued vector!');
end
if(~isscalar(maxlag) || ~isreal(maxlag) || maxlag<=0)
    error('seizmo:daydirs_correlate:badInput',...
        'MAXLAG must be a positive scalar (in seconds)!');
end
if(~isscalar(o) || ~islogical(o))
    error('seizmo:daydirs_correlate:badInput',...
        'OVERWRITE flag must be a scalar logical!');
end

% check directories
if(~ischar(indir) || ~isvector(indir))
    error('seizmo:daydirs_correlate:fileNotString',...
        'INDIR must be a string!');
end
if(~exist(indir,'dir'))
    error('seizmo:daydirs_correlate:dirConflict',...
        ['Input Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(~ischar(outdir) || ~isvector(outdir))
    error('seizmo:daydirs_correlate:fileNotString',...
        'OUTDIR must be a string!');
end
if(exist(outdir,'file'))
    if(~exist(outdir,'dir'))
        error('seizmo:daydirs_correlate:dirConflict',...
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
    error('seizmo:daydirs_correlate:badLayout',...
        'Improper directory layout!');
end
jdays=cell(size(years));
for i=1:nyears
    % get day directories
    dirs=xdir([indir fs num2str(years(i))]);
    dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dir
    jdays{i}=str2double({dirs.name});
    if(any(isnan(jdays{i})))
        error('seizmo:daydirs_correlate:badLayout',...
            'Improper directory layout!');
    end
end

% verbosity (turn it off for the loop)
verbose=seizmoverbose(false);
if(verbose); disp('Correlating record(s)'); end

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
        
        % detail message
        if(verbose); disp(['PROCESSING DAY ' syr '.' sjday]); end
        
        % attempt correlation
        try
            % read verticals
            skip=false;
            try
                data=readseizmo([indir fs syr fs sjday fs '*.LHZ.*'],...
                    [indir fs syr fs sjday fs '*.BHZ.*']);
            catch
                skip=true;
            end
            if(numel(data)<2); skip=true; end
            
            if(~skip)
                % zero padding
                data=cut(data,'z',zpad(1),zpad(2),'fill',true);

                % correlating
                delta=getheader(data(1),'delta');
                data=correlate(data,'lags',(maxlag+4*delta).*[-1 1]);

                % interpolating
                data=interpolate(data,1,[],-maxlag,maxlag);

                % write verticals
                writeseizmo(data,'path',[outdir fs syr fs sjday fs]);
            end
            
            %%%
            
            % read horizontals
            skip=false;
            try
                data=readseizmo([indir fs syr fs sjday fs '*.LHE.*'],...
                    [indir fs syr fs sjday fs '*.LHN.*'],...
                    [indir fs syr fs sjday fs '*.BHE.*'],...
                    [indir fs syr fs sjday fs '*.BHN.*']);
            catch
                skip=true;
            end
            if(numel(data)<2); skip=true; end
            
            if(~skip)
                % zero padding
                data=cut(data,'z',zpad(1),zpad(2),'fill',true);

                % correlating
                delta=getheader(data(1),'delta');
                data=correlate(data,'lags',(maxlag+4*delta).*[-1 1]);

                % interpolating
                data=interpolate(data,1,[],-maxlag,maxlag);

                % write horizontals
                writeseizmo(data,'path',[outdir fs syr fs sjday fs]);
            end
        catch
            % close pool & fix verbosity
            tmp=lasterror;
            warning(tmp.message);
            %matlabpool close;
            seizmoverbose(verbose);
            
            % ???
            error(lasterror);
        end
    end
end

% parallel processing takedown
%matlabpool close;
seizmoverbose(verbose);

end
