function []=daydirs_ampspectra(indir,outdir,cmp,o)
%DAYDIRS_AMPSPECTRA    Stacks amplitude spectra of day directories
%
%    Usage:    daydirs_ampspectra(indir,outdir,cmp)
%              daydirs_ampspectra(indir,outdir,cmp,overwrite)
%
%    Description: DAYDIRS_AMPSPECTRA(INDIR,OUTDIR,CMP) calculates the
%     stacked amplitude spectra for records in a day directory layout under
%     directory INDIR.  The amplitude spectra are output in OUTDIR.  CMP is
%     the record component name CMP and is used as a file pattern.
%     Examples are 'LHZ', 'BHZ', 'LHN', etc.
%
%     DAYDIRS_AMPSPECTRA(INDIR,OUTDIR,OVERWRITE) quietly overwrites
%     pre-existing records in OUTDIR when OVERWRITE is set to TRUE.  By
%     default OVERWRITE is FALSE.
%
%    Notes:
%     - Output is general xy records with the filename YYYY.DDD.CMP.SAC
%
%    Header changes: TOTAL
%
%    Examples:
%
%    See also: DAYDIRS_MERGECUT_25HRS, DAYDIRS_RINST, DAYDIRS_NORMALIZE,
%              DAYDIRS_CORRELATE, DAYDIRS_ROTCORR, DAYDIRS_STACKCORR,
%              DAYDIRS_MAKE, DAYDIRS_RESAMPLE

%     Version History:
%        June 30, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 30, 2010 at 12:55 GMT

% todo:

% check nargin
error(nargchk(3,4,nargin));

% defaults
if(nargin<4 || isempty(o)); o=false; end
if(~ischar(cmp) || ~isequal(size(cmp),[1 3]))
    error('seizmo:daydirs_ampspectra:badInput',...
        'CMP1 must be a 3 character component code!');
end
if(~isscalar(o) || ~islogical(o))
    error('seizmo:daydirs_ampspectra:badInput',...
        'OVERWRITE flag must be a scalar logical!');
end

% check directories
if(~ischar(indir) || ~isvector(indir))
    error('seizmo:daydirs_ampspectra:fileNotString',...
        'INDIR must be a string!');
end
if(~exist(indir,'dir'))
    error('seizmo:daydirs_ampspectra:dirConflict',...
        ['Input Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(~ischar(outdir) || ~isvector(outdir))
    error('seizmo:daydirs_resample:fileNotString',...
        'OUTDIR must be a string!');
end
if(exist(outdir,'file'))
    if(~exist(outdir,'dir'))
        error('seizmo:daydirs_ampspectra:dirConflict',...
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

% get directory separator
fs=filesep;

% parallel processing setup (8 instances)
matlabpool(8);

% get year directories and day directories
dirs=xdir([indir fs]);
dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dirs
years=str2double({dirs.name});
nyears=numel(years);
if(any(isnan(years)))
    error('seizmo:daydirs_ampspectra:badLayout',...
        'Improper directory layout!');
end
jdays=cell(size(years));
for i=1:nyears
    % get day directories
    dirs=xdir([indir fs num2str(years(i))]);
    dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dir
    jdays{i}=str2double({dirs.name});
    if(any(isnan(jdays{i})))
        error('seizmo:daydirs_ampspectra:badLayout',...
            'Improper directory layout!');
    end
end

% verbosity (turn it off for the loop)
verbose=seizmoverbose(false);
if(verbose); disp('Getting amplitude spectra of array record(s)'); end

% loop over years
for i=1:nyears
    % working year
    yr=years(i);
    syr=num2str(yr);
    
    % loop over days
    parfor j=1:numel(jdays{i})
        % working julian day
        jday=jdays{i}(j);
        sjday=num2str(jday,'%03d');
        
        % detail message
        if(verbose); disp(['PROCESSING DAY ' syr '.' sjday]); end
        
        % attempt resample
        try
            % read in data
            try
                data=readseizmo([indir fs syr fs sjday fs '*' cmp '*']);
            catch
                % empty day
                continue;
            end
            
            % get stacked amp spectra
            nrecs=numel(data);
            data=divide(stack(keepam(dft(data))),nrecs);
            
            % clear out header info that will confuse
            data=changeheader(data,'st',nan,'ev',nan,'cmp',nan,...
                'kname',nan,'scale',nrecs);
            
            % write
            writeseizmo(data,'path',outdir,...
                'name',[syr '.' sjday '.' cmp '.SAC']);
        catch
            % close pool & fix verbosity
            tmp=lasterror;
            warning(tmp.message);
            matlabpool close;
            seizmoverbose(verbose);
            
            % ???
            error(lasterror);
        end
    end
end

% parallel processing takedown
matlabpool close;
seizmoverbose(verbose);

end
