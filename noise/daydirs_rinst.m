function []=daydirs_rinst(indir,outdir,varargin)
%DAYDIRS_RINST    Remove instrument response from records in day dirs
%
%    Usage:    daydirs_rinst(indir,outdir)
%              daydirs_rinst(...,'option',value,...)
%              daydirs_rinst(...,overwrite)
%
%    Description: DAYDIRS_RINST(INDIR,OUTDIR) removes the response from
%     records in day directory layout under directory INDIR.  A highpass
%     taper is applied with fullstop at 250s and fullpass at 150s to
%     stabilize the deconvolution.  Processed records are output in an
%     equivalent directory layout under OUTDIR.
%
%     DAYDIRS_RINST(...,'OPTION',VALUE,...) passes options to REMOVESACPZ.
%     See that function for available options.
%
%     DAYDIRS_RINST(...,OVERWRITE) quietly overwrites pre-existing records
%     in OUTDIR when OVERWRITE is set to TRUE.  By default OVERWRITE is
%     FALSE.  OVERWRITE must be the final argument.
%
%    Notes:
%
%    Header changes: DEP*, IDEP
%
%    Examples:
%
%    See also: DAYDIRS_MERGECUT_25HRS, DAYDIRS_RESAMPLE, DAYDIRS_NORMALIZE,
%              DAYDIRS_CORRELATE, DAYDIRS_ROTCORR, DAYDIRS_STACKCORR,
%              DAYDIRS_MAKE

%     Version History:
%        June 18, 2010 - initial version
%        Sep. 21, 2010 - commented out parallel processing lines
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 21, 2010 at 11:15 GMT

% todo:

% check nargin
error(nargchk(2,inf,nargin));

% extract overwrite option
if(mod(nargin,2))
    o=varargin{end};
    varargin(end)=[];
else
    o=false;
end

% check o
if(~isscalar(o) || ~islogical(o))
    error('seizmo:daydirs_rinst:badInput',...
        'OVERWRITE flag must be a scalar logical!');
end

% check directories
if(~ischar(indir) || ~isvector(indir))
    error('seizmo:daydirs_rinst:fileNotString',...
        'INDIR must be a string!');
end
if(~exist(indir,'dir'))
    error('seizmo:daydirs_rinst:dirConflict',...
        ['Input Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(~ischar(outdir) || ~isvector(outdir))
    error('seizmo:daydirs_rinst:fileNotString',...
        'OUTDIR must be a string!');
end
if(exist(outdir,'file'))
    if(~exist(outdir,'dir'))
        error('seizmo:daydirs_rinst:dirConflict',...
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
    error('seizmo:daydirs_rinst:badLayout',...
        'Improper directory layout!');
end
jdays=cell(size(years));
for i=1:nyears
    % get day directories
    dirs=xdir([indir fs num2str(years(i))]);
    dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dir
    jdays{i}=str2double({dirs.name});
    if(any(isnan(jdays{i})))
        error('seizmo:daydirs_rinst:badLayout',...
            'Improper directory layout!');
    end
end

% verbosity (turn it off for the loop)
verbose=seizmoverbose(false);
if(verbose); disp('Removing instrument response from record(s)'); end

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
        
        % attempt response removal
        try
            % read in data
            try
                data=readseizmo([indir fs syr fs sjday fs]);
            catch
                % empty day
                continue;
            end
            
            % remove response
            data=removesacpz(data,'f',[1/250 1/150],varargin{:});
            
            % skip if no rotated output
            if(~numel(data)); continue; end
            
            % write
            writeseizmo(data,'pathchange',{indir outdir});
        catch
            % close pool & fix verbosity
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
