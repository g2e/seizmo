function []=daydirs_rotcorr(indir,outdir,o)
%DAYDIRS_ROTCORR    Rotates correlograms in day directories
%
%    Usage:    daydirs_rotcorr(indir,outdir)
%              daydirs_rotcorr(indir,outdir,overwrite)
%
%    Description: DAYDIRS_ROTCORR(INDIR,OUTDIR) rotates EE, EN, NE, NN
%     correlograms under INDIR to RR, RT, TR, & TT under OUTDIR.  R & T
%     correspond to the radial & transverse directions with respect to the
%     interstation azimuths (radial points from the master station to the
%     slave station and transverse leads the radial component by 90
%     degrees).  See Lin et al 2008 for more details.
%
%     DAYDIRS_ROTCORR(INDIR,OUTDIR,OVERWRITE)  quietly overwrites
%     pre-existing records in OUTDIR when OVERWRITE is set to TRUE.  By
%     default OVERWRITE is FALSE.
%
%    Notes:
%
%    Header changes: USER0, USER1, KCMPNM, KT3
%                    Master & Slave field info may be switched
%
%    Examples:
%
%    See also: DAYDIRS_MERGECUT_25HRS, DAYDIRS_RESAMPLE, DAYDIRS_NORMALIZE,
%              DAYDIRS_CORRELATE, DAYDIRS_STACKCORR, DAYDIRS_RINST,
%              DAYDIRS_MAKE

%     Version History:
%        June 20, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 20, 2010 at 12:55 GMT

% todo:

% check nargin
error(nargchk(2,3,nargin));

% defaults
if(nargin<3 || isempty(o)); o=false; end
if(~isscalar(o) || ~islogical(o))
    error('seizmo:daydirs_rotcorr:badInput',...
        'OVERWRITE flag must be a scalar logical!');
end

% check directories
if(~ischar(indir) || ~isvector(indir))
    error('seizmo:daydirs_rotcorr:fileNotString',...
        'INDIR must be a string!');
end
if(~exist(indir,'dir'))
    error('seizmo:daydirs_rotcorr:dirConflict',...
        ['Input Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(~ischar(outdir) || ~isvector(outdir))
    error('seizmo:daydirs_rotcorr:fileNotString',...
        'OUTDIR must be a string!');
end
if(exist(outdir,'file'))
    if(~exist(outdir,'dir'))
        error('seizmo:daydirs_rotcorr:dirConflict',...
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
matlabpool(8);

% get year directories and day directories
dirs=xdir([indir fs]);
dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dirs
years=str2double({dirs.name});
nyears=numel(years);
if(any(isnan(years)))
    error('seizmo:daydirs_rotcorr:badLayout',...
        'Improper directory layout!');
end
jdays=cell(size(years));
for i=1:nyears
    % get day directories
    dirs=xdir([indir fs num2str(years(i))]);
    dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dir
    jdays{i}=str2double({dirs.name});
    if(any(isnan(jdays{i})))
        error('seizmo:daydirs_rotcorr:badLayout',...
            'Improper directory layout!');
    end
end

% verbosity (turn it off for the loop)
verbose=seizmoverbose(false);
if(verbose); disp('Rotating Correlograms'); end

% loop over years
for i=1:nyears
    % working year
    yr=years(i);
    syr=num2str(yr);
    
    % loop over days
    parfor j=1:numel(jdays)
        % working julian day
        jday=jdays{i}(j);
        sjday=num2str(jday,'%03d');
        
        % only for restarting
        %if(yr<2006 || (yr==2006 && jday<256))
        %    continue;
        %end
        
        % detail message
        if(verbose); disp(['PROCESSING DAY ' syr '.' sjday]); end
        
        % attempt rotation
        try
            % read
            try
                data=readseizmo([indir fs syr fs sjday fs]);
            catch
                % empty day
                continue;
            end
            
            % rotate horizontal correlations
            data=rotate_correlations(data);
            
            % skip if no rotated output
            if(~numel(data)); continue; end
            
            % write
            writeseizmo(data,'pathchange',{indir outdir});
        catch
            % close pool & fix verbosity
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
