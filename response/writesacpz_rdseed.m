function []=writesacpz_rdseed(pz,filename,oflag)
%WRITESACPZ_RDSEED    Writes SAC polezero files equivalent to RDSEED
%
%    Usage:    writesacpz_rdseed(pz)
%              writesacpz_rdseed(pz,filename)
%              writesacpz_rdseed(...,overwrite)
%
%    Description:
%     WRITESACPZ_RDSEED(PZ) writes the SAC polezero information contained
%     in the structure PZ to file(s).  PZ must be formatted as from
%     READSACPZ_RDSEED.  This attempts to replicate the RDSEED SAC polezero
%     output format (in particular the comment block) and uses the .path &
%     .name fields to determine which file(s) are written.  Note that this
%     function will ask to overwrite EACH file that already exists (see the
%     last usage form to overwrite without user input).
%
%     WRITESACPZ_RDSEED(PZ,FILENAME) write to output file FILENAME instead
%     of writing to the file(s) as determined from the .path & .name fields
%     in PZ.  This is useful for making a single SAC polezero file.
%
%     WRITESACPZ_RDSEED(...,OVERWRITE) quietly overwrites pre-existing SAC
%     PoleZero files without confirmation when OVERWRITE is set to TRUE.
%     By default OVERWRITE is FALSE.
%
%    Notes:
%     - WRITESACPZ_RDSEED replaced DB2SACPZ.
%
%    Examples:
%     % Grab ANMO & CMB polezero info and write to a single file:
%     url=['http://service.iris.edu/irisws/sacpz/1/' ...
%          'query?net=IU&loc=*&cha=*&sta=ANMO'];
%     pz=readsacpz_rdseed(urlread(url),true);
%     url=['http://service.iris.edu/irisws/sacpz/1/' ...
%          'query?net=BK&loc=*&cha=*&sta=CMB'];
%     pz=[pz; readsacpz_rdseed(urlread(url),true)];
%     writesacpz_rdseed(pz,'SAC_PZs_IU.ANMO_BK.CMB');
%
%    See also: READSACPZ_RDSEED, WRITESACPZ, READSACPZ, REMOVESACPZ,
%              APPLYSACPZ, MAKESACPZDB, GENSACPZNAME, PARSE_SACPZ_FILENAME,
%              GETSACPZ, ISSCAPZ_RDSEED, FIX_OLD_SACPZ, SSIDX, SSCAT

%     Version History:
%        Feb. 25, 2014 - initial version
%        Mar.  6, 2014 - update See also section
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  6, 2014 at 15:05 GMT

% todo:

% check number of inputs
error(nargchk(1,3,nargin));

% check pz struct
if(~issacpz_rdseed(pz))
    error('seizmo:writesacpz_rdseed:badInput',...
        'PZ must be a struct as from READSACPZ_RDSEED!');
end

% check for overwrite flag
if(nargin<3 || isempty(oflag)); oflag=false; end
if(nargin==2 && islogical(filename) && isscalar(filename))
    oflag=filename;
    filename=[];
end
if(~isscalar(oflag) || ~islogical(oflag))
    error('seizmo:writesacpz_rdseed:badInput',...
        'OVERWRITE flag must be a scalar logical!');
end

% check filename input
if(nargin==1); filename=[]; end
if(nargin>1 && ~isempty(filename) && ~isstring(filename))
    error('seizmo:writesacpz_rdseed:badInput',...
        'FILENAME must be a string!');
end

% override filenames if desired
if(~isempty(filename))
    [path,name,ext]=fileparts(filename);
    if(isempty(path)); path='.'; end
    name=[name ext];
    [pz.path{:}]=deal([path filesep]);
    [pz.name{:}]=deal(name);
end

% get filenames
filenames=strcat(pz.path,pz.name);

% group by filename
[filename,first,grpid]=unique(filenames);
nfiles=numel(filename);

% loop over each group and write
for i=1:nfiles
    % check for conflict
    if(exist(filename{i},'file'))
        if(exist(filename{i},'dir'))
            error('seizmo:writesacpz_rdseed:dirConflict',...
                'SAC PoleZero File: %s\nIs A Directory!',filename{i});
        end
        if(~oflag)
            fprintf('SAC PoleZero File: %s\nFile Exists!\n',filename{i});
            reply=input('Overwrite? Y/N [N]: ','s');
            if(isempty(reply) || ~strncmpi(reply,'y',1))
                disp('Not overwriting!');
                continue;
            end
            disp('Overwriting!');
        end
    end
    
    % open file for writing as ascii
    fid=fopen(filename{i},'wt');
    
    % check if file is openable
    if(fid<0)
        error('seizmo:writesacpz_rdseed:cannotOpenFile',...
            'SAC PoleZero File: %s\nNot Openable!',file);
    end
    
    % who is in this group
    in=find(i==grpid)';
    
    % loop over polezeros writing the info to file
    for j=in
        % write comment block
        fprintf(fid,'* **********************************\n');
        fprintf(fid,'* NETWORK   (KNETWK): %s\n',upper(pz.knetwk{j}));
        fprintf(fid,'* STATION    (KSTNM): %s\n',upper(pz.kstnm{j}));
        fprintf(fid,'* LOCATION   (KHOLE): %s\n',upper(pz.khole{j}));
        fprintf(fid,'* CHANNEL   (KCMPNM): %s\n',upper(pz.kcmpnm{j}));
        fprintf(fid,'* CREATED           : %s\n',datestr(...
            [doy2cal(pz.created(j,1:2)) pz.created(j,3:5)],...
            'yyyy-mm-ddTHH:MM:SS'));
        fprintf(fid,'* START             : %s\n',datestr(...
            [doy2cal(pz.b(j,1:2)) pz.b(j,3:5)],'yyyy-mm-ddTHH:MM:SS'));
        % DATESTR fails when sum(caltime)-2000>500
        % - this happens around when the year exceeds 2300 which is
        %   common for the end time of many polezero timespans
        % - we just attempt to write the string using fprintf instead
        fprintf(fid,...
            '* END               : %04d-%02d-%02dT%02d:%02d:%02d\n',...
            [doy2cal(pz.e(j,1:2)) pz.e(j,3:5)]);
        fprintf(fid,'* DESCRIPTION       : %s\n',pz.description{j});
        fprintf(fid,'* LATITUDE          : %.6f\n',pz.stla(j));
        fprintf(fid,'* LONGITUDE         : %.6f\n',pz.stlo(j));
        fprintf(fid,'* ELEVATION         : %.6f\n',pz.stel(j));
        fprintf(fid,'* DEPTH             : %.6f\n',pz.stdp(j));
        fprintf(fid,'* DIP               : %.6f\n',pz.cmpinc(j));
        fprintf(fid,'* AZIMUTH           : %.6f\n',pz.cmpaz(j));
        fprintf(fid,'* SAMPLE RATE       : %.6f\n',pz.sr(j));
        fprintf(fid,'* INPUT UNIT        : %s\n',pz.input{j});
        fprintf(fid,'* OUTPUT UNIT       : %s\n',pz.output{j});
        fprintf(fid,'* INSTTYPE          : %s\n',pz.insttype{j});
        if(~isnan(pz.instgain(j)))
            fprintf(fid,'* INSTGAIN          : %e %s\n',...
                pz.instgain(j),pz.instgainunits{j});
        else
            fprintf(fid,'* INSTGAIN          : \n');
        end
        fprintf(fid,'* COMMENT           : %s\n',pz.comment{j});
        if(~isnan(pz.sensitivity(j)))
            fprintf(fid,'* SENSITIVITY       : %e %s\n',...
                pz.sensitivity(j),pz.sensitivityunits{j});
        else
            fprintf(fid,'* SENSITIVITY       : \n');
        end
        fprintf(fid,'* A0                : %e\n',pz.a0(j));
        fprintf(fid,'* **********************************\n');
        
        % extract poles & zeros
        z=pz.z{j};
        p=pz.p{j};
        
        % get total number of poles/zeros
        nz=numel(z);
        np=numel(p);
        
        % remove poles/zeros at origin
        z(abs(z)==0)=[];
        p(abs(p)==0)=[];
        
        % get new number of poles/zeros
        nnz=numel(z);
        nnp=numel(p);
        
        % write ZPK info
        fprintf(fid,'ZEROS %d\n',nz);
        if(nnz)
            fprintf(fid,'\t%+e\t%+e\n',[real(z(:).'); imag(z(:).')]);
        end
        fprintf(fid,'POLES %d\n',np);
        if(nnp)
            fprintf(fid,'\t%+e\t%+e\n',[real(p(:).'); imag(p(:).')]);
        end
        fprintf(fid,'CONSTANT %e\n\n\n',pz.k(j));
    end
    
    % close file
    fclose(fid);
end

end
