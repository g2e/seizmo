function []=srtm30plus2mat(varargin)
%SRTM30PLUS2MAT    Create mat-files from SRTM30PLUS SRTM files
%
%    Usage:    srtm30plus2mat(filelist)
%
%    Description:
%     SRTM30PLUS2MAT(FILELIST1,...,FILELISTN) converts SRTM30plus files
%     into a Matlab .mat file for simpler loading.  Accepts char arrays of
%     filenames (one filename per row) and/or cell arrays of filenames (one
%     filename per cell).  Wildcards are valid.  Passing an empty list of
%     files or no argument at all will bring up a graphical file selection
%     menu.  The .mat file is saved as 'srtm30plus.mat' in the current
%     working directory.  SRTM grids are saved as variables according to
%     the string preceeding the first '.' in their name.  So
%     'w020n40.Bathymetry.srtm' would be saved under the variable name
%     'w020n40'.  The string is always lowercased.
%
%    Notes:
%     - Watch out for variable name conflicts (what are you doing?)!
%
%    Examples:
%     % Convert SRTM files in the current directory to a mat-file:
%     srtm30plus2mat('*.srtm');
%
%     % Now load just the grids of Africa:
%     load srtm30plus w020n40 e020n40 w020s10 e020s10
%
%     % Plot Africa:
%     imagesc([w020n40 e020n40; w020s10 e020s10])
%
%    See also: SRTM30PLUS_TILE, TOPO_POINTS, TOPO_REGION

%     Version History:
%        Feb. 14, 2010 - initial version
%        Apr.  2, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  2, 2012 at 16:50 GMT

% todo:

% set filterspec appropriately
global SEIZMO
SEIZMO.ONEFILELIST.FILTERSPEC=...
    {'*.srtm;*.SRTM' 'SRTM Files (*.srtm,*.SRTM)'};

% compile file lists
varargin=onefilelist(varargin{:});
nfiles=numel(varargin);

% error if no files
if(nfiles<1)
    error('seizmo:srtm2mat:noFilesToConvert',...
        'No files to convert!');
end

% verbosity
verbose=seizmoverbose;
if(verbose)
    disp('Converting SRTM File(s) to MAT');
    print_time_left(0,nfiles);
end

% loop for each file
for i=1:nfiles
    % get name
    fullname=[varargin(i).path varargin(i).name];
    name=getwords(varargin(i).name,'.');
    name=lower(name{1});
    
    % check number of bytes
    if(varargin(i).bytes==57600000)
        ncols=4800;
    elseif(varargin(i).bytes==51840000)
        ncols=7200;
    else
        warning('seizmo:srtm30plus2mat:badBytes',...
            'File not properly formatted SRTM30plus File: %s !',fullname);
        % detail message
        if(verbose); print_time_left(i,nfiles,true); end
        continue;
    end
    
    % open file as big-endian
    fid=fopen(fullname,'r','ieee-be');
    
    % fid check
    if(fid<0)
        % bad permissions?
        warning('seizmo:readheader:badFID',...
            'File not openable, %s !',fullname);
        % detail message
        if(verbose); print_time_left(i,nfiles,true); end
        continue;
    end
    
    % read file
    srtm30plus.(name)=fread(fid,inf,'*int16');
    
    % close file
    fclose(fid);
    
    % reshape based on size
    srtm30plus.(name)=reshape(srtm30plus.(name),ncols,[])';
    
    % detail message
    if(verbose); print_time_left(i,nfiles); end
end

% save to srtm30plus.mat in current directory
if(verbose); disp('Saving To ./srtm30plus.mat'); end
save ./srtm30plus.mat -struct srtm30plus

end
