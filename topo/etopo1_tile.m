function []=etopo1_tile(version)
%ETOPO1_TILE    Tiles ETOPO1 into 10x10 tiles
%
%    Usage:    etopo1_tile(version)
%
%    Description:
%     ETOPO1_TILE(VERSION) reads the ETOPO1 mat files 'etopo1_*.mat' in the
%     current working directory and reformats them to be in 10x10 degree
%     tiles.  This makes loads faster and lighter.  The output mat-files
%     are 'etopo1_*10.mat' in the current working directory.
%
%    Notes:
%     - Due to grid registration the tiles will have overlap.  Use
%       TOPO_REGION to extract the region properly.
%
%    Examples:
%     % Tile (setting version to 1):
%     etopo1_tile(1);
%
%    See also: SRTM30PLUS_TILE

%     Version History:
%        Feb. 16, 2010 - initial version
%        Feb. 25, 2010 - minor doc update
%        Feb. 11, 2011 - mass nargchk fix
%        Apr.  2, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  2, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% verbosity
verbose=seizmoverbose;

% new tile size (must divide 90 evenly)
w=10;
nntiles=360*180/w^2;

% pixels per degree
n=60;

% set other fields
etopo1_bed.version=version;
etopo1_bed.registration='grid';
etopo1_bed.pixelsperdegree=n;
etopo1_bed.latname=[strcat('n',cellstr(num2str((90:-w:w)','%02d'))); ...
    strcat('s',cellstr(num2str((0:w:90-w)','%02d')))];
etopo1_bed.lonname=[strcat('w',cellstr(num2str((180:-w:w)','%03d'))); ...
    strcat('e',cellstr(num2str((0:w:180-w)','%03d')))];
save ./etopo1_bed10.mat -struct etopo1_bed
save ./etopo1_ice10.mat -struct etopo1_bed
save ./etopo1_thk10.mat -struct etopo1_bed

% loop over tiles
tile={'bed' 'ice' 'thk'};
for c=1:3
    % detail message
    if(verbose)
        disp(['Loading ETOPO1_' upper(tile{c})]);
    end
    
    % load tile
    s=load(['./etopo1_' tile{c} '.mat']);
    s=s.(['etopo1_' tile{c}]);
    nx=(size(s,2)-1)/(n*w);
    ny=(size(s,1)-1)/(n*w);
    
    % detail message
    if(verbose)
        disp(['Tiling ETOPO1_' upper(tile{c})]);
        cnt=0;
        print_time_left(cnt,nntiles);
    end
    
    % make new tiles
    for i=1:nx
        for j=1:ny
            % extract tile (note that there will be tile overlap)
            tmp.(strcat(etopo1_bed.lonname{i}, ...
                etopo1_bed.latname{j}))=...
                s((j-1)*n*w+(1:n*w+1),(i-1)*n*w+(1:n*w+1));
            
            % save and clear to keep memory footprint low
            save(['./etopo1_' tile{c} '10.mat'],'-struct','tmp','-append');
            clear tmp
            
            % detail message
            if(verbose); cnt=cnt+1; print_time_left(cnt,nntiles); end
        end
    end
end

end
