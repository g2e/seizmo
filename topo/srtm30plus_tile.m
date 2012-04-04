function []=srtm30plus_tile(version)
%SRTM30PLUS_TILE    Tiles SRTM30plus into 10x10 tiles
%
%    Usage:    srtm30plus_tile(version)
%
%    Description:
%     SRTM30PLUS_TILE(VERSION) reads the SRTM mat file 'srtm30plus.mat'
%     produced by SRTM30PLUS2MAT in the current working directory and
%     reformats it to be in 10x10 degree tiles.  This makes loads faster,
%     lighter, and simpler.  Output mat-file is 'srtm30plus10' in the
%     current working directory.
%
%    Notes:
%
%    Examples:
%     % Convert and tile (setting version to 6):
%     srtm30plus2mat('*.srtm');
%     srtm30plus_tile(6);
%
%    See also: SRTM30PLUS2MAT

%     Version History:
%        Feb. 14, 2010 - initial version
%        Feb. 16, 2010 - significant memory use reduction
%        Feb. 11, 2011 - mass nargchk fix, fix etopo1_bed bug
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
n=120;

% make tilenames
tilename={'w180n90' 'w140n90' 'w100n90' 'w060n90' 'w020n90' 'e020n90' ...
    'e060n90' 'e100n90' 'e140n90' 'w180n40' 'w140n40' 'w100n40' ...
    'w060n40' 'w020n40' 'e020n40' 'e060n40' 'e100n40' 'e140n40' ...
    'w180s10' 'w140s10' 'w100s10' 'w060s10' 'w020s10' 'e020s10' ...
    'e060s10' 'e100s10' 'e140s10' 'w180s60' 'w120s60' 'w060s60' ...
    'w000s60' 'e060s60' 'e120s60'};
ntiles=numel(tilename);

% offsets (not sure how to do this automatically)
lonoff=[0:4:32 0:4:32 0:4:32 0:6:30];
latoff=[zeros(1,9) 5*ones(1,9) 10*ones(1,9) 15*ones(1,6)];

% set other fields
srtm30plus.version=version;
srtm30plus.registration='pixel';
srtm30plus.pixelsperdegree=n;
srtm30plus.latname=[strcat('n',cellstr(num2str((90:-w:w)','%02d'))); ...
    strcat('s',cellstr(num2str((0:w:90-w)','%02d')))];
srtm30plus.lonname=[strcat('w',cellstr(num2str((180:-w:w)','%03d'))); ...
    strcat('e',cellstr(num2str((0:w:180-w)','%03d')))];
save ./srtm30plus10.mat -struct srtm30plus

% detail message
if(verbose)
    disp('Tiling SRTM30+');
    cnt=0;
    print_time_left(cnt,nntiles);
end

% loop over tiles
for c=1:ntiles
    % load tile
    s=load('./srtm30plus.mat',tilename{c});
    s=s.(char(fieldnames(s)));
    nx=size(s,2)/(n*w);
    ny=size(s,1)/(n*w);
    
    % make new tiles
    for i=1:nx
        for j=1:ny
            % extract tile
            tmp.(strcat(srtm30plus.lonname{i+lonoff(c)}, ...
                srtm30plus.latname{j+latoff(c)}))=...
                s((j-1)*n*w+(1:n*w),(i-1)*n*w+(1:n*w));
            
            % save and clear to keep memory footprint low
            save ./srtm30plus10.mat -struct tmp -append
            clear tmp
            
            % detail message
            if(verbose); cnt=cnt+1; print_time_left(cnt,nntiles); end
        end
    end
end

end
