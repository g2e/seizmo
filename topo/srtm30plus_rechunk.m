function []=srtm30plus_rechunk(version)
%SRTM30PLUS_RECHUNK    Rechunks SRTM30plus into 10x10 blocks
%
%    Usage:    srtm30plus_rechunk(version)
%
%    Description: SRTM30PLUS_RECHUNK(VERSION) reads the SRTM mat file
%     'srtm30plus.mat' produced by SRTM30PLUS2MAT in the current working
%     directory and reformats it to be in 10x10 degree chunks.  This makes
%     loads faster, lighter, and simpler.  Output mat-file is
%     'srtm30plus10' in the current working directory.
%
%    Notes:
%
%    Examples:
%     Convert and rechunk (setting version to 6):
%      srtm30plus2mat('*.srtm');
%      srtm30plus_rechunk(6);
%
%    See also: SRTM30PLUS2MAT

%     Version History:
%        Feb. 14, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 14, 2010 at 18:50 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% verbosity
verbose=seizmoverbose;

% new chunk size (must divide 90 evenly)
w=10;
nnchunks=360*180/w^2;

% nodes per degree
n=120;

% make chunknames
chunkname={'w180n90' 'w140n90' 'w100n90' 'w060n90' 'w020n90' 'e020n90' ...
    'e060n90' 'e100n90' 'e140n90' 'w180n40' 'w140n40' 'w100n40' ...
    'w060n40' 'w020n40' 'e020n40' 'e060n40' 'e100n40' 'e140n40' ...
    'w180s10' 'w140s10' 'w100s10' 'w060s10' 'w020s10' 'e020s10' ...
    'e060s10' 'e100s10' 'e140s10' 'w180s60' 'w120s60' 'w060s60' ...
    'w000s60' 'e060s60' 'e120s60'};
nchunks=numel(chunkname);

% offsets (not sure how to do this automatically)
lonoff=[0:4:32 0:4:32 0:4:32 0:6:30];
latoff=[zeros(1,9) 5*ones(1,9) 10*ones(1,9) 15*ones(1,6)];

% set other fields
srtm30plus.version=version;
srtm30plus.latname=[strcat('n',cellstr(num2str((90:-w:w)','%02d'))); ...
    strcat('s',cellstr(num2str((0:w:90-w)','%02d')))];
srtm30plus.lonname=[strcat('w',cellstr(num2str((180:-w:w)','%03d'))); ...
    strcat('e',cellstr(num2str((0:w:180-w)','%03d')))];

% loop over new chunks
if(verbose)
    disp('Rechunking');
    cnt=0;
    print_time_left(cnt,nnchunks);
end
for c=1:nchunks
    % load chunk
    s=load('./srtm30plus.mat',chunkname{c});
    s=s.(char(fieldnames(s)));
    nx=size(s,2)/(n*w);
    ny=size(s,1)/(n*w);
    
    % make new chunks
    for i=1:nx
        for j=1:ny
            srtm30plus.(strcat(srtm30plus.lonname{i+lonoff(c)}, ...
                srtm30plus.latname{j+latoff(c)}))=...
                s((j-1)*1200+(1:1200),(i-1)*1200+(1:1200));
            if(verbose); cnt=cnt+1; print_time_left(cnt,nnchunks); end
        end
    end
end

% save new srtm30plus
if(verbose); disp('Saving New Chunks to mat-file'); end
save ./srtm30plus10.mat -struct srtm30plus

end
