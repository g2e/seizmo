function [geofs]=fs_phase2latlon(fs,phase,llstep)
%FS_PHASE2LATLON    Converts a fs spectra to geofs using a given phase
%
%    Usage:    geofs=fs_phase2latlon(fs,phase)
%              geofs=fs_phase2latlon(fs,phase,llstep)
%
%    Description:
%     GEOFS=FS_PHASE2LATLON(FS,PHASE) converts the frequency-slowness power
%     spectra FS into a geographic frequency-slowness power spectra GEOFS
%     using the body-wave phases in PHASE to convert slowness & azimuth to
%     position (latitude & longitude).  The GEOFS structure has the
%     .horzslow field set to nan because the spectra has variable slowness
%     across the map.
%
%     GEOFS=FS_PHASE2LATLON(FS,PHASE,LLSTEP) uses the latitude/longitude
%     step size given by LLSTEP.  LLSTEP is in degrees and by default is 1.
%
%    Notes:
%
%    Examples:
%     % Convert using P & PKPbc or PP:
%     geofs=fs_phase2latlon(fs,{'P' 'PKPbc'});
%     geofs=fs_phase2latlon(fs,'PP');
%
%     % Plotting:
%     plotgeofs(geofs);
%
%    See also: FS_CART2POL, PLOTGEOFS, FSSPECTRA

%     Version History:
%        Apr. 24, 2012 - initial version
%        May  18, 2012 - hack for arf plotting
%        Aug. 10, 2012 - remove debugging output
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 10, 2012 at 16:05 GMT

% todo:

% check nargin
error(nargchk(2,3,nargin));

% default llstep
if(nargin<3 || isempty(llstep)); llstep=1; end

% check phase
validphase={'P','PP','PKPab','PKPbc','PKPdf','PKiKP','PKIKP'};
if(ischar(phase)); phase=cellstr(phase); end
if(~iscellstr(phase) || any(~ismember(phase,validphase)))
    error('seizmo:fs_phase2latlon:badPhase',...
        ['PHASE must be one of the following strings:\n' ...
        sprintf('%s ',validphase{:})]);
end
nphase=numel(phase);

% check fk struct
if(isfkarfstruct(fs))
    fs.stel=fs.stla*0;
    fs.stdp=fs.stla*0;
    fs.freq=fs.f;
    fs.butc=zeros(1,5);
    fs.eutc=zeros(1,5);
    fs.volume=false;
    fs.npts=1;
    fs.delta=1;
end
error(chkfkstruct(fs));
nfs=numel(fs);

% check lat/lon step size
if(~isscalar(llstep) || ~isnumeric(llstep) || ~isreal(llstep) || llstep<=0)
    error('seizmo:fs_phase2latlon:badStep',...
        'LLSTEP must be a real-valued positive scalar!');
end

% lat/lon grid
[lat,lon]=meshgrid(-90+llstep/2:llstep:90,-180+llstep/2:llstep:180);
nll=numel(lat);

% create output structure
[geofs(1:nfs).nsta]=deal(fs.nsta);
[geofs(1:nfs).stla]=deal(fs.stla);
[geofs(1:nfs).stlo]=deal(fs.stlo);
[geofs(1:nfs).stel]=deal(fs.stel);
[geofs(1:nfs).stdp]=deal(fs.stdp);
[geofs(1:nfs).butc]=deal(fs.butc);
[geofs(1:nfs).eutc]=deal(fs.eutc);
[geofs(1:nfs).delta]=deal(fs.delta);
[geofs(1:nfs).npts]=deal(fs.npts);
[geofs(1:nfs).volume]=deal(false(1,2)); % altered in the loop below
[geofs(1:nfs).latlon]=deal([lat(:) lon(:)]);
[geofs(1:nfs).horzslow]=deal(nan); % variable slowness so nan is better
[geofs(1:nfs).npairs]=deal(fs.npairs);
[geofs(1:nfs).method]=deal(fs.method);
[geofs(1:nfs).center]=deal(fs.center);
[geofs(1:nfs).weights]=deal(fs.weights);
[geofs(1:nfs).freq]=deal(fs.freq);
[geofs(1:nfs).beam]=deal([]); % created in loop below
[geofs(1:nfs).normdb]=deal(fs.normdb);

% verbosity
verbose=seizmoverbose;    

% convert cartesian grid to polar
if(any(~[fs.polar]))
    % detail message
    if(verbose); disp('Converting fs spectra from cartesian to polar'); end
    arecart=~[fs.polar];
    fs(arecart)=fkcart2pol(fs(arecart));
end

% loop over each element in map (separate times or something like that)
for a=1:nfs
    % detail message
    if(verbose && nfs>1)
        fprintf('Working on fs spectra frame %d of %d\n',a,nfs);
    end
    
    % get deg/baz for all lat/lon positions
    [deg,baz]=sphericalinv(fs(a).center(1),fs(a).center(2),lat(:),lon(:));
    
    % initialize output
    nfreq=size(fs(a).beam,3);
    geofs(a).volume=[false fs(a).volume];
    geofs(a).beam=nan([nll 1 nfreq]);
    
    % loop over phases
    for b=nphase:-1:1
        % detail message
        if(verbose); fprintf('Getting %s geofs spectra\n',phase{b}); end
        
        % get slowness given distance
        slow=deg2slowness(phase{b},deg);
        notnans=find(~isnan(slow));
        
        % detail message
        if(verbose); print_time_left(0,nfreq); end
        
        % loop over frequencies
        for c=1:nfreq
            % interpolate dB given baz/slowness
            % x == cols == baz
            % y == rows == slow
            if(fs(a).x(end)==360)
                geofs(a).beam(notnans+nll*(c-1))=interp2(...
                    fs(a).x,fs(a).y,fs(a).beam(:,:,c),...
                    baz(notnans),slow(notnans));
            else
                % need to help with wraparound
                geofs(a).beam(notnans+nll*(c-1))=interp2(...
                    [fs(a).x 360],fs(a).y,fs(a).beam(:,[1:end end],c),...
                    baz(notnans),slow(notnans));
            end
            
            % detail message
            if(verbose); print_time_left(c,nfreq); end
        end
    end
end

end
