function [topo]=topo_points(lat,lon,toponame)
%TOPO_POINTS    Returns topography data at the points indicated
%
%    Usage:    topo=topo_points(lat,lon)
%              topo=topo_points(lat,lon,toponame)
%
%    Description:
%     TOPO=TOPO_POINTS(LAT,LON) grabs SRTM30+ topography in meters for the
%     positions given by LAT & LON.  LAT & LON must be in degrees and must
%     be scalar or equal-sized arrays.
%
%     TOPO=TOPO_POINTS(LAT,LON,TOPONAME) changes the dataset from which the
%     topography is pulled.  TOPONAME must be a string and currently
%     valid topography models are:
%      'srtm30+'    - SRTM30+ topography (30 arc-second resolution)
%                     (http://topex.ucsd.edu/WWW_html/srtm30_plus.html)
%      'etopo1_bed' - ETOPO1 bedrock topography (1 arc-minute resolution)
%      'etopo1_ice' - ETOPO1 ice topography (1 arc-minute resolution)
%                     (http://www.ngdc.noaa.gov/mgg/global/global.html)
%      'crust1.0'   - Crustal model at 1x1deg blocks (w/ topography)
%      'crust2.0'   - Crustal model at 2x2deg blocks (w/ topography)
%      'crust5.1'   - Crustal model at 5x5deg blocks (w/ topography)
%                     (http://igppweb.ucsd.edu/~gabi/crust1.html)
%
%    Notes:
%     - SRTM30+ elevation is read from the "topo30" binary file.
%     - ETOPO1* elevation is read from the "etopo1_*_g_i2.bin" binary file.
%     - CRUST*.* uses GETCRUST to get the elevation values.
%     - For either SRTM30+ or ETOPO1 you need to put the binary files on
%       the path (download them from the websites above).  For CRUST*.* the
%       zip file seizmo_3d_models.zip downloaded by INSTALL_SEIZMO must
%       have been expanded on the path.
%
%    Examples:
%     % Plot up an equatorial profile:
%     lon=-180:0.1:180;
%     figure;
%     plot(lon,topo_points(0,lon))
%
%     % Compare equatorial profiles:
%     lon=-180:0.1:180;
%     figure;
%     plot(lon,topo_points(0,lon,'srtm30+'),'r',...
%          lon,topo_points(0,lon,'etopo1_bed'),'g--',...
%          lon,topo_points(0,lon,'etopo1_ice'),'b-.',...
%          lon,topo_points(0,lon,'crust1.0'),'k:',...
%          lon,topo_points(0,lon,'crust2.0'),'y',...
%          lon,topo_points(0,lon,'crust5.1'),'m--')
%     legend({'srtm30+' 'etopo1_bed' 'etopo1_ice' ...
%             'crust1.0' 'crust2.0' 'crust5.1'},'interpreter','none');
%
%    See also: TOPO_REGION, TOPO_COLORMAP

%     Version History:
%        Feb. 16, 2010 - initial version
%        May  17, 2010 - minor doc touch
%        May  19, 2010 - return nothing when given nothing
%        May  20, 2010 - added scrollbar (cause it is slow)
%        Feb. 11, 2011 - mass nargchk fix
%        Apr.  2, 2012 - minor doc update
%        Jan. 28, 2014 - read from int16 binary files directly, drop
%                        etopo1_thk (just do two calls and do the math),
%                        added CrustX.X by redirecting to GETCRUST
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 28, 2014 at 15:05 GMT

% todo:

% check nargin
error(nargchk(2,3,nargin));

% valid topo list (add new topo here)
validtopo={'srtm30+' 'etopo1_bed' 'etopo1_ice' ...
    'crust1.0' 'crust2.0' 'crust5.1'};

% check inputs
if(nargin==2 || isempty(toponame)); toponame='srtm30+'; end
if(~isreal(lat) || ~isreal(lon))
    error('seizmo:topo_points:badInput',...
        'LAT/LON must be real-valued!');
elseif(~isscalar(lat) && ~isscalar(lon) ...
        && ~isequal(size(lat),size(lon)))
    error('seizmo:topo_points:badInput',...
        'LAT/LON must be scalar or equal sized!');
elseif(~any(strcmpi(toponame,validtopo)))
    error('seizmo:topo_points:badInput',...
        ['TOPONAME must be one of the following:\n' ...
        sprintf('%s ',validtopo{:})]);
end

% expand scalars
if(isscalar(lat)); lat=lat(ones(size(lon))); end
if(isscalar(lon)); lon=lon(ones(size(lat))); end
npts=numel(lat);

% quick return if no points
if(~npts); topo=[]; return; end

% take care of wrap-around
[lat,lon]=fixlatlon(lat,lon);

% proceed by topo
switch lower(toponame)
    case 'srtm30+'
        type='cell';
        rows=21600;
        cols=43200;
        slat=90;
        latsign=-1;
        slon=0;
        lonsign=1;
        inclat1st=false;
        bytesize=2;
        bytetype='int16';
        endian='ieee-be';
        width=1/120;
    case {'etopo1_bed' 'etopo1_ice'}
        type='grid';
        rows=10801;
        cols=21601;
        slat=90;
        latsign=-1;
        slon=-180;
        lonsign=1;
        inclat1st=false;
        bytesize=2;
        bytetype='int16';
        endian='ieee-le';
        width=1/60;
    case {'crust1.0' 'crust2.0' 'crust5.1'}
        topo=getcrust(lat,lon,toponame);
        topo=reshape(topo.top(:,2),size(lat)).*1000;
        return;
end

% offset cell registration
if(strcmp(type,'cell'))
    slat=slat+width/2*latsign;
    slon=slon+width/2*lonsign;
end

% convert lat/lon to nearest row/column
row=round(latsign.*(lat-slat)/width)+1;
col=round(lonsign.*(lon-slon)/width)+1;

% deal with pole issues (rounding issue for cell registration)
if(any(row==0)); row(row==0)=1; end
if(any(row==rows+1)); row(row==rows+1)=rows; end

% deal with negative columns
if(any(col<=0))
    switch type
        case 'cell'
            col(col<=0)=cols+col(col<=0);
        case 'grid'
            col(col<=0)=cols+col(col<=0)-1;
    end
end

% deal with cols+1 values (rounding issue for cell registration)
if(any(col==cols+1)); col(col==cols+1)=cols; end

% binary filename
switch toponame
    case 'srtm30+'
        filename='topo30';
    case 'etopo1_ice'
        filename='etopo1_ice_g_i2.bin';
    case 'etopo1_bed'
        filename='etopo1_bed_g_i2.bin';
end

% open file
fid=fopen(filename,'r',endian);

% fid check
if(fid<0)
    % bad permissions?
    error('seizmo:topo_region:badFID',...
        'File not openable, %s !',filename);
end

% going by rows or columns?
topo=nan(size(lat));
if(inclat1st)
    % seek locations
    seek=bytesize.*(rows.*(col-1)+row);
    
    % optimize their order for a speed boost
    [seek,idx]=sort(seek(:));
    
    % goes across latitudes first
    for i=1:numel(col)
        fail=fseek(fid,seek(i),'bof');
        if(fail)
            error('seizmo:topo_region:badSEEK',...
                'Seek failed on %s !',filename);
        end
        topo(idx(i))=fread(fid,1,bytetype);
    end
else
    % seek locations
    seek=bytesize.*(cols.*(row-1)+col);
    
    % optimize their order for a speed boost
    [seek,idx]=sort(seek(:));
    
    % goes across longitudes first
    for i=1:numel(col)
        fail=fseek(fid,seek(i),'bof');
        if(fail)
            error('seizmo:topo_region:badSEEK',...
                'Seek failed on %s !',filename);
        end
        topo(idx(i))=fread(fid,1,bytetype);
    end
end

% close file
fclose(fid);

end
