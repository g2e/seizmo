function [topo,lat,lon]=topo_region(latrange,lonrange,toponame)
%TOPO_REGION    Returns topography data for the region indicated
%
%    Usage:    [topo,lat,lon]=topo_region(latrange,lonrange)
%              [topo,lat,lon]=topo_region(latrange,lonrange,toponame)
%
%    Description:
%     [TOPO,LAT,LON]=TOPO_REGION(LATRANGE,LONRANGE) returns a matrix of the
%     topography for the region bounded by LATRANGE & LONRANGE.  LATRANGE &
%     LONRANGE must in degrees and be 1x2 real-valued arrays.  LATRANGE
%     values should be between +/-90deg, LONRANGE values between +/-180deg.
%     If LONRANGE(1)>LONRANGE(2) the region extends the long way around the
%     Earth.  The topography is SRTM30+ and is in meters.
%
%     [TOPO,LAT,LON]=TOPO_REGION(LATRANGE,LONRANGE,TOPONAME) allows setting
%     the dataset that the topography (z-values) is pulled from.  TOPONAME
%     must be a string and the current possibilities are:
%      'srtm30+'    - SRTM30+ topography (30 arc-second resolution)
%                     (http://topex.ucsd.edu/WWW_html/srtm30_plus.html)
%      'etopo1_bed' - ETOPO1 bedrock topography (1 arc-minute resolution)
%      'etopo1_ice' - ETOPO1 ice topography (1 arc-minute resolution)
%                     (http://www.ngdc.noaa.gov/mgg/global/global.html)
%      'crust1.0'   - Crustal model at 1x1deg blocks (w/ topography)
%      'crust2.0'   - Crustal model at 2x2deg blocks (w/ topography)
%      'crust5.1'   - Crustal model at 5x5deg blocks (w/ topography)
%                     (http://igppweb.ucsd.edu/~gabi/crust1.html)
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
%     % The Bering Sea:
%     [topo,lat,lon]=topo_region([64 49],[162 194],'etopo1_bed');
%     figure; imagesc(lon,lat,topo);
%     topo_colormap(topo); axis equal tight xy;
%
%    See also: TOPO_POINTS, TOPO_COLORMAP

%     Version History:
%        Feb. 17, 2010 - initial version
%        May  17, 2010 - slight improvement to example
%        Feb. 11, 2011 - mass nargchk fix
%        Mar. 24, 2012 - minor doc update
%        Jan. 28, 2014 - read from int16 binary files directly, drop
%                        etopo1_thk (just do two calls and do the math),
%                        added CrustX.X by redirecting to GETCRUST
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 28, 2014 at 15:05 GMT

% todo:

% check nargin
error(nargchk(0,3,nargin));

% valid topo list
validtopo={'srtm30+' 'etopo1_bed' 'etopo1_ice' ...
    'crust1.0' 'crust2.0' 'crust5.1'};

% check inputs
if(nargin==2 || isempty(toponame)); toponame='srtm30+'; end
if(~isreal(latrange) || ~isreal(lonrange))
    error('seizmo:topo_region:badInput',...
        'LATRANGE/LONRANGE must be real-valued!');
elseif(~isequal(size(latrange),[1 2]) ...
        && ~isequal(size(lonrange),[1 2]))
    error('seizmo:topo_region:badInput',...
        'LATRANGE/LONRANGE must be 1x2!');
elseif(any(abs(latrange)>90))
    error('seizmo:topo_region:badInput',...
        'LATRANGE must be in the range +/-90deg!');
elseif(~any(strcmpi(toponame,validtopo)))
    error('seizmo:topo_region:badInput',...
        ['TOPONAME must be one of the following:\n' ...
        sprintf('%s ',validtopo{:})]);
end

% take care of wrap-around
longway=~isequal(lonrange,sort(lonrange));
[latrange,lonrange]=fixlatlon(latrange,lonrange);
latrange=sort(latrange);

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
    case 'crust1.0'
        step=1;
    case 'crust2.0'
        step=2;
    case 'crust5.1'
        step=5;
end

% CrustX.X shortcut
switch lower(toponame)
    case {'crust1.0' 'crust2.0' 'crust5.1'}
        latrange=fix(latrange/step)*step+sign(latrange)*step/2;
        latrange(latrange>90)=latrange(latrange>90)-step;
        latrange(latrange<-90)=latrange(latrange<-90)+step;
        lat=latrange(1):step:latrange(2);
        lonrange=fix(lonrange/step)*step+sign(lonrange)*step/2;
        if(longway)
            if(lonrange(1)>lonrange(2))
                lon=lonrange(2):step:lonrange(1);
            else
                lon=lonrange(2):step:lonrange(1)+360;
            end
        else
            if(lonrange(1)>lonrange(2))
                lon=lonrange(1):step:lonrange(2)+360;
            else
                lon=lonrange(1):step:lonrange(2);
            end
        end
        [lon0,lat0]=meshgrid(lon,lat);
        topo=getcrust(lat0,lon0,toponame);
        topo=reshape(topo.top(:,2),size(lat0)).*1000;
        return;
end

% offset cell registration
if(strcmp(type,'cell'))
    slat=slat+width/2*latsign;
    slon=slon+width/2*lonsign;
end

% convert box limits to nearest row/column
rowrange=sort(round(latsign.*(latrange-slat)/width)+1);
colrange=round(lonsign.*(lonrange-slon)/width)+1;

% deal with pole issues (rounding issue for cell registration)
if(any(rowrange==0)); rowrange(rowrange==0)=1; end
if(any(rowrange==rows+1)); rowrange(rowrange==rows+1)=rows; end

% deal with negative columns
if(any(colrange<=0))
    switch type
        case 'cell'
            colrange(colrange<=0)=cols+colrange(colrange<=0);
        case 'grid'
            colrange(colrange<=0)=cols+colrange(colrange<=0)-1;
    end
end

% deal with cols+1 values (rounding issue for cell registration)
if(any(colrange==cols+1)); colrange(colrange==cols+1)=cols; end

% get lat/lon output
lat=slat+latsign.*((rowrange(1):rowrange(2))-1).*width;
lon=slon+lonsign.*colrange.*width;
if(longway)
    if(lon(2)>lon(1))
        lon=lon(2):width:lon(1)+360;
    else
        lon=lon(2):width:lon(1);
    end
else
    if(lon(2)>lon(1))
        lon=lon(1):width:lon(2);
    else
        lon=lon(1):width:lon(2)+360;
    end
end

% sort colrange based on ordering
if(longway)
    if(colrange(1)>colrange(2))
        % one read
        colrange=colrange([2 1]);
    else
        % two read
        colrange=[colrange(2) cols; 1 colrange(1)];
        if(strcmp(type,'grid')); colrange(2,1)=2; end
    end
else
    if(colrange(1)>colrange(2))
        % two read
        colrange=[colrange(1) cols; 1 colrange(2)];
        if(strcmp(type,'grid')); colrange(2,1)=2; end
    else
        % one read
        %colrange=colrange;
    end
end

% binary filename
switch toponame
    case 'srtm30+'
        filename='topo30';
    case 'etopo1_ice'
        filename='etopo1_ice_g_i2.bin';
    case 'etopo1_bed'
        filename='etopo1_bed_g_i2.bin';
end

% debug
%numel(lat)
%numel(lon)
%[rowrange rows]
%colrange
%cols
%die

% open file
fid=fopen(filename,'r',endian);

% fid check
if(fid<0)
    % bad permissions?
    error('seizmo:topo_region:badFID',...
        'File not openable, %s !',filename);
end

% going by rows or columns?
topo=cell(1,size(colrange,1));
if(inclat1st)
    % goes across latitudes first
    for i=1:size(colrange,1)
        total=(diff(rowrange)+1).*(diff(colrange(i,:))+1);
        strips=diff(rowrange)+1;
        skips=bytesize*(rows-rowrange(2)+rowrange(1)-1);
        fail=fseek(fid,...
            bytesize*(rows*(colrange(i,1)-1)+rowrange(1)),'bof');
        if(fail)
            error('seizmo:topo_region:badSEEK',...
                'Seek failed on %s !',filename);
        end
        topo{i}=fread(fid,total,[num2str(strips) '*' bytetype],skips);
        
        % reshape/flip/permute
        topo{i}=reshape(topo{i},[strips diff(colrange(i,:))+1]);
        if(latsign==1); topo{i}=flipud(topo{i}); end
        if(lonsign==-1); topo{i}=fliplr(topo{i}); end
    end
    
    % combine
    if(numel(topo)>1); topo=cat(2,topo{:}); else topo=topo{1}; end
else
    % goes across longitudes first
    for i=1:size(colrange,1)
        total=(diff(rowrange)+1).*(diff(colrange(i,:))+1);
        strips=diff(colrange(i,:))+1;
        skips=bytesize*(cols-colrange(i,2)+colrange(i,1)-1);
        fail=fseek(fid,...
            bytesize*(cols*(rowrange(1)-1)+colrange(i,1)),'bof');
        if(fail)
            error('seizmo:topo_region:badSEEK',...
                'Seek failed on %s !',filename);
        end
        topo{i}=fread(fid,total,[num2str(strips) '*' bytetype],skips);
        
        % reshape/flip/permute
        topo{i}=reshape(topo{i},[strips diff(rowrange)+1]);
        topo{i}=topo{i}.';
        if(latsign==1); topo{i}=flipud(topo{i}); end
        if(lonsign==-1); topo{i}=fliplr(topo{i}); end
    end
    
    % combine
    if(numel(topo)>1); topo=cat(2,topo{:}); else topo=topo{1}; end
end

% close file
fclose(fid);

end
