function [area]=sph_poly_area(lat,lon,method,ellipsoid)
%SPH_POLY_AREA    Finds area of a polygon on a sphere
%
%    Usage:    area=sph_poly_area(lat,lon)
%              area=sph_poly_area(lat,lon,method)
%              area=sph_poly_area(lat,lon,method,[a f])
%
%    Description:
%     AREA=SPH_POLY_AREA(LAT,LON) calculates the fractional area on a
%     sphere for a polygon whose vertices are given by LAT/LON.  The method
%     used is based on the angles at each of the vertices.  LAT & LON may
%     be equal sized vectors or matrices.  In the case that LAT & LON are
%     matrices, each column is considered a separate polygon (so the areas
%     of multiple polygons with an equal number of vertices may be found).
%     The vertices of the polygon should proceed in a clockwise manner to
%     retreive the area within the polygon.  Counter-clockwise progression
%     will return the area outside a polygon.  A polygon's sides should not
%     cross (separate it into distinct polygons).
%
%     AREA=SPH_POLY_AREA(LAT,LON,METHOD) specifies the method to use in
%     finding the fractional area of a polygon on a sphere.  Available
%     methods are 'LINES', 'PLANES', & 'AZIMUTHS'.  The default method is
%     'AZIMUTHS'.  All methods have some drawbacks.  'LINES' uses a line
%     contour to integrate the area but suffers from inaccuracy if the
%     polygon is not well sampled (points on polygon should be within a few
%     kilometers for decent accuracy).  It always returns the internal area
%     of a polygon no matter which direction the vertices proceed.
%     'PLANES' does not require any extra sampling of the polygon if the
%     sides are great circle arcs but does fail for polygons whose internal
%     angles are greater than 180 degrees or for multiple vertices at the
%     same position.  It always returns the internal area of a polygon no
%     matter which direction the vertices proceed.  'AZIMUTHS' has the
%     benefit of the accuracy of the 'PLANES' method without any of its
%     failings (except that of not handling polygons who cross themselves).
%     The area returned for 'AZIMUTHS' depends on the progression of the
%     vertices, which is not like either of the other two methods.
%
%     AREA=SPH_POLY_AREA(LAT,LON,METHOD,[A F]) specifies the ellipsoid used
%     in the area calculation.  Specifying the ellipsoid causes AREA to be
%     in actual units (based on the units of the semimajor axis A).  F is
%     the flattening.
%
%    Notes:
%     - Are you getting weird areas?  Your polygon is probably crossing
%       itself (not allowed with the default method -- use 'lines' if you
%       want a reasonable area for such unreasonable polygons).  
%
%    Examples:
%     % Area of a polygon on the WGS-84 ellipsoid vs spherical earth:
%     sph_poly_area([0 90 0],[0 0 90],[],[6378.137 1/298.257223563])
%     sph_poly_area([0 90 0],[0 0 90],[],[6371 0])
%
%     % Compare the area from the default method with that of 'lines' (the
%     % true area is 0.125):
%     sph_poly_area([0 90 0],[0 0 90])
%     sph_poly_area([0 90 0],[0 0 90],'lines')
%
%    See also: SPH_POLY_IN, SPH_TRI_INIT, SPH_TRI_AUTO, SPH_TRI_SPLIT

%     Version History:
%        Nov. 16, 2009 - initial version
%        Feb.  7, 2011 - azimuths method works much better now
%        Feb. 18, 2011 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 18, 2011 at 19:30 GMT

% todo:

% check nargin
error(nargchk(2,4,nargin));

% check lat/lon
if(isvector(lat)); lat=lat(:); end
if(isvector(lon)); lon=lon(:); end
if(~isequal(size(lat),size(lon)))
    error('seizmo:sph_poly_area:badInput',...
        'LAT & LON must have the same size!');
end
if(~isreal(lat) || ~isreal(lon))
    error('seizmo:sph_poly_area:badInput',...
        'LAT & LON must be real valued!');
end
v=size(lat,1);

% default/check method
if(nargin==2 || isempty(method)); method='azimuths'; end
if(~ischar(method))
    error('seizmo:sph_poly_area:badInput','METHOD must be a string!');
end

% handle ellipsoid
if(nargin==4 && ~isempty(ellipsoid))
    if(~isreal(ellipsoid) || ~numel(ellipsoid)==2 || ellipsoid(2)>=1)
        error('seizmo:sph_poly_area:badEllipsoid',...
            ['Ellipsoid must be a 2 element vector specifying:\n'...
            '[equatorial_km_radius flattening(<1)]']);
    end
    ecc=sqrt(ellipsoid(2)*(2-ellipsoid(2)));
    lat=geographic2authaliclat(lat,ecc);
    radius=mean_ellipsoid_radius('authalic',ellipsoid);
end

% choose method
switch lower(method)
    case 'lines'
        % close the polygon
        lat(end+1,:)=lat(1,:);
        lon(end+1,:)=lon(1,:);

        % get relative locations
        [dist,az]=sphericalinv(0,0,lat,lon);

        % convert to radians
        D2R=pi/180;
        dist=dist*D2R;
        az=az*D2R;

        % width of rectangle (difference in azimuth)
        width=diff(az,1,1);

        % height of rectangle (average of distances)
        height=dist(1:end-1,:)+diff(dist,1,1)/2;

        % handle azimuth wrap-around
        width(width>pi)=width(width>pi)-2*pi;
        width(width<-pi)=width(width<-pi)+2*pi;

        % integrate
        % - fraction of unit sphere
        % - abs to handle direction around polygon
        % - force area to be the lesser
        area=abs(sum(width.*(1-cos(height)),1))/(4*pi);
        area=min(area,1-area);
    case 'planes'
        % wrap the polygon (for angles)
        lat=[lat(end,:); lat; lat(1,:)];
        lon=[lon(end,:); lon; lon(1,:)];
        
        % next dimension
        nd=ndims(lat)+1;
        
        % convert to xyz
        A=cat(nd,cosd(lon).*cosd(lat),sind(lon).*cosd(lat),sind(lat));
        
        % get perpendiculars to planes formed by sides of polygon
        B=cross(submat(A,1,1:v+1),submat(A,1,2:v+2),nd);
        
        % normalize perpendiculars
        N=vecnorm(B,nd);
        B=B./submat(N,nd,[1 1 1]);
        
        % get interior angles
        angles=pi-acos(dot(submat(B,1,1:v),submat(B,1,2:v+1),nd));
        
        % get area
        area=(sum(angles,1)-(v-2)*pi)/(4*pi);
    case 'azimuths'
        % loops are easy (leave me alone!)
        area=nan(1,size(lat,2));
        for i=1:size(lat,2)
            % remove duplicate points
            lat0=lat(:,i); lon0=lon(:,i);
            dist=sphericalinv(lat0,lon0,lat0([2:end 1]),lon0([2:end 1]));
            bad=dist<10*eps('single');
            lat0(bad)=[]; lon0(bad)=[];
            
            % 0 area cheats (numel = 0 or 2)
            if(numel(lat0)<3)
                area(1,i)=0;
                continue;
            end
            
            % wrap the polygon (for angles)
            lat0=[lat0(end); lat0; lat0(1)];
            lon0=[lon0(end); lon0; lon0(1)];
            
            % get azimuths
            [az,az]=sphericalinv(lat0(2:end-1),lon0(2:end-1),...
                lat0(3:end),lon0(3:end));
            [baz,baz]=sphericalinv(lat0(2:end-1),lon0(2:end-1),...
                lat0(1:end-2),lon0(1:end-2));
            
            % force az to be less than baz (requires clockwise progression)
            bad0=baz<az;
            if(any(bad0(:))); az(bad0)=az(bad0)-360; end
            
            % get area
            area(1,i)=(pi/180*sum(baz-az,1)-(v-2-sum(bad))*pi)/(4*pi);
        end
    otherwise
        error('seizmo:sph_poly_area:badMethod',...
            'Unknown METHOD: %s',method);
end

% convert to actual area if ellipsoid given
if(nargin==4 && ~isempty(ellipsoid))
    area=area*4*pi*radius*radius;
end

end
