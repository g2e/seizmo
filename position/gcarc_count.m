function [iso,idx,az,azidx]=gcarc_count(eq,st,minor,pt,r,n)
%GCARC_COUNT    Great-circle arc histogram count
%
%    Usage:    iso=gcarc_count(eq,st,minor,pt,r)
%              [iso,idx]=gcarc_count(...)
%              [iso,idx,az]=gcarc_count(eq,st,minor,pt,r,n)
%              [iso,idx,az,azidx]=gcarc_count(eq,st,minor,pt,r,n)
%
%    Description:
%     ISO=GCARC_COUNT(EQ,ST,MINOR,PT,R) counts the number of great-circle
%     arcs given by the positions in EQ, ST & MINOR that come within R
%     kilometers of the points given in PT.  Thus ISO is a column vector
%     with as many elements as there are positions in PT.  EQ & ST provide
%     the start/end positions of the great-cicle-arcs and must be formatted
%     as 1x2 or Ax2 arrays of [LAT LON].  MINOR is a logical indicating if
%     the arc is a minor arc or a major arc ('long way round') and should
%     be scalar or an A-length vector (one entry per great-circle arc).  PT
%     must be a 1x2 or Bx2 array of [LAT LON].  R is the bin radius in
%     kilometers and must be scalar.  The default for R is 125km and may be
%     omitted.
%
%     [ISO,IDX]=GCARC_COUNT(...) returns a Bx1 cell array of indices
%     corresponding to the great-cicles-arcs counted for each point.  For
%     instance, EQ(IDX{3},:) & ST(IDX{3},:) give the great-circle-arcs
%     found near PT(3,:).
%
%     [ISO,IDX,AZ]=GCARC_COUNT(EQ,ST,MINOR,PT,R,N) also returns the local
%     azimuth counts for each location in PT as a BxN array AZ where N is
%     the number of azimuthal bins from 0-180 degrees.  The default for N
%     is 10 (18deg bins) and may be omitted.
%
%     [ISO,IDX,AZ,AZIDX]=GCARC_COUNT(EQ,ST,MINOR,PT,R,N) returns a BxN cell
%     array of indices corresponding to the great-circle-arcs counted for
%     each local azimuth bin at each point.
%
%    Notes:
%     - References:
%        Barmin, Ritzwoller, & Levshin 2001, PAG, Vol. 158, pp. 1351-1375
%         A Fast and Reliable Method for Surface Wave Tomography
%
%    Examples:
%     % Make a random dataset:
%     eq=randlatlon(100);
%     st=randlatlon(100);
%     pt=randlatlon(1000);
%     [iso,~,az]=gcarc_count(eq,st,true,pt,1000); % 1000km radius bins
%
%     % Plot the arcs & points:
%     ax=mmap;
%     hold(ax,'on');
%     [lat,lon]=gcarc2latlon(eq(:,1),eq(:,2),st(:,1),st(:,2));
%     lon=unwrap(lon*pi/180,[],2)*180/pi; % avoid wraparound streak
%     m_line(lon',lat','linewi',3);
%     m_line(pt(:,2),pt(:,1),'marker','o',...
%         'markerfacecolor','y','linestyle','none');
%     m_grid('linestyle','none','box','fancy','tickdir','out');
%     title(ax,'Great Circle Arcs');
%     hold(ax,'off');
%
%     % Now plot the isotropic count values (black=best, white=worst):
%     ax=mmap;
%     hold(ax,'on');
%     m_scatter(ax,pt(:,2),pt(:,1),[],...
%         (max(iso)-iso(:,[1 1 1]))/(max(iso)-min(iso)),...
%         'filled','markeredgecolor','k');
%     m_grid('linestyle','none','box','fancy','tickdir','out');
%     title(ax,'Great Circle Arc Density');
%     hold(ax,'off');
%
%     % Now plot the azimuthal distribution (black=uniform, white=poor):
%     ax=mmap;
%     hold(ax,'on');
%     m_scatter(ax,pt(:,2),pt(:,1),[],...
%         1-repmat(sum(az,2)./(10*max(az,[],2)),[1 3]),...
%         'filled','markeredgecolor','k');
%     m_grid('linestyle','none','box','fancy','tickdir','out');
%     title(ax,'Great Circle Arc Azimuthal Distribution');
%     hold(ax,'off');
%
%    See also: DEGDIST_FROM_GC, CLOSEST_POINT_ON_GC, HAVERSINE,
%              SPHERICALINV, HISTC

%     Version History:
%        Feb. 11, 2012 - initial version
%        Oct.  9, 2012 - fixed bug that included arcs 2*R from point,
%                        return indices for azimuth bins
%        Oct. 11, 2012 - add progress bar output
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 11, 2012 at 19:00 GMT

% todo:

% check nargin
error(nargchk(2,6,nargin));

% defaults
if(nargin<3 || isempty(minor)); minor=true; end % minor arc paths
if(nargin<4 || isempty(pt)); pt=[0 0]; end % equator
if(nargin<5 || isempty(r)); r=125; end % 125km radius spatial bins
if(nargin<6 || isempty(n)); n=10; end % 10 azimuth bins from 0-180deg

% force minor as column vector
minor=minor(:)';

% check inputs
[A(1),two(1)]=size(eq); [A(2),two(2)]=size(st); A(3)=size(minor,1);
[Bv,two(3)]=size(pt);
Av=unique(A(A~=1)); if(isempty(Av)); Av=1; end
if(any(two~=2))
    error('seizmo:gcarc_count:badInput',...
        'EQ, ST, PT must be Nx2 as [Lat Lon] !');
elseif(numel(Av)>1)
    error('seizmo:gcarc_count:badInput',...
        'EQ, ST, MINOR do not correspond to the same number of points!');
elseif(~isreal(eq) || ~isreal(st) || ~isreal(pt))
    error('seizmo:gcarc_count:badInput',...
        'EQ, ST, & PT must be real-valued arrays!');
elseif(~islogical(minor))
    error('seizmo:gcarc_count:badInput',...
        'MINOR must be a logical array!');
elseif(numel(r)~=1 || numel(n)~=1 || ~isreal(r) || ~isreal(n) ...
        || r<=0 || n<=0)
    error('seizmo:gcarc_count:badInput',...
        'R & N must be real-valued positive scalars!');
end

% fix lats & lons
[eq(:,1),eq(:,2)]=fixlatlon(eq(:,1),eq(:,2));
[st(:,1),st(:,2)]=fixlatlon(st(:,1),st(:,2));
[pt(:,1),pt(:,2)]=fixlatlon(pt(:,1),pt(:,2));

% expand scalars
if(A(1)==1); eq=eq(ones(Av,1),:); end
if(A(2)==1); st=st(ones(Av,1),:); end
if(A(3)==1); minor=minor(ones(Av,1),1); end

% arc distances
d=haversine(eq(:,1),eq(:,2),st(:,1),st(:,2));
r=r/(6371*pi/180);

% bins
bins=linspace(0,180,n+1);

% verbosity
verbose=seizmoverbose;

% this shortcut for no azi output and if Av<Bv
if(nargout<3 && Av<Bv)
    % detail message
    if(verbose); print_time_left(0,Av); end
    
    % loop over each gc
    iso=zeros(Bv,1);
    if(nargout>1); idx=cell(Bv,1); end
    for i=1:Av
        % find (pt) within (r) of great circle (eq,st)
        in=find(degdist_from_gc(eq(i,1),eq(i,2),st(i,1),st(i,2),...
            pt(:,1),pt(:,2))<=r);
        if(numel(in)==0); continue; end
        
        % get closest points (cp) on that great circle (eq,st)
        [cplat,cplon]=closest_point_on_gc(...
            eq(i,1),eq(i,2),st(i,1),st(i,2),pt(in,1),pt(in,2));
        
        % are closest points (cp) on the minor arc (eq,st)?
        dcpeq=haversine(eq(i,1),eq(i,2),cplat,cplon);
        dcpst=haversine(st(i,1),st(i,2),cplat,cplon);
        cponminor=abs(dcpeq+dcpst-d(i))<100*eps(d(i));
        
        % specified arcs (eq,st) within (r) of (pt)
        ok=false(size(in));
        ok(minor(i) & cponminor)=true;
        ok(~minor(i) & ~cponminor)=true;
        deq=haversine(eq(i,1),eq(i,2),pt(in,1),pt(in,2));
        dst=haversine(st(i,1),st(i,2),pt(in,1),pt(in,2));
        ok(deq<r | dst<r)=true;
        
        % trim
        in=in(ok);
        iso(in)=iso(in)+1;
        if(verbose); print_time_left(i,Av); end
        if(nargout<2); continue; end
        for j=1:numel(in); idx{in(j)}=[idx{in(j)}; i]; end
    end
    return;
end

% loop over each point in pt
iso=zeros(Bv,1);
if(nargout>1); idx=cell(Bv,1); end
if(nargout>2); az=zeros(Bv,n+1); end
if(nargout>3); azidx=cell(Bv,n); end
if(verbose); print_time_left(0,Bv); end
for i=1:Bv
    % find great circles (eq,st) within (r) of (pt)
    in=find(degdist_from_gc(eq(:,1),eq(:,2),...
        st(:,1),st(:,2),pt(i,1),pt(i,2))<=r);
    if(numel(in)==0); continue; end
    
    % get closest point (cp) on those great circles (eq,st)
    [cplat,cplon]=closest_point_on_gc(...
        eq(in,1),eq(in,2),st(in,1),st(in,2),pt(i,1),pt(i,2));
    
    % is closest point (cp) on the minor arc (eq,st)?
    dcpeq=haversine(eq(in,1),eq(in,2),cplat,cplon);
    dcpst=haversine(st(in,1),st(in,2),cplat,cplon);
    cponminor=abs(dcpeq+dcpst-d(in))<100*eps(max(d(in)));
    
    % specified arcs (eq,st) within (r) of (pt)
    ok=false(size(in));
    ok(minor(in) & cponminor)=true;
    ok(~minor(in) & ~cponminor)=true;
    deq=haversine(eq(in,1),eq(in,2),pt(i,1),pt(i,2));
    dst=haversine(st(in,1),st(in,2),pt(i,1),pt(i,2));
    ok(deq<r | dst<r)=true;
    
    % trim
    cplat=cplat(ok);
    cplon=cplon(ok);
    in=in(ok);
    
    % have # of great circle arcs (eq,st) within (r) of (pt)
    if(verbose); print_time_left(i,Bv); end
    iso(i)=numel(in);
    if(nargout<2); continue; end
    idx{i}=in;
    if(nargout<3); continue; end
    
    % get local azimuths
    % -> az(pt,cp)+90 (tangent!)
    [dcp,azcp]=sphericalinv(pt(i,1),pt(i,2),cplat,cplon);
    localaz=azcp+90;
    if(any(dcp==0))
        % cp is pt, use az(pt,st)
        same=dcp==0;
        [dst,azst]=sphericalinv(pt(i,1),pt(i,2),...
            st(in(same),1),st(in(same),2)); %#ok<ASGLU>
        localaz(same)=azst;
    end
    
    % force from 0-180
    localaz(localaz<0)=localaz(localaz<0)+180;
    localaz(localaz>180)=localaz(localaz>180)-180;
    
    % bin them
    [az(i,:),binmem]=histc(localaz,bins);
    
    % get indexing
    if(nargout<4); continue; end
    azidx{i,1}=in(binmem==1 | binmem==n+1); % exactly 180deg is in 1st bin
    for j=2:n; azidx{i,j}=in(binmem==j); end
end

% az(:,end) is the exactly 180deg count
if(nargout<3); return; end
az(:,1)=az(:,1)+az(:,end);
az(:,end)=[];

end


% resolution from resolution matrix R
% - row of R is a map
% - spatial resolution: fit with cone
%   - radius of base of cone is resolution sigma
%   - 2l where l is node spacing is best resolution (in their inversion)
%     - resolution smaller than 2l: sigma=2l
% - amplitude bias: multiply unit cylinder of radius sigma by row
%   - avg amplitude of output within sigma gives bias (they got +/-10%)
% - amplitude bias is better for reliability assessment
