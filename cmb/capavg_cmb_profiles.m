function [varargout]=capavg_cmb_profiles(pf,step,r)
%CAPAVG_CMB_PROFILES    Cap-average CMB profiles over a lat/lon grid
%
%    Usage:    images=capavg_cmb_profiles(pf,step,r)
%              [images,iso,isoidx]=capavg_cmb_profiles(pf,step,r)
%              capavg_cmb_profiles(...)
%
%    Description:
%     IMAGES=CAPAVG_CMB_PROFILES(PF,STEP,R) localizes CMB profile data in
%     PF (as returned by SLOWDECAYPAIRS) in a lat/lon grid with nodes every
%     STEP degrees and bins of radius R degrees.  IMAGES is a struct
%     with fields for the average of each measurement at every location,
%     the average of the error at every location, and the standard
%     deviation of the measurements at each location.  Also the fields
%     IMAGES.lat & IMAGES.lon allow mapping those fields using MMAP.
%
%     [IMAGES,ISO,ISOIDX]=CAPAVG_CMB_PROFILES(PF,STEP,R) also returns the
%     profile count and indexing for other operations.  See GCARC_COUNT for
%     more info.
%
%     CAPAVG_CMB_PROFILES(...) will make plots for all the measurements
%     that would be returned in IMAGES (see the above usage forms).
%
%    Notes:
%
%    Examples:
%     % Plot just the corrected slowness results:
%     im=capavg_cmb_profiles(pf,4,10);
%     mmap('image',{im.lat im.lon im.cslow},...
%          'fg','k','l',false,'o',false);
%
%     % No outputs will create maps for everything:
%     capavg_cmb_profiles(pf,4,10);
%
%    See also: SLOWDECAYPAIRS, SLOWDECAYPROFILES, MAP_CMB_PROFILES

%     Version History:
%        Oct. 16, 2012 - initial version
%        Nov.  2, 2012 - r is now in degrees
%        Nov. 23, 2012 - now capavg_...
%        Jan. 10, 2013 - documented no output plotting
%        Oct. 11, 2013 - minor example fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 11, 2013 at 13:30 GMT

% todo:

% check number of inputs
error(nargchk(3,3,nargin));

% check profiles struct
reqfields={'gcdist','azwidth','slow','slowerr','decay','decayerr',...
    'cslow','cslowerr','cdecay','cdecayerr','cluster','kname','st','ev',...
    'delaz','corrcoef','freq','phase','runname','dirname',...
    'time'};
if(~isstruct(pf) || any(~isfield(pf,reqfields)))
    error('seizmo:capavg_cmb_profiles:badInput',...
        ['PF must be a struct with the fields:\n' ...
        sprintf('''%s'' ',reqfields{:}) '!']);
end

% check inputs
if(~isscalar(step) || ~isnumeric(step) || ~isreal(step) || step<=0)
    error('seizmo:capavg_cmb_profiles:badInput',...
        'STEP must be a positive real-valued scalar in degrees!');
elseif(~isscalar(r) || ~isnumeric(r) || ~isreal(r) || r<=0)
    error('seizmo:capavg_cmb_profiles:badInput',...
        'R must be a positive real-valued scalar in degrees!');
end

% convert R to surface kilometers
r=r*6371*pi/180;

% make lat/lon grid
[im.lat,im.lon]=meshgrid(-90+step/2:step:90-step/2,...
    -180+step/2:step:180-step/2);

% extract info for profile location
delaz=cat(1,pf.delaz);
delaz1=delaz(1:2:end,:);
delaz2=delaz(2:2:end,:);
ev=cat(1,pf.ev);

% get profile location
[lat1,lon1]=sphericalfwd(ev(:,1),ev(:,2),...
    delaz1(:,1)-27,azmean([delaz1(:,2)';delaz2(:,2)'])');
[lat2,lon2]=sphericalfwd(ev(:,1),ev(:,2),...
    delaz2(:,1)-27,azmean([delaz1(:,2)';delaz2(:,2)'])');

% binning profiles
[iso,isoidx]=gcarc_count([lat1 lon1],[lat2 lon2],true,...
    [im.lat(:) im.lon(:)],r);

% map values
[im.slow,im.cslow,im.decay,im.cdecay]=deal(nan(size(im.lat)));
[im.slowerr,im.cslowerr,im.decayerr,im.cdecayerr]=deal(nan(size(im.lat)));
[im.slowstd,im.cslowstd,im.decaystd,im.cdecaystd]=deal(nan(size(im.lat)));
for i=1:numel(iso)
    im.slow(i)=mean(cat(1,pf(isoidx{i}).slow));
    im.cslow(i)=mean(cat(1,pf(isoidx{i}).cslow));
    im.decay(i)=mean(cat(1,pf(isoidx{i}).decay));
    im.cdecay(i)=mean(cat(1,pf(isoidx{i}).cdecay));
    im.slowerr(i)=sqrt(mean(cat(1,pf(isoidx{i}).slowerr).^2));
    im.cslowerr(i)=sqrt(mean(cat(1,pf(isoidx{i}).cslowerr).^2));
    im.decayerr(i)=sqrt(mean(cat(1,pf(isoidx{i}).decayerr).^2));
    im.cdecayerr(i)=sqrt(mean(cat(1,pf(isoidx{i}).cdecayerr).^2));
    im.slowstd(i)=std(cat(1,pf(isoidx{i}).slow));
    im.cslowstd(i)=std(cat(1,pf(isoidx{i}).cslow));
    im.decaystd(i)=std(cat(1,pf(isoidx{i}).decay));
    im.cdecaystd(i)=std(cat(1,pf(isoidx{i}).cdecay));
end

% output or plot?
if(nargout>0); varargout={im iso isoidx}; return; end

% plotting everything
ax(1)=mmap('image',{im.lat im.lon reshape(log10(iso),size(im.lat))},'fg','k','l',false,'o',false);
title(ax(1),'count'); colorbar('peer',ax(1));
ax(2)=mmap('image',{im.lat im.lon im.slow},'fg','k','l',false,'o',false);
title(ax(2),'slow'); colorbar('peer',ax(2));
ax(3)=mmap('image',{im.lat im.lon im.cslow},'fg','k','l',false,'o',false);
title(ax(3),'cslow'); colorbar('peer',ax(3));
ax(4)=mmap('image',{im.lat im.lon im.decay},'fg','k','l',false,'o',false);
title(ax(4),'decay'); colorbar('peer',ax(4));
ax(5)=mmap('image',{im.lat im.lon im.cdecay},'fg','k','l',false,'o',false);
title(ax(5),'cdecay'); colorbar('peer',ax(5));
ax(6)=mmap('image',{im.lat im.lon im.slowerr},'fg','k','l',false,'o',false);
title(ax(6),'slowerr'); colorbar('peer',ax(6));
ax(7)=mmap('image',{im.lat im.lon im.cslowerr},'fg','k','l',false,'o',false);
title(ax(7),'cslowerr'); colorbar('peer',ax(7));
ax(8)=mmap('image',{im.lat im.lon im.decayerr},'fg','k','l',false,'o',false);
title(ax(8),'decayerr'); colorbar('peer',ax(8));
ax(9)=mmap('image',{im.lat im.lon im.cdecayerr},'fg','k','l',false,'o',false);
title(ax(9),'cdecayerr'); colorbar('peer',ax(9));
ax(10)=mmap('image',{im.lat im.lon im.slowstd},'fg','k','l',false,'o',false);
title(ax(10),'slowstd'); colorbar('peer',ax(10));
ax(11)=mmap('image',{im.lat im.lon im.cslowstd},'fg','k','l',false,'o',false);
title(ax(11),'cslowstd'); colorbar('peer',ax(11));
ax(12)=mmap('image',{im.lat im.lon im.decaystd},'fg','k','l',false,'o',false);
title(ax(12),'decaystd'); colorbar('peer',ax(12));
ax(13)=mmap('image',{im.lat im.lon im.cdecaystd},'fg','k','l',false,'o',false);
title(ax(13),'cdecaystd'); colorbar('peer',ax(13));

end
