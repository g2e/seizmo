function [pf]=capavg_cmb_profiles_specific(pf,pts,r)
r=r*6371*pi/180;
delaz=cat(1,pf.delaz);
delaz1=delaz(1:2:end,:);
delaz2=delaz(2:2:end,:);
ev=cat(1,pf.ev);
[lat1,lon1]=sphericalfwd(ev(:,1),ev(:,2),...
    delaz1(:,1)-27,azmean([delaz1(:,2)';delaz2(:,2)'])');
[lat2,lon2]=sphericalfwd(ev(:,1),ev(:,2),...
    delaz2(:,1)-27,azmean([delaz1(:,2)';delaz2(:,2)'])');
[idx,idx]=gcarc_count([lat1 lon1],[lat2 lon2],true,pts,r);
idx=unique(cat(1,idx{:}));
pf=pf(idx);
end
