function [newv,newvsec,newz,criticalz]=mklvzanalysis(v,vsec,z,rp);
% mklvzanalysis......analyze velocoty profile for low velocity zones
%
% call: [newv,newvsec,newz,criticalz]=mklvzanalysis(v,z,rp);
%
%       v: velocity as function of depth [km/s]
%       vsec: secondary velocity profile [km/s]
%             only the LVZ in V are analyzed. From VSEC, new samples for
%             the other velocity are constructed.
%       z: depths at which velocity is defined [km]
%          v(i) is the velocity at depth z(i)
%      rp: planetary radius [km]
%          needed for falt earth transformation
%
% result: newv: velocity at new depth samples [km/s]
%               These values are not necesarily identical to those
%               MKINTERPMODEL would find, since we interpolate in flat
%               earth coordinates here!
%         newvsec: secondary velocity at new depth sample [km/s]
%                  These values are not necesarily identical to those
%                  MKINTERPMODEL would find, since we interpolate in flat
%                  earth coordinates here!
%         newz: depth of new depth samples [km]
%         criticalz: critical depths: these depths need special
%                    consideration in the shooting process. [km]
%
%
% This routine analyzes one velocity profile (either Vp or Vs) for low
% velocity zones (LVZ). it determines the depth and velocity at the top and
% bottom of the LVZ and constructs new depth samples that should be
% included into the model data structure in order to allow for better ray
% construction. It also constructs a list of critical depths and ray
% parameters that need to be considered in the shooting process.
%
% Only one velocity profile (Vp OR Vs) is analyzed at a time. But since
% both velocities are needed to complete the depth sample, VSEC is the
% OTHER velocity profile.
% You therefore have to carry out two calls:
% - mklvzanalysis(vp,vs,z,rp) analyzes the LVZ in Vp and constructs velocity
%   samples for Vp and Vs at the critical depths of the Vp profile.
% - mklvzanalysis(vs,vp,z,rp) analyzes the LVZ in Vs and constructs
%   velocoty samples for Vs and Vp at the critical depths of the Vs profile.
% Afterwards, you have to rake care that two resulting sets of depth
% samples are mixed in the right way.
%
% It is assumed that V and VSEC together represent both Vp and Vs. If V
% contains zero velocities, it is assumed that it represents Vs and the
% zeroes denote liquid zones (like the outer core of the earth). To analyze
% for LVZ, the zeroes are replaced by VSEC (i.e. Vp) velocities in order to
% obtain the behaviour of S-P-converted waves.
%
% Martin Knapmeyer, 20.07.2006, 27.09.2006, 28.09.2006



%%% a magical number
%%% vertical epsilon value for placing new samples around important parts
%%% of the velocity model
zepsilon=0.001; % [km]



%%%%%% some initializations
criticalz=[]; % to collect critical depths detected in spherical co'
newz=[]; % to collect depths of extra samples
newv=[]; % to collect velocities at depths given in NEWZ
newvsec=[]; % to collect secondary velocoties at depths given in NEWZ


%%%%%% zeroth step: test for presence of zeros in V. If there are zeros,
%%%%%%              replace them with the corresponding velocities of VSEC.
%%%%%%              The justification of this is: the P velocities of the
%%%%%%              Earth's outer core together with the S velocities of the 
%%%%%%              lower mantle form an LVZ for S waves. To detect and
%%%%%%              analyze this LVZ, it is necessary to imply an S-to-P
%%%%%%              conversion. This does no harm to the other analyses.
zerolist=find(v==0);
if ~isempty(zerolist)
    %%% zero velocities exist in V, hence we assume it represents Vs.
    %%% and replace the zeros by VSEC values.
    %disp('MKLVZANALYSIS: Zero velocities found in v(z), assuming S-P conversion.');
    v(zerolist)=vsec(zerolist);
end; % if ~isempty(v)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               FLAT EARTH CORRDINATES SECTION                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% the analysis must be done in flat coordinates 


%%%%%% FET: transformation into flat earth
[vf,zf]=mksfer2flat(v,z,rp); % main velocity profile
vsecf=mksfer2flat(vsec,z,rp); % secondary velocity profile


%%%%%% begin of model analysis in flat coordinates


%%%%%% first step: loop over velocity model to identify velocity inversions
%%%%%%             a velocity inversion begins at the first sample where
%%%%%%             velocity is smaller than at the sample immediately
%%%%%%             above, and ends at the first sample where velocity is
%%%%%%             larger than at the last sample above the LVZ.
%%%%%%             We have to look for Vp and Vs inversions separatey, but
%%%%%%             can do that in the same loop.
sampanz=length(z);
inlvz=0; % flag: if 1, we're in LVZ, 0 else
lastabove=[]; % to collect indices of last samples above LVZs
firstbelow=[]; % to collect indices of first samples below LVZs
vtop=[]; % to collect velocities at top of LVZs
vsectop=[]; % to collect secondary velocites at top of LVZ
for sampcnt=2:sampanz % LVZ can not begin at first sample!
    
    if (inlvz==0)&(vf(sampcnt)<vf(sampcnt-1))
       %%% we're currently not in an LVZ, but here starts one!
       inlvz=1; % LVZ started
       lastabove=[lastabove; sampcnt-1]; % record start position: last sample above
       vtop=[vtop; vf(sampcnt-1)]; % velocity on top of LVZ
       vsectop=[vsectop; vsecf(sampcnt-1)];
    end; % if (inlvz==0)&(vf(sampcnt)<vf(sampcnt-1))
    
    if (inlvz==1)&(vf(sampcnt)>=vtop(end))
       %%% we are currently within a LVZ, but it ends here!
       inlvz=0; % LVZ ended
       firstbelow=[firstbelow; sampcnt];
    end; % if  (inlvz==1)&(vf(sampcnt)>vtop(end))
    
end; % for sampcnt


%%%%%% second step: true LVZ depths
%%%%%% the depth of the first sample below the LVZ is too deep: what we really
%%%%%% need is the exact depth where the LVZ top velocity is reached again.
%%%%%% This can be found by simple proportional partitioning between the
%%%%%% identified first-sample-below and its predecessor in depth.


% figure(2);
% clf;
% plot(vf,zf,'.-');
% hold on
% plot(vsecf,zf,'.-');
% plot(vf([lastabove]),zf([lastabove]),'rv');
% plot(vf([firstbelow]),zf([firstbelow]),'r^');
% hold off
% axis ij
% grid on



%%% determine exact bottom depth of LVZ
bottomanz=length(firstbelow); % number of LVZ bottoms
bottomz=zeros(size(firstbelow))+NaN; % to collect true depths
bottomv=bottomz; % to collect true velocities at BOTTOMZ depths
bottomvsec=bottomz; % to collect true secondary velocities at BOTTOMZ depths
for lvzcnt=1:bottomanz
    
    %%% two points defining the velocity gradient across the LVZ bottom
    vbelow=vf(firstbelow(lvzcnt)); % V at first sample below LVZ
    vabove=vf(firstbelow(lvzcnt)-1); % V at deepest sample within LVZ
    vsecbelow=vsecf(firstbelow(lvzcnt)); % VSEC at first sample below LVZ
    vsecabove=vsecf(firstbelow(lvzcnt)-1); % VSEC at deepest sample within LVZ
    zbelow=zf(firstbelow(lvzcnt)); % depth of first sample below LVZ
    zabove=zf(firstbelow(lvzcnt)-1); % depth of deepest sample within LVZ
    
    
%     figure(2);
%     hold on
%     plot([vbelow vabove],[zbelow zabove],'ro');
%     hold off
    
    
    if zabove~=zbelow
        %%% LVZ bottom is not a first order discontionuity
        %%% interpolate exact depth of LVZ bottom from thw two velocities and
        %%% the velocity at the LVZ top (linear interpolation)
        bottomz(lvzcnt)=interp1([vabove vbelow],...
                                [zabove zbelow],...
                                vtop(lvzcnt),...
                                'linear');
        bottomv(lvzcnt)=vtop(lvzcnt);
    else
        %%% LVZ bottom is a first order discontinuity
        %%% VPBOTTOMZ remains NaN to avoid multiple samples at
        %%% discontinuity undersides
    end; % if zabove~=zbelow
    
%     figure(2);
%     hold on
%     plot(bottomv(lvzcnt),bottomz(lvzcnt),'k^');
%     hold off

    %%% primary velocity profile done
    
    %%% now interpolate the secodnary velocity
    %%% two points defining the velocity gradient across the LVZ bottom
    vsecbelow=vsecf(firstbelow(lvzcnt)); % V at first sample below LVZ
    vsecabove=vsecf(firstbelow(lvzcnt)-1); % V at deepest sample within LVZ
    
    
    if zabove~=zbelow
        %%% LVZ bottom depth is already known, just store the secondary 
        %%% bottom velocity
 
        %bottomvsec(lvzcnt)=vsectop(lvzcnt);
        bottomvsec(lvzcnt)=interp1([zabove zbelow],...
                                   [vsecabove vsecbelow],...
                                   bottomz(lvzcnt),...
                                   'linear');
       
    else
        %%% LVZ bottom is a first order discontinuity
        %%% VPBOTTOMZ remains NaN to avoid multiple samples at
        %%% discontinuity undersides
    end; % if zabove~=zbelow
    
    %%% secondary velocity profile done

end; % for lvzcnt


%%% remove NaN elements from bottoms list
remain=find(~isnan(bottomz));
bottomz=bottomz(remain);
bottomv=bottomv(remain);
bottomvsec=bottomvsec(remain);


%%%%%% third step: collect LVZ top depths
topz=zf(lastabove);


% %%%%%% 3.5th step: introduce new samples slightly above LVZ top MK10102006
% %%%%%%             these new samples have to be introduced in NEWZ and
% %%%%%%             its companions.
%% this was intended to include the ray that touches the top side of the
%% LVZ without entering it. But didnt work. I leave this issue open for the
%% moment, since I do not have the impression that it is very important.
%% The fuiture will show...
% abovetopz=topz-zepsilon;
% [originalz,sorter]=sort([zf(lastabove); zf(lastabove-1)]);
% originalvf=[vf(lastabove); vf(lastabove-1)];
% originalvf=originalvf(sorter);
% originalvfsec=[vf(lastabove); vf(lastabove-1)];
% originalvfsec=originalvfsec(sorter);
% abovetopv=interp1(originalz,...
%                   originalvf,...
%                   abovetopz,...
%                   'linear');
% abovetopvsec=interp1(originalz,...
%                   originalvfsec,...
%                   abovetopz,...
%                   'linear');
                            
              
              
%%%%%% fourth step: construct new depth samples (depth and velocity)
newz=[bottomz]; %[zf(lastabove);...
      %zf(firstbelow);...
      %bottomz];
newv=[bottomv]; %[vf(lastabove);...
      %vf(firstbelow);...
      %bottomv];
newvsec=[bottomvsec];
              

%%%%%% end of analysis in flat coordinates

% % figure(2);
% % clf;
% % plot(vf,zf,'.-');
% % axis ij;
% % hold on;
% % plot(vf,zf,'gv');
% % plot(vsecf,zf,'gv');
% % plot(newv,newz,'ro');
% % plot(newvsec,newz,'rs');
% % hold off


% figure(2);
% clf;
% plot(vf,zf,'y.-'); % primary velocity  profile
% axis ij;
% hold on;
% plot(vsecf,zf,'c-'); % secondary velocity profile
% %%% primary velocity LVZ marker
% plot(vf(lastabove),topz,'bv'); % last sample above LVZ
% plot(newv,newz,'b^'); % bottom of LVZ
% %%% secondary velocity LVZ marker
% %plot(vsecf(lastabove),topz,'rv');
% %plot(newvsec,newz,'r^');
% hold off



%%%%%% INVERSE FET: transform the back into spherical earth
%%% transform LVZ top and bottom depths
[dmy,topz]=mkflat2sfer(topz,topz,rp);
[dmy,bottomz]=mkflat2sfer(bottomz,bottomz,rp);
%%% transform velocities
dmy=newz; % needed after back transform of primary profile
[newv,newz]=mkflat2sfer(newv,newz,rp);
[newvsec]=mkflat2sfer(newvsec,dmy,rp);




% figure(2);
% clf;
% plot(v,z,'.-');
% axis ij;
% hold on;
% plot(vsec,z,'.-');
% plot(v(lastabove),topz,'gv');
% plot(vsec(lastabove),topz,'gv');
% plot(newv,newz,'r^');
% plot(newvsec,newz,'m^');
% hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               END OF FLAT EARTH CORRDINATES SECTION                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% shift bottom depths VERY slightly downwards to distinguish the bottom
%%% from the top in terms of (otherwise identical) ray parameters. This is
%%% necessary to be able to send a ray to the LVZ bottom at all because
%%% MKX4P uses only the shallowest possible turning depth! MK27092006
bottomz=bottomz+zepsilon;

%%%%%% collect all critical depths
criticalz=[topz; bottomz];



%%% return results