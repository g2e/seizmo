function [newdepths,criticalz]=mkderivsmp(oldmodel);
% mkderivsmp......analyze first derivative of v(z) to find triplications
%
% call: [newdepths,prayp,srayp]=mkderivsmp(oldmodel);
%
%       model: MODEL structure as returned by MKREADND or MKCLR2MODEL
%
% result: newdepths: depths of the depth samples to be added to the model [km]
%                    For each discontinuity, more than one (namely two) new
%                    depth samples are generated.
%
%
%         criticalz: list of the critical depths [km]
%                    This list contains all the critical depths: the actual
%                    depth samples of the discontinuities themselves.
%                    Critical rays can must derived from these.
%
% Strong changes in the velocity gradient produce triplications in the
% travel time. The beginning (=highest ray parameter) of these
% triplications can be derived from the velocity model, whereas their
% endings (=smallest ray parameter) depend on seismic phase and focal
% depth.
% The analysis has to be done in flat earth corrdinates.
%
% The new depths generated here are intended for insertion into the
% velocity sampling by a call of MKINTERPMODEL.
%
% Martin Knapmeyer, 23.06.2003



%%% computing the velocity derivatives produces divide-by-zero warnings at
%%% each discontinuity. We do not want these.
warning off MatLab:divideByZero


%%%%%% zeroth step: some initializations
criticalz=[]; % to collect critical depths detected in spherical co'
newdepths=[]; % to collect all extra samples


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               FLAT EARTH CORRDINATES SECTION                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% the analysis should be done in spherical coordinates and
%%% another part in flat earth coordinates. using the COORDSYS flag, the
%%% flat earth transformation can be switched off for testing.
%%% This option was introduced for testing purposes, but remained in the code.
coordsys='flat';


%%%%%% FET: transformation into flat earth (can be dactivated)
%%%%%%      since TTBOX does its ray tracing in flat earth, the properties of
%%%%%%      the flat earth velocity distribution have to be analyzed!
switch coordsys
    case {'flat'}
        [vp]=mksfer2flat(oldmodel.vp,oldmodel.z,oldmodel.rp);
        [vs,z]=mksfer2flat(oldmodel.vs,oldmodel.z,oldmodel.rp); %note that depth sampling is same as for Vp
    case {'spherical'}
        vp=oldmodel.vp;
        vs=oldmodel.vs;
        z=oldmodel.z;
end; % switch coord


%%%%%% begin of mdoel analysis in flat coordinates
             
%%%%%% third step: identify depths at which the velocity gradient changes,
%%%%%%             rapidly. This has to be done for Vp and Vs separately.
[dvp,dz]=mkderive(z,vp,1); % derivative of Vp
[dvs]=mkderive(z,vs,1); % derivative of Vs


%%%%%% fourth step: identify extrema of velocity derivative.
%%%%%%              depth samples will be added if either the P or the S
%%%%%%              velocity derivative has an extremum.
%%%%%%              To be able to treat the first and last sample without
%%%%%%              special code, two zeros are appended.
dvp=abs([0; dvp; 0]);
dvs=abs([0; dvs; 0]);
dz=[0; dz; 0];
dvpdz=dvp./dz; % P velocoty derivative
dvsdz=dvs./dz; % S velocity derivative
derivlen=length(dz); % number of samples in the derivatives
dvextrema=[];
dvextremextension=[];


%%% identify jumps in the derivative
%%% we have to look for upward and downward jumps in both P and S velocity
% pepsilon=1.05;
% pratio=dvpdz(1:(end-1))./dvpdz(2:end);
% sratio=dvsdz(1:(end-1))./dvsdz(2:end);
% 
% % pindy=find((pratio~=0)&(~isinf(pratio))&((pratio>pepsilon)|(pratio<1/pepsilon)));
% % sindy=find((sratio~=0)&(~isinf(sratio))&((sratio>pepsilon)|(sratio<1/pepsilon)));
% pindy=find(((pratio>pepsilon)|(pratio<1/pepsilon)));
% sindy=find(((sratio>pepsilon)|(sratio<1/pepsilon)));
% 
% dvextrema=[dvextrema; dz(unique([pindy; sindy]))];
% 
% figure(17);
% clf;
% plot(dz(1:(end-1)),pratio,'.-');
% hold on
% plot(dz(1:(end-1)),sratio,'r.-');
% plot(dz(pindy),pratio(pindy),'mv');
% plot(dz(sindy),sratio(sindy),'m^');
% hold off





for indy=2:derivlen-1 % loop only over original samples, not the extra zeros
      
    
    isextremum=((mkisextreme(dvp((-1:1)+indy))==1) || (mkisextreme(dvs((-1:1)+indy))==1));
    if isextremum
       %%% the current sample of either d2vp or d2vs is larger than its two neighbours
       %%% put it into the collection of critical depths
       dvextrema=[dvextrema; dz(indy)];
       
       %disp(num2str([dz(indy) mkdepth2rayp(oldmodel,dz(indy)) dvp(indy)/dvp(indy-1) dvp(indy-1)/dvp(indy) dvs(indy)/dvs(indy-1) dvs(indy-1)/dvs(indy)]));
       
    end; % if (dvp(indy)>dvp(indy-1))&&(dvp(indy)>dvp(indy+1))...
    
    
%         figure(17);
%         hold on
%         rayp=mkdepth2rayp(oldmodel,dz(indy));
%         handle=plot([1 1 NaN 1 1]*rayp,...
%              [dvp(indy)/dvp(indy-1) dvp(indy-1)/dvp(indy),...
%               NaN,...
%               dvs(indy)/dvs(indy-1) dvs(indy-1)/dvs(indy)],...
%              '.-');
% %         handle=plot([dvp(indy)/dvp(indy+1) dvp(indy)/dvp(indy-1),...
% %                      NaN,...
% %                      dvs(indy)/dvs(indy+1) dvs(indy)/dvs(indy-1)],...
% %                      [1 1 NaN 1 1]*dz(indy-1),'.-');
%         if isextremum
%             set(handle,'Color','g');
%         else
%             set(handle,'Color','r');
%         end;
%         hold off
      
end; % for indy

%%% append new samples to samples lists
%disp(['MKDERIVSMP: ' int2str(length(dvextrema))...
%      ' critical depths from extrema of velocity derivative identified.' ]);
criticalz=[criticalz; dvextrema];
%dvextremextension=[dvextrema-zepsilon; dvextrema+zepsilon]; adding new depth samples not necessary in this case!!


%%%%%% end of analysis in flat coordinates


%%%%%% INVERSE FET: transform the identified depths back into spherical
%%%%%%              earth coordinates
flatnewz=[dvextremextension];
switch coordsys
    case {'flat'}
        [dmy,flatnewz]=mkflat2sfer(flatnewz,flatnewz,oldmodel.rp);
        [dmy,criticalz]=mkflat2sfer(criticalz,criticalz,oldmodel.rp);
    case {'spherical'}
        %%% nothing happens, since the depths are still in spherical coords!
end; % switch coord



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               END OF FLAT EARTH CORRDINATES SECTION                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% switch the warning on againg - may be we want to get warnings from
%%% other parts of TTBOX.
warning on MatLab:divideByZero
