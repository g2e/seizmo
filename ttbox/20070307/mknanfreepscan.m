function pscannew=mknanfreepscan(pscan);
% mknanfreepscan.......remove NaN and inf entries from pscan structure
%
% call: pscannew=mknanfreepscan(pscan);
%
%       pscan: pscan structure as described in the hel lines of MKPSAMPLER
%
% result: pscannew: as pscan, but entries with NaN or inf distance are
%                   removed.
%
% Martin Knapmeyer, 19.12.2006

%%% init result
pscannew=pscan;

%%% test for existence of NaN and inf entries
nandists=find((isnan(pscannew.dist))|(isinf(pscannew.dist)));

%%% if there are any, remove them
if ~isempty(nandists)
   %disp(['MKNANFREEPSCAN: PSCAN contains ' int2str(length(nandists))...
   %      ' NaN or inf distances for phase ' pscan.phase ' at h='...
   %      num2str(pscan.h) 'km source depth!']);
     
   %%% remove the NaN elements
   nonnandists=find((~isnan(pscannew.dist))&(~isinf(pscannew.dist)));
   pscannew.dist=pscannew.dist(nonnandists);
   pscannew.angles=pscannew.angles(nonnandists);
   pscannew.p=pscannew.p(nonnandists);
   %disp(['MKNANFREEPSCAN: PSCAN elements with NaN or inf distances removed.']);
   
end; % if ~isempty(nandists)

return; 