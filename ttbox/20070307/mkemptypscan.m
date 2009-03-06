function pscan=mkemptypscan(cnt);
% MKEMPTYPSCAN......create empty pscan structure for MKSCANP
%
% call: pscan=mkemptypscan(cnt);
%
%             cnt: number of phases to be scanned
%
% result: pscan: [struct]
%                a structure to contain the output of a MKSCANP call. Ths
%                structure has the following fields:
%
%
%                .phase: [string]
%                        Name of the phase for which this struct is valid
%                .h: [km]
%                    focal depth for which this struct is valid
%                .angles: [deg]
%                         list of take off angles
%                .p: [sec/deg]
%                    list of ray parameters correpsonding to the take off
%                    angles at depth h
%                .dist: [deg]
%                    list of epicentral distances corresponding to p.
%                .vp: [km/s]
%                    P wave velocity at focal depth
%                .vs: [km/s]
%                    S wave velocity at focal depth
%                .starts: [index]
%                    positions where continuous pieces of dist(p) begin
%                .ends: [index]
%                    positions where continuous pieces of dist(p) end
%                    The i-th continuous piece begins at starts(i) and ends at
%                    ends(i).
%
% Martin Knapmeyer, 30.05.2005

pscan=repmat(struct('phase','',...
                    'h',0,...
                    'angles',[],...
                    'p',[],...
                    'dist',[],...
                    'vp',0,...
                    'vs',0,...
                    'starts',[],...
                    'ends',[]),...
             cnt,1);
