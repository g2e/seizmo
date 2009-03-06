function newclr=mksortclr(oldclr);
% mksortclr........sort layers in CLR structure by depth
%
% call: newclr=mksortclr(oldclr);
%
%       oldclr: CLR structure as returned by MKMAKECLR
%               in which layers are in no specific order
%
% result: newclr: CLR structure as returned by MKMAKECLR, representing
%                 the same velocity structure as OLDCLR but with layers
%                 ordered by depth.
%
% Martin Knapmeyer, 10.11.2003


%%% init result
newclr=oldclr;


%%% collect upper boundary depths of all layers
upperdep=zeros(newclr.lyrcnt,1);
for indy=1:newclr.lyrcnt
    upperdep(indy)=min(newclr.layers(indy).depth);
end; % for indy

%%% sort layers according to upper boudnary depth
[upperdep,sorter]=sort(upperdep);
newclr.layers=newclr.layers(sorter);
