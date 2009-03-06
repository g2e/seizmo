function p=mkraydepthinv(rpd,v,r);
% mkraydepthinv....... find ray parameter by vertex depth
%
% call: p=mkraydepthinv(rpd,v,r);
%
%         rpd: ray penetration depth in terms of radius [km]
%           v: velocity distribution with radius [km/s]
%           r: radii at which velocities are defined [km]
%
% result: p: ray parameter at which a ray has its vertex at depth RPD [s/rad]
%            NaN if no such ray parameter exists.
%
% This is the inverse function to MKRAYDEPTH.
%
% The routine is designed specifically to find ray parameters for
% touching discontinuities of the model.
% It is therefore possible that it is less suitable for other purposes.
%
% Martin Knapmeyer 15.05.2002, 16.11.2006

%%% 16.11.2006: if RPD is NaN, NaN must be returned as result.
%%%             if find(r>=rpd) is empty, NaN is returned.


%%% init result
res=NaN;

if (~isnan(rpd))&(rpd~=0)
    %%% rpd is neither NaN nor zero, so we can start to compute something.

    %%% v and z have to be row vectors
    v=v(:)';
    r=r(:)';


    %%% interpolate velocity at vertex depth from model
    ident=find(r==rpd); % all layer boundaries identical to RPD
                        % length(ident) will be 2 if RPD is a discontinuity
                        %                       1 if RPD is a layer boundary
                        %                       0 if RPD is within a layer
    aboves=find(r>=rpd); % all layer boundaries above or at RPD
    if ~isempty(aboves)
        %%% depth samples above the core exist

        anz=length(aboves);
        %%% find index to top and bottom of layer to which RPD belongs
        switch length(ident)
           case {0}
              top=aboves(anz);
              bottom=aboves(anz)+1;
           case {1}
              top=aboves(anz-1);
              bottom=aboves(anz);
           case {2}
              top=aboves(anz-2);
              bottom=aboves(anz-1);
           otherwise
              error('MKRAYDEPTHINV: unexpected value of IDENT');
        end; % switch

        v1=interp1([r(top) r(bottom)],[v(top) v(bottom)],rpd);


        %%% compute the ray parameter
        if v1==0
           res=NaN;
        else
           res=rpd/v1;
        end; % if v1==0

    else
        %%% there are no depth samples above the core
        res=NaN;
    end; % if ~isempty(aboves)

end; % if (~isnan(rpd))&(rpd~=0)

%%% return result
p=res;