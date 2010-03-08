function [li1,li2,nr,M]=ttrefine(m,std,lag,z,minstd)
%TTREFINE   Rearranges cross-correlation peak info

% check nargin
msg=nargchk(3,4,nargin);
if(~isempty(msg)); error(msg); end

% check inputs

% find misfit
% - weighted by 1/z^2
% - scales with stddev^2 but with a plateau under minstd standard
%   deviations (says we care little about the timing unless it is
%   outside some number of standard deviations)
if(vector)
    [i,j]=ind2sub([nr nr],find(tril(true(nr),-1)));
    M=(1./z.^2).*max(minstd^2,((lv+m(i,1,ones(1,np))-m(j,1,ones(1,np)))./(std(i,1,ones(1,np))+std(j,1,ones(1,np)))).^2);
    factor=1;
else
    
    factor=2;
end

end
