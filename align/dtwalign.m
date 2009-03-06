function [m,Gg]=dtwalign(cv,lv,xcpow,xcw,absw,abstt,absi)
%DTWALIGN    Gets weighted relative/absolute arrival times of signals
%
% INPUTS REQUIRED FOR RELATIVE ARRIVAL TIMES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cg/cv - matrix/vector of xc correlation values (used for weighting)
% lg/lv - matrix/vector of xc lag times
% xcpow - power(s) to apply to xc correlation values for weighting (0=no xc value weight)
% xcw   - coefficients to apply to xc correlation values for weighting
%
% ADDITIONAL INPUTS REQUIRED FOR TYING XC TO PICKS (ABSOLUTE TIMES)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% absw  - weight(s) for absolute travel times
% abstt - absolute travel times
% absi  - indices of records corresponding to absolue travel times


% FINDING SUBSCRIPT INDICES AND CONVERT GRIDS TO VECTORS
if(isvector(cv))
    len=length(cv);                                 % NUMBER OF CORRELATIONS
    nr=ceil(sqrt(2*len));                           % NUMBER OF RECORDS (FAST)
    %nr=round(max(roots([1 -1 -2*len])));           % NUMBER OF RECORDS (OLD)
    if((nr^2-nr)/2~=len); error('cv bad size'); end % ASSURE LENGTH IS OK
    if(~isvector(lv)); error('lv not vector'); end  % ASSURE LV IS VECTOR
    if(length(lv)~=len); error('lv bad size'); end  % ASSURE LV LENGTH IS OK
    [j,i]=ind2sub([nr nr],find(tril(true(nr),-1))); % MATRIX SUBSCRIPTS FOR VECTORS
else
    [nr]=size(cv);                                  % GET GRID SIZE
    if(length(nr)>2); error('cg 2D only'); end      % 2D GRIDS ONLY
    if(nr(1)~=nr(2)); error('not square grid'); end % ASSURE GRID IS SQUARE
    if(~isequal(size(lv),nr)); error('lg bad'); end % ASSURE LG MATCHES
    cv=cv.'; lv=lv.';                               % TRANSPOSE TO ACCESS UPPER TRIANGLE
    li=find(triu(true(nr(1)),1));                   % LINEAR INDICES
    cv=cv(li); lv=lv(li);                           % GRID TO VECTOR
    [j,i]=ind2sub(nr,li);                           % MATRIX SUBSCRIPTS FOR VECTORS
    nr=nr(1); len=(nr^2-nr)/2;                      % NUMBER OF CORRELATIONS
end

% ARE WE FINDING RELATIVE ARRIVALS OR ABSOLUTE ARRIVALS?
if(nargin<5)
    % BUILDING KERNEL MATRIX (G)
    totlen=len+1;
    G=sparse([1:len 1:len],[i j],[-ones(1,len) ones(1,len)],totlen,nr,2*len+nr);
    G(totlen,:)=1;
    
    % BUILDING THE WEIGHTING MATRIX (W)
    W=sparse(1:totlen,1:totlen,[xcw(:).*cv(:).^xcpow(:); 1],totlen,totlen);
    
    % GENERALIZED INVERSE (Gg) (OVERDETERMINED CASE)
    Gg=full(inv(G.'*W*G)*G.'*W);
    
    % FINDING LEAST SQUARES RELATIVE ARRIVALS (m)
    m=Gg*[lv(:); 0];
% TIE TO ABSOLUTE ARRIVALS
else
    % BUILDING KERNEL MATRIX (G)
    alen=length(absi); totlen=alen+len;
    G=sparse([1:len 1:len 1:alen],[i j absi(:).'],[-ones(1,len) ones(1,len) ones(1,alen)],totlen,nr,2*len+alen);
    
    % BUILDING THE WEIGHTING MATRIX (W)
    W=sparse(1:totlen,1:totlen,[xcw(:).*cv(:).^xcpow(:); absw(:).*ones(alen,1)],totlen,totlen);
    
    % GENERALIZED INVERSE (Gg) (OVERDETERMINED CASE)
    Gg=full(inv(G.'*W*G)*G.'*W);
    
    % FINDING LEAST SQUARES ABSOLUTE ARRIVALS (m)
    m=Gg*[lv(:); abstt(:)];
end

end
