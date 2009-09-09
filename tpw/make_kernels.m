
% need freq
% totwin -2000 7000
% delta 4
% window from avg - start at 2000
% tprlen from avg - 0.2 * avg win
% need vel
% lambda is 110 - 2 * node spacing

% this must be changed
bank=flipud(filter_bank([0.0055 0.055],'variable',0.2,0.1));

% this must be changed
wlen=[ 463.14
        450.4
       479.82
       548.04
       524.07
       547.09
       552.15
       557.55
       581.55
       617.31
        644.3
       713.39
       754.56
       823.97
       920.38
       993.66
       1111.3
         1200
       1317.7
       1440.4
       1570.7
       1732.2
       1899.4
         2089
       2340.9];

% this must be changed
local=[ % simple fit to 1D results
 18.5000    3.6000
 35.0000    3.8800
170.0000    4.3500];
movein=interp1(local(:,1),local(:,2),1./bank(:,1),[],'extrap');

% this must be changed
for i=1:25
    [f,amp]=getpeak(bank(i,1),[-2000 7000],4,[2000 2000+wlen(i)],0.2*wlen(i));
    [x,y,p,a]=sensitivity(f,amp,movein(i));
    [sp,sa]=smooth2d(x,y,p,a,110);
    writekernel(['sens' sprintf('%03d',round(1/bank(i,1))) 's110km.dat'],x,y,sp,sa);
end

