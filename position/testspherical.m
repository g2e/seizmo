function []=testspherical(n,tol)

% random start/end points
lat1=rand(n,1)*180-90;
lon1=rand(n,1)*360-180;
lat2=rand(n,1)*180-90;
lon2=rand(n,1)*360-180;

% get relative positioning (forward and backward)
[dist1,az1,baz1]=sphericalinv(lat1,lon1,lat2,lon2);
[dist2,az2,baz2]=sphericalinv(lat2,lon2,lat1,lon1);

% murder antipode pairs
kill=dist1>20000;
dist1(kill)=[]; dist2(kill)=[];
az1(kill)=[]; az2(kill)=[];
baz1(kill)=[]; baz2(kill)=[];
lat1(kill)=[]; lat2(kill)=[];
lon1(kill)=[]; lon2(kill)=[];

% fix azimuths to always be 0<=az<360
az1=mod(az1,360); az2=mod(az2,360);
baz1=mod(baz1,360); baz2=mod(baz2,360);

% check that these match (within tolerance)
if(any(abs(dist1-dist2)>tol))
    disp(max(abs(dist1-dist2)))
    error('distances dont match')
end
if(any(abs(az1-baz2)>tol) || any(abs(baz1-az2)>tol))
    disp(max(abs(az1-baz2)))
    disp(max(abs(baz1-az2)))
    error('azimuths dont match')
end

% get new positioning using relative positioning (forward and backward)
[lat3,lon3,baz3]=sphericalfwd(lat1,lon1,dist1,az1);
[lat4,lon4,baz4]=sphericalfwd(lat2,lon2,dist1,baz1);
[lat5,lon5,baz5]=sphericalfwd(lat1,lon1,dist2,baz2);
[lat6,lon6,baz6]=sphericalfwd(lat2,lon2,dist2,az2);

% fix azimuths to always be 0<=az<360
baz3=mod(baz3,360); baz4=mod(baz4,360);
baz5=mod(baz5,360); baz6=mod(baz6,360);

% cross check
if(any(abs(lat3-lat2)>tol) || any(abs(lon3-lon2)>tol) ||...
        any(abs(lat5-lat2)>tol) || any(abs(lon5-lon2)>tol))
    [blah,i]=max(abs(lat3-lat2));
    [blah,j]=max(abs(lon3-lon2));
    disp([max(abs(lat3-lat2)),dist1(i)])
    disp([max(abs(lon3-lon2)),dist1(j)])
    disp(max(abs(lat5-lat2)))
    disp(max(abs(lon5-lon2)))
    error('crosscheck 1')
end
if(any(abs(lat4-lat1)>tol) || any(abs(lon4-lon1)>tol) ||...
        any(abs(lat6-lat1)>tol) || any(abs(lon6-lon1)>tol))
    disp(max(abs(lat4-lat1)))
    disp(max(abs(lon4-lon1)))
    disp(max(abs(lat6-lat1)))
    disp(max(abs(lon6-lon1)))
    error('crosscheck 2')
end
if(any(abs(baz3-baz1)>tol) || any(abs(baz6-baz2)>tol) ||...
        any(abs(baz4-az1)>tol) || any(abs(baz5-az2)>tol))
    disp(max(abs(baz3-baz1)))
    disp(max(abs(baz6-baz2)))
    disp(max(abs(baz4-az1)))
    disp(max(abs(baz5-az2)))
    error('az crosscheck')
end

end
