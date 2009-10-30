function [anis]=convert_anis(anis)

sumsq=anis(:,3).^2+anis(:,5).^2;
vanis=100*sqrt(sumsq)./anis(:,2);
stdanis=100*sqrt(((anis(:,3).*anis(:,4)).^2+(anis(:,5).*anis(:,6)).^2)./sumsq)./anis(:,2);
azim=0.5*180/pi*atan2(anis(:,5),anis(:,3));
stdazim=0.5*180/pi*sqrt((anis(:,3).*anis(:,6)).^2+(anis(:,5).*anis(:,4)).^2)./sumsq;
anis(:,3:6)=[vanis stdanis azim stdazim];

end