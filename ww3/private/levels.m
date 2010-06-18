%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% levels                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lev=levels(oct10,oct11,oct12)

k=bitshift2(oct11,oct12);

switch(oct10)
   case 1,lev='surface';
   case 2,lev= 'cloud base';
   case 3,lev= 'cloud top';
   case 4,lev= '0C isotherm';
   case 5,lev= 'cond lev';
   case 6,lev= 'max wind lev';
   case 7,lev= 'tropopause';
   case 8,lev= 'nom. top';
   case 9,lev= 'sea bottom';
   case 10,lev='atmos col';
   case 100, lev=sprintf('%d mb',k);
   case 101, lev=sprintf('%d-%d mb',oct11*10,oct12*10);
   case 102, lev=sprintf('MSL');
   case 103, lev=sprintf('%d m above MSL',k);
   case 104, lev=sprintf('%d-%d m above msl',oct11*100,oct12*100);
   case 105, lev=sprintf('%d m above gnd',k);
   case 106, lev=sprintf('%d-%d m above gnd',oct11*100,oct12*100);
   case 107, lev=sprintf('sigma=%.4f',k/10000.0);
   case 108, lev=sprintf('sigma %.2f-%.2f',oct11/100.0,oct12/100.0);
   case 109, lev=sprintf('hybrid lev %d',k);
   case 110, lev=sprintf('hybrid %d-%d',oct11,oct12);
   case 111, lev=sprintf('%d cm down',k);
   case 112, lev=sprintf('%d-%d cm down',oct11,oct12);
   case 113, lev=sprintf('%dK',k);
   case 114, lev=sprintf('%d-%dK',475-oct11,475-oct12);
   case 115, lev=sprintf('%d mb above gnd',k);
   case 116, lev=sprintf('%d-%d mb above gnd',oct11,oct12);
   case 121, lev=sprintf('%d-%d mb',1100-oct11,1100-oct12);
   case 212, lev=sprintf('low cld bot');
   case 213, lev=sprintf('low cld top');
   case 214, lev=sprintf('low cld lay');
   case 222, lev=sprintf('mid cld bot');
   case 223, lev=sprintf('mid cld top');
   case 224, lev=sprintf('mid cld lay');
   case 232, lev=sprintf('high cld bot');
   case 233, lev=sprintf('high cld top');
   case 234, lev=sprintf('high cld lay');
   otherwise,lev='unknown';   
end
