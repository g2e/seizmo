function [varargout]=bathy_micro_excite(ph,z,f,vs)
%BATHY_MICRO_EXCITE    Bathymetry 2ndary microseism excitation coefficients
%
%    Usage:    c=bathy_micro_excite(ph,z,f,vs)
%              [c0,...,c6]=bathy_micro_excite('rayleigh',z,f,vs)
%
%    Description:
%     C=BATHY_MICRO_EXCITE(PH,Z,F,Vs) returns the bathymetry secondary 
%     microseism excitation coefficient for the specified phase as
%     determined in Ardhuin & Herbers (2013).  These coefficients modulate
%     the bottom pressure field from wave-wave interaction of ocean gravity
%     waves to give the microseism energy (Longuet-Higgins 1950; Hasselmann
%     1963).  PH is the phase name and may be one of the following:
%        'RAYLEIGH', 'P', 'S', 'S-ALL'.
%     Z is the bathymetry and may be an array (C is the same size as Z).  F
%     is the frequency (not the angular frequency) and Vs is the shear
%     velocity of the sea bed in meters per second.  The default values are
%     1/(2pi) Hz (about 6.3s period) for F and 3000 m/s for Vs.
%
%     [C0,...,C6]=BATHY_MICRO_EXCITE('RAYLEIGH',Z,F,VS) gives individual
%     Rayleigh excitation components from C0 to C6.  These can be used to
%     calculate the former usage output C^2=C0^2+C1^2+C2^2+C3^2+...
%
%    Notes:
%     - References:
%        Longuet-Higgins 1950, Theory of the Origin of Microseisms
%        Hasselmann 1963, A Statistical Analysis of the Generation of
%         Microseisms
%        Ardhuin & Herbers 2013, Noise Generation in the Solid Earth,
%         Oceans and Atmosphere, from Nonlinear Interacting Surface Gravity
%         Waves in Finite Depth
%
%    Examples:
%     % The effect of frequency on the Rayleigh
%     % wave excitation coefficient with depth:
%     figure;
%     plot(1:6500,[bathy_micro_excite('R',1:6500,1/10);
%                  bathy_micro_excite('R',1:6500,1/7.5);
%                  bathy_micro_excite('R',1:6500,1/5);
%                  bathy_micro_excite('R',1:6500,1/2.5)],'linewidth',2);
%     set(gca,'fontweight','bold');
%     text(800,.85,'2.5s','fontweight','bold');
%     text(1800,.85,'5s','fontweight','bold');
%     text(2700,.85,'7.5s','fontweight','bold');
%     text(3650,.85,'10s','fontweight','bold');
%     xlabel('Bathymetric Depth (m)','fontweight','bold');
%     ylabel('Rayleigh Wave Bathymetric Excitation Coefficient',...
%         'fontweight','bold');
%     title(['Effect of Wave Period on the Rayleigh Wave ' ...
%         'Bathymetric Excitation Coefficient versus Depth'],...
%         'fontweight','bold');
%     axis tight;
%     ylim([0 1]);
%
%     % The effect of Crustal Vs:
%     figure;
%     plot(1:6500,[bathy_micro_excite('R',1:6500,1/7.5,2000);
%                  bathy_micro_excite('R',1:6500,1/7.5,2400);
%                  bathy_micro_excite('R',1:6500,1/7.5,2800);
%                  bathy_micro_excite('R',1:6500,1/7.5,3200);
%                  bathy_micro_excite('R',1:6500,1/7.5,3600)],...
%         'linewidth',2);
%     set(gca,'fontweight','bold');
%     text(2000,.85,'2.0km/s','fontweight','bold','rotation',45);
%     text(2400,.85,'2.4km/s','fontweight','bold','rotation',45);
%     text(2800,.85,'2.8km/s','fontweight','bold','rotation',45);
%     text(3200,.85,'3.2km/s','fontweight','bold','rotation',45);
%     text(3600,.85,'3.6km/s','fontweight','bold','rotation',45);
%     xlabel('Bathymetric Depth (m)','fontweight','bold');
%     ylabel('Rayleigh Bathymetric Excitation Coefficient',...
%         'fontweight','bold');
%     title(['Effect of Crustal Shear Velocity ' ...
%         'on Rayleigh Bathymetric Excitation Coefficient'],...
%         'fontweight','bold');
%     axis tight;
%     ylim([0 1]);
%
%     % Rayleigh Bathymetric Excitation Map:
%     [lon,lat]=meshgrid(-179.5:179.5,-89.5:89.5);
%     mod=getcrust(lat,lon);
%     mod.top(mod.top(:,2)>0,2)=0; % mask out land
%     mod.vs(mod.vs==0)=1; % avoid divide by zero
%     vsavg=sum(mod.thk(:,2:8),2)./sum(mod.thk(:,2:8)./mod.vs(:,2:8),2);
%     c=bathy_micro_excite('R',-mod.top(:,2)*1000,1/7.5,vsavg*1000);
%     ax=plotbathyexcite('R',reshape(c,[180 360]),lat,lon);
%     title(ax,{[] ...
%         'Rayleigh Bathymetric Excitation Map for' ...
%         'Wave-Wave Interference of Ocean Waves' ...
%         'Seismic Period: 7.5s' []},'color','w');
%
%    See also: PLOTBATHYEXCITE, FS_PHASE2LATLON, SLOWNESS2DEG, DEG2SLOWNESS

%     Version History:
%        Sep. 22, 2010 - initial version
%        May   5, 2012 - added couple more examples and some doc fixes
%        May  18, 2012 - minor touches
%        Jan. 14, 2014 - updated for Ardhuin & Herbers and forcing phase
%                        choice because it matters, fixed example
%        Jan. 23, 2014 - update to example for recent changes
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 23, 2014 at 10:35 GMT

% todo

% check nargin
error(nargchk(2,4,nargin));

% defaults
if(nargin<3 || isempty(f)); f=1/(2*pi); end
if(nargin<4 || isempty(vs)); vs=3000; end

% check inputs
if(~ischar(ph) || ndims(ph)>2 || size(ph,1)~=1)
    error('seizmo:bathy_micro_excite:badInput',...
        'PH must be a string!');
elseif(~isreal(z) || any(z(:)<0))
    error('seizmo:bathy_micro_excite:badInput',...
        'Z must be the bathymetry (positive reals in meters)!');
elseif(~isreal(f) || any(f(:)<0))
    error('seizmo:bathy_micro_excite:badInput',...
        'F must be the frequency (positive reals in Hz)!');
elseif(~isreal(vs) || any(vs(:)<100))
    error('seizmo:bathy_micro_excite:badInput',...
        'Vs must be the sea bed shear velocity (positive reals in m/s)!');
end

% normalized
x=f.*z./vs;

% bathymetry modulated excitation coefficients
% Column 1 = freq * bathymetry / shear velocity of sea floor
% Column 2+= coefficients of excitation
%        x        cR        c0        c1        c2        c3        c4        c5        c6        cP        cS     cSall
xc=[0.0000         0       NaN       NaN       NaN       NaN       NaN       NaN       NaN    0.0000    0.0000    0.0000
    0.0021    0.4388    0.4388       NaN       NaN       NaN       NaN       NaN       NaN    0.2906    0.0784    0.4319
    0.0042    0.4403    0.4403       NaN       NaN       NaN       NaN       NaN       NaN    0.2917    0.0774    0.4674
    0.0064    0.4419    0.4419       NaN       NaN       NaN       NaN       NaN       NaN    0.2923    0.0770    0.4797
    0.0085    0.4436    0.4436       NaN       NaN       NaN       NaN       NaN       NaN    0.2930    0.0769    0.4863
    0.0106    0.4455    0.4455       NaN       NaN       NaN       NaN       NaN       NaN    0.2939    0.0769    0.4906
    0.0127    0.4475    0.4475       NaN       NaN       NaN       NaN       NaN       NaN    0.2951    0.0771    0.4941
    0.0149    0.4496    0.4496       NaN       NaN       NaN       NaN       NaN       NaN    0.2964    0.0773    0.4971
    0.0170    0.4518    0.4518       NaN       NaN       NaN       NaN       NaN       NaN    0.2980    0.0776    0.5000
    0.0191    0.4542    0.4542       NaN       NaN       NaN       NaN       NaN       NaN    0.2999    0.0780    0.5030
    0.0212    0.4567    0.4567       NaN       NaN       NaN       NaN       NaN       NaN    0.3019    0.0785    0.5061
    0.0234    0.4594    0.4594       NaN       NaN       NaN       NaN       NaN       NaN    0.3042    0.0790    0.5094
    0.0255    0.4623    0.4623       NaN       NaN       NaN       NaN       NaN       NaN    0.3068    0.0796    0.5129
    0.0276    0.4653    0.4653       NaN       NaN       NaN       NaN       NaN       NaN    0.3097    0.0803    0.5166
    0.0297    0.4684    0.4684       NaN       NaN       NaN       NaN       NaN       NaN    0.3128    0.0811    0.5207
    0.0318    0.4717    0.4717       NaN       NaN       NaN       NaN       NaN       NaN    0.3163    0.0821    0.5251
    0.0340    0.4753    0.4753       NaN       NaN       NaN       NaN       NaN       NaN    0.3200    0.0831    0.5299
    0.0361    0.4789    0.4789       NaN       NaN       NaN       NaN       NaN       NaN    0.3241    0.0841    0.5351
    0.0382    0.4828    0.4828       NaN       NaN       NaN       NaN       NaN       NaN    0.3285    0.0852    0.5406
    0.0403    0.4869    0.4869       NaN       NaN       NaN       NaN       NaN       NaN    0.3333    0.0864    0.5465
    0.0425    0.4912    0.4912       NaN       NaN       NaN       NaN       NaN       NaN    0.3384    0.0877    0.5529
    0.0446    0.4957    0.4957       NaN       NaN       NaN       NaN       NaN       NaN    0.3440    0.0891    0.5597
    0.0467    0.5005    0.5005       NaN       NaN       NaN       NaN       NaN       NaN    0.3500    0.0906    0.5670
    0.0488    0.5054    0.5054       NaN       NaN       NaN       NaN       NaN       NaN    0.3565    0.0922    0.5749
    0.0509    0.5107    0.5107       NaN       NaN       NaN       NaN       NaN       NaN    0.3634    0.0939    0.5833
    0.0531    0.5162    0.5162       NaN       NaN       NaN       NaN       NaN       NaN    0.3709    0.0958    0.5924
    0.0552    0.5219    0.5219       NaN       NaN       NaN       NaN       NaN       NaN    0.3791    0.0978    0.6020
    0.0573    0.5280    0.5280       NaN       NaN       NaN       NaN       NaN       NaN    0.3878    0.0999    0.6124
    0.0594    0.5344    0.5344       NaN       NaN       NaN       NaN       NaN       NaN    0.3973    0.1023    0.6235
    0.0616    0.5411    0.5411       NaN       NaN       NaN       NaN       NaN       NaN    0.4075    0.1048    0.6353
    0.0637    0.5481    0.5481       NaN       NaN       NaN       NaN       NaN       NaN    0.4186    0.1076    0.6481
    0.0658    0.5555    0.5555       NaN       NaN       NaN       NaN       NaN       NaN    0.4306    0.1106    0.6618
    0.0679    0.5633    0.5633       NaN       NaN       NaN       NaN       NaN       NaN    0.4436    0.1138    0.6764
    0.0701    0.5715    0.5715       NaN       NaN       NaN       NaN       NaN       NaN    0.4577    0.1173    0.6922
    0.0722    0.5802    0.5802       NaN       NaN       NaN       NaN       NaN       NaN    0.4731    0.1210    0.7091
    0.0743    0.5892    0.5892       NaN       NaN       NaN       NaN       NaN       NaN    0.4900    0.1251    0.7274
    0.0764    0.5988    0.5988       NaN       NaN       NaN       NaN       NaN       NaN    0.5084    0.1296    0.7471
    0.0785    0.6089    0.6089       NaN       NaN       NaN       NaN       NaN       NaN    0.5287    0.1345    0.7683
    0.0807    0.6195    0.6195       NaN       NaN       NaN       NaN       NaN       NaN    0.5510    0.1399    0.7913
    0.0828    0.6307    0.6307       NaN       NaN       NaN       NaN       NaN       NaN    0.5757    0.1459    0.8163
    0.0849    0.6426    0.6426       NaN       NaN       NaN       NaN       NaN       NaN    0.6030    0.1524    0.8434
    0.0870    0.6550    0.6550       NaN       NaN       NaN       NaN       NaN       NaN    0.6336    0.1597    0.8729
    0.0892    0.6682    0.6682       NaN       NaN       NaN       NaN       NaN       NaN    0.6678    0.1679    0.9052
    0.0913    0.6820    0.6820       NaN       NaN       NaN       NaN       NaN       NaN    0.7063    0.1770    0.9405
    0.0934    0.6966    0.6966       NaN       NaN       NaN       NaN       NaN       NaN    0.7500    0.1873    0.9794
    0.0955    0.7120    0.7120       NaN       NaN       NaN       NaN       NaN       NaN    0.7999    0.1990    1.0223
    0.0977    0.7282    0.7282       NaN       NaN       NaN       NaN       NaN       NaN    0.8573    0.2123    1.0699
    0.0998    0.7452    0.7452       NaN       NaN       NaN       NaN       NaN       NaN    0.9238    0.2276    1.1228
    0.1019    0.7630    0.7630       NaN       NaN       NaN       NaN       NaN       NaN    1.0017    0.2454    1.1820
    0.1040    0.7816    0.7816       NaN       NaN       NaN       NaN       NaN       NaN    1.0938    0.2662    1.2487
    0.1061    0.8010    0.8010       NaN       NaN       NaN       NaN       NaN       NaN    1.2041    0.2908    1.3242
    0.1083    0.8211    0.8211       NaN       NaN       NaN       NaN       NaN       NaN    1.3376    0.3202    1.4103
    0.1104    0.8417    0.8417       NaN       NaN       NaN       NaN       NaN       NaN    1.5013    0.3558    1.5092
    0.1125    0.8627    0.8627       NaN       NaN       NaN       NaN       NaN       NaN    1.7039    0.3995    1.6237
    0.1146    0.8840    0.8840       NaN       NaN       NaN       NaN       NaN       NaN    1.9551    0.4536    1.7529
    0.1168    0.9050    0.9050       NaN       NaN       NaN       NaN       NaN       NaN    2.2613    0.5206    1.9097
    0.1189    0.9257    0.9257       NaN       NaN       NaN       NaN       NaN       NaN    2.6118    0.6023    2.0943
    0.1210    0.9453    0.9453       NaN       NaN       NaN       NaN       NaN       NaN    2.9513    0.6965    2.3099
    0.1231    0.9636    0.9636       NaN       NaN       NaN       NaN       NaN       NaN    3.1679    0.7911    2.5548
    0.1252    0.9799    0.9799       NaN       NaN       NaN       NaN       NaN       NaN    3.1637    0.8601    2.8194
    0.1274    0.9939    0.9939       NaN       NaN       NaN       NaN       NaN       NaN    2.9498    0.8710    3.0865
    0.1295    1.0050    1.0050       NaN       NaN       NaN       NaN       NaN       NaN    2.6192    0.8050    3.3415
    0.1316    1.0129    1.0129       NaN       NaN       NaN       NaN       NaN       NaN    2.2730    0.7026    3.5727
    0.1337    1.0174    1.0174       NaN       NaN       NaN       NaN       NaN       NaN    1.9671    0.6028    3.7512
    0.1359    1.0186    1.0186       NaN       NaN       NaN       NaN       NaN       NaN    1.7145    0.5186    3.8912
    0.1380    1.0165    1.0165       NaN       NaN       NaN       NaN       NaN       NaN    1.5102    0.4509    4.0026
    0.1401    1.0114    1.0114       NaN       NaN       NaN       NaN       NaN       NaN    1.3450    0.3968    4.0960
    0.1422    1.0037    1.0037       NaN       NaN       NaN       NaN       NaN       NaN    1.2102    0.3533    4.1831
    0.1444    0.9937    0.9937       NaN       NaN       NaN       NaN       NaN       NaN    1.0989    0.3179    4.2775
    0.1465    0.9820    0.9820       NaN       NaN       NaN       NaN       NaN       NaN    1.0060    0.2887    4.3985
    0.1486    0.9688    0.9688       NaN       NaN       NaN       NaN       NaN       NaN    0.9275    0.2644    4.5809
    0.1507    0.9545    0.9545       NaN       NaN       NaN       NaN       NaN       NaN    0.8605    0.2438    4.9071
    0.1528    0.9395    0.9395       NaN       NaN       NaN       NaN       NaN       NaN    0.8027    0.2263    5.3960
    0.1550    0.9241    0.9241       NaN       NaN       NaN       NaN       NaN       NaN    0.7525    0.2111    4.2535
    0.1571    0.9125    0.9084    0.3340       NaN       NaN       NaN       NaN       NaN    0.7085    0.1980    3.3316
    0.1592    0.9000    0.8926    0.3819       NaN       NaN       NaN       NaN       NaN    0.6697    0.1864    2.7811
    0.1613    0.8874    0.8769    0.4124       NaN       NaN       NaN       NaN       NaN    0.6353    0.1762    2.4062
    0.1635    0.8747    0.8613    0.4329       NaN       NaN       NaN       NaN       NaN    0.6046    0.1672    2.1303
    0.1656    0.8620    0.8460    0.4470       NaN       NaN       NaN       NaN       NaN    0.5770    0.1591    1.9170
    0.1677    0.8494    0.8310    0.4567       NaN       NaN       NaN       NaN       NaN    0.5522    0.1519    1.7464
    0.1698    0.8368    0.8164    0.4634       NaN       NaN       NaN       NaN       NaN    0.5298    0.1454    1.6064
    0.1719    0.8244    0.8021    0.4678       NaN       NaN       NaN       NaN       NaN    0.5094    0.1394    1.4892
    0.1741    0.8122    0.7882    0.4706       NaN       NaN       NaN       NaN       NaN    0.4909    0.1341    1.3896
    0.1762    0.8002    0.7748    0.4721       NaN       NaN       NaN       NaN       NaN    0.4740    0.1292    1.3038
    0.1783    0.7885    0.7617    0.4727       NaN       NaN       NaN       NaN       NaN    0.4585    0.1247    1.2291
    0.1804    0.7770    0.7490    0.4725       NaN       NaN       NaN       NaN       NaN    0.4443    0.1207    1.1635
    0.1826    0.7659    0.7367    0.4718       NaN       NaN       NaN       NaN       NaN    0.4312    0.1169    1.1054
    0.1847    0.7551    0.7248    0.4706       NaN       NaN       NaN       NaN       NaN    0.4192    0.1135    1.0537
    0.1868    0.7445    0.7133    0.4692       NaN       NaN       NaN       NaN       NaN    0.4081    0.1103    1.0073
    0.1889    0.7343    0.7021    0.4675       NaN       NaN       NaN       NaN       NaN    0.3978    0.1074    0.9655
    0.1911    0.7244    0.6913    0.4656       NaN       NaN       NaN       NaN       NaN    0.3883    0.1046    0.9277
    0.1932    0.7148    0.6808    0.4636       NaN       NaN       NaN       NaN       NaN    0.3795    0.1021    0.8933
    0.1953    0.7054    0.6706    0.4616       NaN       NaN       NaN       NaN       NaN    0.3714    0.0998    0.8620
    0.1974    0.6964    0.6608    0.4594       NaN       NaN       NaN       NaN       NaN    0.3638    0.0976    0.8333
    0.1995    0.6876    0.6512    0.4573       NaN       NaN       NaN       NaN       NaN    0.3568    0.0956    0.8070
    0.2017    0.6792    0.6420    0.4551       NaN       NaN       NaN       NaN       NaN    0.3503    0.0938    0.7828
    0.2038    0.6709    0.6330    0.4530       NaN       NaN       NaN       NaN       NaN    0.3443    0.0920    0.7605
    0.2059    0.6630    0.6243    0.4509       NaN       NaN       NaN       NaN       NaN    0.3387    0.0904    0.7399
    0.2080    0.6553    0.6158    0.4489       NaN       NaN       NaN       NaN       NaN    0.3336    0.0889    0.7209
    0.2102    0.6478    0.6076    0.4469       NaN       NaN       NaN       NaN       NaN    0.3288    0.0876    0.7032
    0.2123    0.6406    0.5996    0.4449       NaN       NaN       NaN       NaN       NaN    0.3244    0.0863    0.6868
    0.2144    0.6336    0.5918    0.4430       NaN       NaN       NaN       NaN       NaN    0.3203    0.0851    0.6716
    0.2165    0.6269    0.5842    0.4412       NaN       NaN       NaN       NaN       NaN    0.3165    0.0840    0.6575
    0.2187    0.6203    0.5769    0.4395       NaN       NaN       NaN       NaN       NaN    0.3131    0.0830    0.6443
    0.2208    0.6140    0.5697    0.4378       NaN       NaN       NaN       NaN       NaN    0.3099    0.0821    0.6320
    0.2229    0.6079    0.5628    0.4362       NaN       NaN       NaN       NaN       NaN    0.3071    0.0812    0.6206
    0.2250    0.6019    0.5560    0.4347       NaN       NaN       NaN       NaN       NaN    0.3045    0.0804    0.6100
    0.2271    0.5962    0.5494    0.4332       NaN       NaN       NaN       NaN       NaN    0.3021    0.0797    0.5979
    0.2293    0.5907    0.5430    0.4319       NaN       NaN       NaN       NaN       NaN    0.3000    0.0791    0.5888
    0.2314    0.5853    0.5367    0.4306       NaN       NaN       NaN       NaN       NaN    0.2982    0.0785    0.5802
    0.2335    0.5801    0.5306    0.4295       NaN       NaN       NaN       NaN       NaN    0.2966    0.0780    0.5723
    0.2356    0.5751    0.5246    0.4284       NaN       NaN       NaN       NaN       NaN    0.2952    0.0775    0.5649
    0.2378    0.5703    0.5188    0.4273       NaN       NaN       NaN       NaN       NaN    0.2940    0.0771    0.5581
    0.2399    0.5657    0.5131    0.4264       NaN       NaN       NaN       NaN       NaN    0.2930    0.0768    0.5518
    0.2420    0.5612    0.5076    0.4256       NaN       NaN       NaN       NaN       NaN    0.2923    0.0765    0.5459
    0.2441    0.5568    0.5021    0.4248       NaN       NaN       NaN       NaN       NaN    0.2918    0.0762    0.5405
    0.2462    0.5527    0.4968    0.4241       NaN       NaN       NaN       NaN       NaN    0.2915    0.0761    0.5356
    0.2484    0.5487    0.4917    0.4236       NaN       NaN       NaN       NaN       NaN    0.2914    0.0759    0.5310
    0.2505    0.5448    0.4866    0.4231       NaN       NaN       NaN       NaN       NaN    0.2915    0.0759    0.5269
    0.2526    0.5411    0.4816    0.4226       NaN       NaN       NaN       NaN       NaN    0.2918    0.0758    0.5231
    0.2547    0.5375    0.4768    0.4223       NaN       NaN       NaN       NaN       NaN    0.2923    0.0759    0.5197
    0.2569    0.5341    0.4721    0.4221       NaN       NaN       NaN       NaN       NaN    0.2930    0.0760    0.5167
    0.2590    0.5309    0.4674    0.4219       NaN       NaN       NaN       NaN       NaN    0.2940    0.0761    0.5140
    0.2611    0.5277    0.4629    0.4218       NaN       NaN       NaN       NaN       NaN    0.2951    0.0763    0.5117
    0.2632    0.5248    0.4584    0.4219       NaN       NaN       NaN       NaN       NaN    0.2965    0.0765    0.5097
    0.2654    0.5219    0.4541    0.4220       NaN       NaN       NaN       NaN       NaN    0.2982    0.0768    0.5081
    0.2675    0.5192    0.4498    0.4221       NaN       NaN       NaN       NaN       NaN    0.3000    0.0772    0.5068
    0.2696    0.5167    0.4456    0.4224       NaN       NaN       NaN       NaN       NaN    0.3021    0.0776    0.5057
    0.2717    0.5143    0.4415    0.4228       NaN       NaN       NaN       NaN       NaN    0.3045    0.0781    0.5051
    0.2738    0.5120    0.4375    0.4232       NaN       NaN       NaN       NaN       NaN    0.3071    0.0786    0.5047
    0.2760    0.5099    0.4336    0.4238       NaN       NaN       NaN       NaN       NaN    0.3099    0.0792    0.5047
    0.2781    0.5079    0.4297    0.4244       NaN       NaN       NaN       NaN       NaN    0.3131    0.0799    0.5049
    0.2802    0.5060    0.4259    0.4251       NaN       NaN       NaN       NaN       NaN    0.3165    0.0806    0.5055
    0.2823    0.5043    0.4222    0.4259       NaN       NaN       NaN       NaN       NaN    0.3203    0.0814    0.5064
    0.2845    0.5027    0.4185    0.4268       NaN       NaN       NaN       NaN       NaN    0.3244    0.0823    0.5077
    0.2866    0.5013    0.4149    0.4278       NaN       NaN       NaN       NaN       NaN    0.3288    0.0832    0.5092
    0.2887    0.5000    0.4114    0.4289       NaN       NaN       NaN       NaN       NaN    0.3336    0.0843    0.5111
    0.2908    0.4988    0.4080    0.4301       NaN       NaN       NaN       NaN       NaN    0.3388    0.0854    0.5134
    0.2930    0.4978    0.4046    0.4314       NaN       NaN       NaN       NaN       NaN    0.3444    0.0866    0.5160
    0.2951    0.4970    0.4012    0.4328       NaN       NaN       NaN       NaN       NaN    0.3504    0.0879    0.5189
    0.2972    0.4963    0.3979    0.4343       NaN       NaN       NaN       NaN       NaN    0.3570    0.0893    0.5222
    0.2993    0.4957    0.3947    0.4359       NaN       NaN       NaN       NaN       NaN    0.3640    0.0908    0.5259
    0.3014    0.4953    0.3915    0.4376       NaN       NaN       NaN       NaN       NaN    0.3716    0.0925    0.5301
    0.3036    0.4950    0.3884    0.4394       NaN       NaN       NaN       NaN       NaN    0.3797    0.0942    0.5346
    0.3057    0.4949    0.3854    0.4413       NaN       NaN       NaN       NaN       NaN    0.3886    0.0961    0.5395
    0.3078    0.4949    0.3824    0.4433       NaN       NaN       NaN       NaN       NaN    0.3981    0.0982    0.5449
    0.3099    0.4951    0.3794    0.4455       NaN       NaN       NaN       NaN       NaN    0.4085    0.1004    0.5508
    0.3121    0.4955    0.3765    0.4477       NaN       NaN       NaN       NaN       NaN    0.4197    0.1028    0.5572
    0.3142    0.4960    0.3736    0.4501       NaN       NaN       NaN       NaN       NaN    0.4318    0.1053    0.5641
    0.3163    0.4967    0.3708    0.4526       NaN       NaN       NaN       NaN       NaN    0.4450    0.1081    0.5716
    0.3184    0.4975    0.3680    0.4552       NaN       NaN       NaN       NaN       NaN    0.4593    0.1111    0.5797
    0.3205    0.4986    0.3653    0.4579       NaN       NaN       NaN       NaN       NaN    0.4750    0.1144    0.5884
    0.3227    0.4998    0.3626    0.4608       NaN       NaN       NaN       NaN       NaN    0.4921    0.1179    0.5978
    0.3248    0.5011    0.3600    0.4638       NaN       NaN       NaN       NaN       NaN    0.5109    0.1218    0.6080
    0.3269    0.5027    0.3573    0.4670       NaN       NaN       NaN       NaN       NaN    0.5315    0.1260    0.6189
    0.3290    0.5044    0.3548    0.4703       NaN       NaN       NaN       NaN       NaN    0.5543    0.1306    0.6307
    0.3312    0.5064    0.3522    0.4737       NaN       NaN       NaN       NaN       NaN    0.5796    0.1356    0.6435
    0.3333    0.5085    0.3498    0.4773       NaN       NaN       NaN       NaN       NaN    0.6076    0.1412    0.6572
    0.3354    0.5108    0.3473    0.4810       NaN       NaN       NaN       NaN       NaN    0.6391    0.1473    0.6721
    0.3375    0.5133    0.3449    0.4849       NaN       NaN       NaN       NaN       NaN    0.6744    0.1541    0.6883
    0.3397    0.5161    0.3425    0.4890       NaN       NaN       NaN       NaN       NaN    0.7143    0.1616    0.7058
    0.3418    0.5190    0.3401    0.4932       NaN       NaN       NaN       NaN       NaN    0.7597    0.1701    0.7239
    0.3439    0.5221    0.3378    0.4976       NaN       NaN       NaN       NaN       NaN    0.8118    0.1796    0.7447
    0.3460    0.5255    0.3355    0.5022       NaN       NaN       NaN       NaN       NaN    0.8721    0.1903    0.7674
    0.3481    0.5291    0.3333    0.5069       NaN       NaN       NaN       NaN       NaN    0.9425    0.2025    0.7922
    0.3503    0.5329    0.3311    0.5118       NaN       NaN       NaN       NaN       NaN    1.0255    0.2165    0.8196
    0.3524    0.5369    0.3289    0.5169       NaN       NaN       NaN       NaN       NaN    1.1246    0.2327    0.8499
    0.3545    0.5412    0.3267    0.5222       NaN       NaN       NaN       NaN       NaN    1.2443    0.2517    0.8836
    0.3566    0.5457    0.3246    0.5277       NaN       NaN       NaN       NaN       NaN    1.3908    0.2740    0.9214
    0.3588    0.5504    0.3225    0.5334       NaN       NaN       NaN       NaN       NaN    1.5716    0.3007    0.9639
    0.3609    0.5554    0.3204    0.5393       NaN       NaN       NaN       NaN       NaN    1.7951    0.3329    1.0123
    0.3630    0.5606    0.3184    0.5454       NaN       NaN       NaN       NaN       NaN    2.0655    0.3722    1.0679
    0.3651    0.5660    0.3164    0.5517       NaN       NaN       NaN       NaN       NaN    2.3704    0.4200    1.1322
    0.3672    0.5717    0.3144    0.5582       NaN       NaN       NaN       NaN       NaN    2.6626    0.4777    1.2071
    0.3694    0.5776    0.3124    0.5649       NaN       NaN       NaN       NaN       NaN    2.8682    0.5447    1.2939
    0.3715    0.5838    0.3105    0.5717       NaN       NaN       NaN       NaN       NaN    2.9412    0.6175    1.3926
    0.3736    0.5901    0.3085    0.5788       NaN       NaN       NaN       NaN       NaN    2.8908    0.6899    1.5019
    0.3757    0.5967    0.3067    0.5860       NaN       NaN       NaN       NaN       NaN    2.7525    0.7537    1.6197
    0.3779    0.6035    0.3048    0.5935       NaN       NaN       NaN       NaN       NaN    2.5598    0.7997    1.7434
    0.3800    0.6105    0.3029    0.6010       NaN       NaN       NaN       NaN       NaN    2.3366    0.8178    1.8710
    0.3821    0.6176    0.3011    0.6087       NaN       NaN       NaN       NaN       NaN    2.0987    0.7972    2.0010
    0.3842    0.6249    0.2993    0.6166       NaN       NaN       NaN       NaN       NaN    1.8575    0.7281    2.1350
    0.3864    0.6324    0.2975    0.6245       NaN       NaN       NaN       NaN       NaN    1.6318    0.6309    2.2869
    0.3885    0.6399    0.2958    0.6325       NaN       NaN       NaN       NaN       NaN    1.4398    0.5406    2.4115
    0.3906    0.6475    0.2941    0.6405       NaN       NaN       NaN       NaN       NaN    1.2825    0.4666    2.5204
    0.3927    0.6550    0.2923    0.6484       NaN       NaN       NaN       NaN       NaN    1.1543    0.4078    2.6180
    0.3948    0.6626    0.2907    0.6564       NaN       NaN       NaN       NaN       NaN    1.0489    0.3611    2.7073
    0.3970    0.6700    0.2890    0.6642       NaN       NaN       NaN       NaN       NaN    0.9612    0.3235    2.7901
    0.3991    0.6773    0.2873    0.6718       NaN       NaN       NaN       NaN       NaN    0.8873    0.2928    2.8676
    0.4012    0.6845    0.2857    0.6792       NaN       NaN       NaN       NaN       NaN    0.8244    0.2674    2.9406
    0.4033    0.6913    0.2841    0.6863       NaN       NaN       NaN       NaN       NaN    0.7702    0.2461    3.0099
    0.4055    0.6978    0.2825    0.6930       NaN       NaN       NaN       NaN       NaN    0.7232    0.2280    3.0761
    0.4076    0.7039    0.2809    0.6994       NaN       NaN       NaN       NaN       NaN    0.6820    0.2124    3.1397
    0.4097    0.7095    0.2794    0.7052       NaN       NaN       NaN       NaN       NaN    0.6457    0.1989    3.2017
    0.4118    0.7145    0.2778    0.7104       NaN       NaN       NaN       NaN       NaN    0.6135    0.1872    3.2628
    0.4140    0.7190    0.2763    0.7151       NaN       NaN       NaN       NaN       NaN    0.5847    0.1768    3.3243
    0.4161    0.7228    0.2748    0.7190       NaN       NaN       NaN       NaN       NaN    0.5588    0.1677    3.3878
    0.4182    0.7260    0.2733    0.7223       NaN       NaN       NaN       NaN       NaN    0.5356    0.1595    3.4553
    0.4203    0.7284    0.2718    0.7249       NaN       NaN       NaN       NaN       NaN    0.5145    0.1522    3.5298
    0.4224    0.7301    0.2704    0.7267       NaN       NaN       NaN       NaN       NaN    0.4954    0.1456    3.6161
    0.4246    0.7311    0.2689    0.7277       NaN       NaN       NaN       NaN       NaN    0.4780    0.1396    3.7212
    0.4267    0.7313    0.2675    0.7280       NaN       NaN       NaN       NaN       NaN    0.4620    0.1342    3.8582
    0.4288    0.7308    0.2661    0.7276       NaN       NaN       NaN       NaN       NaN    0.4474    0.1293    4.0531
    0.4309    0.7296    0.2647    0.7265       NaN       NaN       NaN       NaN       NaN    0.4341    0.1248    4.3751
    0.4331    0.7278    0.2633    0.7247       NaN       NaN       NaN       NaN       NaN    0.4217    0.1207    5.1858
    0.4352    0.7254    0.2620    0.7223       NaN       NaN       NaN       NaN       NaN    0.4104    0.1170    4.8875
    0.4373    0.7269    0.2606    0.7194    0.2869       NaN       NaN       NaN       NaN    0.3999    0.1135    3.2573
    0.4394    0.7276    0.2593    0.7159    0.3382       NaN       NaN       NaN       NaN    0.3902    0.1103    2.6032
    0.4415    0.7275    0.2579    0.7120    0.3700       NaN       NaN       NaN       NaN    0.3813    0.1074    2.2084
    0.4437    0.7266    0.2566    0.7077    0.3920       NaN       NaN       NaN       NaN    0.3730    0.1046    1.9360
    0.4458    0.7249    0.2553    0.7030    0.4079       NaN       NaN       NaN       NaN    0.3653    0.1021    1.7336
    0.4479    0.7226    0.2541    0.6980    0.4197       NaN       NaN       NaN       NaN    0.3582    0.0998    1.5760
    0.4500    0.7197    0.2528    0.6928    0.4287       NaN       NaN       NaN       NaN    0.3516    0.0976    1.4492
    0.4522    0.7163    0.2515    0.6874    0.4354       NaN       NaN       NaN       NaN    0.3455    0.0956    1.3444
    0.4543    0.7125    0.2503    0.6818    0.4406       NaN       NaN       NaN       NaN    0.3398    0.0938    1.2507
    0.4564    0.7083    0.2491    0.6761    0.4443       NaN       NaN       NaN       NaN    0.3346    0.0920    1.1761
    0.4585    0.7039    0.2478    0.6702    0.4471       NaN       NaN       NaN       NaN    0.3297    0.0904    1.1115
    0.4607    0.6993    0.2466    0.6643    0.4490       NaN       NaN       NaN       NaN    0.3252    0.0889    1.0550
    0.4628    0.6945    0.2454    0.6584    0.4502       NaN       NaN       NaN       NaN    0.3211    0.0876    1.0052
    0.4649    0.6895    0.2443    0.6524    0.4508       NaN       NaN       NaN       NaN    0.3173    0.0863    0.9608
    0.4670    0.6845    0.2431    0.6465    0.4510       NaN       NaN       NaN       NaN    0.3138    0.0851    0.9212
    0.4691    0.6794    0.2419    0.6405    0.4508       NaN       NaN       NaN       NaN    0.3106    0.0840    0.8855
    0.4713    0.6742    0.2408    0.6346    0.4503       NaN       NaN       NaN       NaN    0.3077    0.0830    0.8532
    0.4734    0.6691    0.2396    0.6287    0.4496       NaN       NaN       NaN       NaN    0.3051    0.0821    0.8240
    0.4755    0.6639    0.2385    0.6229    0.4486       NaN       NaN       NaN       NaN    0.3027    0.0812    0.7973
    0.4776    0.6587    0.2374    0.6171    0.4475       NaN       NaN       NaN       NaN    0.3006    0.0804    0.7729
    0.4798    0.6536    0.2363    0.6114    0.4463       NaN       NaN       NaN       NaN    0.2987    0.0797    0.7506
    0.4819    0.6485    0.2352    0.6057    0.4450       NaN       NaN       NaN       NaN    0.2971    0.0791    0.7301
    0.4840    0.6435    0.2341    0.6002    0.4436       NaN       NaN       NaN       NaN    0.2957    0.0785    0.7112
    0.4861    0.6385    0.2330    0.5947    0.4421       NaN       NaN       NaN       NaN    0.2945    0.0780    0.6938
    0.4883    0.6336    0.2320    0.5893    0.4406       NaN       NaN       NaN       NaN    0.2935    0.0775    0.6777
    0.4904    0.6288    0.2309    0.5840    0.4391       NaN       NaN       NaN       NaN    0.2928    0.0771    0.6628
    0.4925    0.6241    0.2299    0.5787    0.4376       NaN       NaN       NaN       NaN    0.2922    0.0768    0.6491
    0.4946    0.6194    0.2288    0.5736    0.4361       NaN       NaN       NaN       NaN    0.2919    0.0765    0.6363
    0.4967    0.6148    0.2278    0.5685    0.4346       NaN       NaN       NaN       NaN    0.2918    0.0763    0.6245
    0.4989    0.6103    0.2268    0.5636    0.4332       NaN       NaN       NaN       NaN    0.2919    0.0761    0.6135
    0.5010    0.6059    0.2258    0.5587    0.4318       NaN       NaN       NaN       NaN    0.2922    0.0760    0.6033
    0.5031    0.6016    0.2248    0.5539    0.4304       NaN       NaN       NaN       NaN    0.2927    0.0759    0.5939
    0.5052    0.5974    0.2238    0.5492    0.4290       NaN       NaN       NaN       NaN    0.2935    0.0759    0.5852
    0.5074    0.5932    0.2228    0.5445    0.4277       NaN       NaN       NaN       NaN    0.2944    0.0760    0.5771
    0.5095    0.5892    0.2219    0.5400    0.4264       NaN       NaN       NaN       NaN    0.2956    0.0760    0.5696
    0.5116    0.5853    0.2209    0.5355    0.4252       NaN       NaN       NaN       NaN    0.2970    0.0762    0.5627
    0.5137    0.5814    0.2199    0.5311    0.4241       NaN       NaN       NaN       NaN    0.2987    0.0764    0.5563
    0.5158    0.5777    0.2190    0.5268    0.4230       NaN       NaN       NaN       NaN    0.3005    0.0766    0.5505
    0.5180    0.5740    0.2181    0.5226    0.4220       NaN       NaN       NaN       NaN    0.3027    0.0770    0.5451
    0.5201    0.5704    0.2171    0.5184    0.4210       NaN       NaN       NaN       NaN    0.3050    0.0773    0.5402
    0.5222    0.5669    0.2162    0.5144    0.4201       NaN       NaN       NaN       NaN    0.3077    0.0777    0.5358
    0.5243    0.5636    0.2153    0.5103    0.4192       NaN       NaN       NaN       NaN    0.3106    0.0782    0.5318
    0.5265    0.5603    0.2144    0.5064    0.4185       NaN       NaN       NaN       NaN    0.3138    0.0788    0.5282
    0.5286    0.5571    0.2135    0.5025    0.4177       NaN       NaN       NaN       NaN    0.3173    0.0794    0.5250
    0.5307    0.5540    0.2126    0.4987    0.4171       NaN       NaN       NaN       NaN    0.3211    0.0801    0.5222
    0.5328    0.5509    0.2117    0.4950    0.4165       NaN       NaN       NaN       NaN    0.3252    0.0808    0.5198
    0.5350    0.5480    0.2108    0.4913    0.4160       NaN       NaN       NaN       NaN    0.3297    0.0816    0.5178
    0.5371    0.5452    0.2100    0.4877    0.4156       NaN       NaN       NaN       NaN    0.3346    0.0825    0.5161
    0.5392    0.5424    0.2091    0.4842    0.4152       NaN       NaN       NaN       NaN    0.3398    0.0835    0.5148
    0.5413    0.5398    0.2082    0.4807    0.4149       NaN       NaN       NaN       NaN    0.3455    0.0845    0.5139
    0.5434    0.5372    0.2074    0.4772    0.4147       NaN       NaN       NaN       NaN    0.3516    0.0857    0.5133
    0.5456    0.5347    0.2066    0.4739    0.4145       NaN       NaN       NaN       NaN    0.3583    0.0869    0.5131
    0.5477    0.5323    0.2057    0.4706    0.4144       NaN       NaN       NaN       NaN    0.3654    0.0883    0.5132
    0.5498    0.5300    0.2049    0.4673    0.4144       NaN       NaN       NaN       NaN    0.3732    0.0897    0.5137
    0.5519    0.5278    0.2041    0.4641    0.4144       NaN       NaN       NaN       NaN    0.3815    0.0913    0.5145
    0.5541    0.5257    0.2033    0.4609    0.4146       NaN       NaN       NaN       NaN    0.3905    0.0929    0.5158
    0.5562    0.5237    0.2024    0.4578    0.4148       NaN       NaN       NaN       NaN    0.4003    0.0947    0.5174
    0.5583    0.5217    0.2016    0.4548    0.4151       NaN       NaN       NaN       NaN    0.4109    0.0967    0.5194
    0.5604    0.5199    0.2009    0.4517    0.4154       NaN       NaN       NaN       NaN    0.4223    0.0988    0.5217
    0.5626    0.5181    0.2001    0.4488    0.4158       NaN       NaN       NaN       NaN    0.4348    0.1010    0.5245
    0.5647    0.5165    0.1993    0.4459    0.4163       NaN       NaN       NaN       NaN    0.4484    0.1035    0.5277
    0.5668    0.5149    0.1985    0.4430    0.4169       NaN       NaN       NaN       NaN    0.4632    0.1061    0.5308
    0.5689    0.5134    0.1977    0.4402    0.4176       NaN       NaN       NaN       NaN    0.4793    0.1090    0.5349
    0.5710    0.5120    0.1970    0.4374    0.4183       NaN       NaN       NaN       NaN    0.4971    0.1121    0.5395
    0.5732    0.5107    0.1962    0.4346    0.4191       NaN       NaN       NaN       NaN    0.5165    0.1155    0.5446
    0.5753    0.5095    0.1955    0.4319    0.4200       NaN       NaN       NaN       NaN    0.5381    0.1191    0.5502
    0.5774    0.5084    0.1947    0.4293    0.4209       NaN       NaN       NaN       NaN    0.5619    0.1231    0.5564
    0.5795    0.5074    0.1940    0.4266    0.4220       NaN       NaN       NaN       NaN    0.5884    0.1275    0.5632
    0.5817    0.5064    0.1932    0.4241    0.4231       NaN       NaN       NaN       NaN    0.6180    0.1323    0.5706
    0.5838    0.5056    0.1925    0.4215    0.4243       NaN       NaN       NaN       NaN    0.6513    0.1376    0.5788
    0.5859    0.5049    0.1918    0.4190    0.4256       NaN       NaN       NaN       NaN    0.6889    0.1434    0.5877
    0.5880    0.5043    0.1911    0.4165    0.4270       NaN       NaN       NaN       NaN    0.7318    0.1499    0.5974
    0.5901    0.5038    0.1903    0.4141    0.4285       NaN       NaN       NaN       NaN    0.7809    0.1571    0.6081
    0.5923    0.5034    0.1896    0.4117    0.4300       NaN       NaN       NaN       NaN    0.8378    0.1651    0.6198
    0.5944    0.5031    0.1889    0.4093    0.4317       NaN       NaN       NaN       NaN    0.9043    0.1742    0.6326
    0.5965    0.5029    0.1882    0.4069    0.4334       NaN       NaN       NaN       NaN    0.9829    0.1845    0.6467
    0.5986    0.5028    0.1875    0.4046    0.4352       NaN       NaN       NaN       NaN    1.0768    0.1962    0.6623
    0.6008    0.5029    0.1869    0.4023    0.4371       NaN       NaN       NaN       NaN    1.1907    0.2097    0.6796
    0.6029    0.5030    0.1862    0.4001    0.4391       NaN       NaN       NaN       NaN    1.3302    0.2254    0.6988
    0.6050    0.5033    0.1855    0.3979    0.4412       NaN       NaN       NaN       NaN    1.5025    0.2439    0.7205
    0.6071    0.5036    0.1848    0.3957    0.4434       NaN       NaN       NaN       NaN    1.7141    0.2658    0.7450
    0.6093    0.5042    0.1842    0.3935    0.4457       NaN       NaN       NaN       NaN    1.9643    0.2921    0.7731
    0.6114    0.5048    0.1835    0.3914    0.4481       NaN       NaN       NaN       NaN    2.2324    0.3237    0.8055
    0.6135    0.5055    0.1828    0.3893    0.4506       NaN       NaN       NaN       NaN    2.4693    0.3615    0.8431
    0.6156    0.5064    0.1822    0.3872    0.4532       NaN       NaN       NaN       NaN    2.6256    0.4054    0.8864
    0.6177    0.5074    0.1815    0.3851    0.4559       NaN       NaN       NaN       NaN    2.6881    0.4543    0.9357
    0.6199    0.5086    0.1809    0.3831    0.4587       NaN       NaN       NaN       NaN    2.6745    0.5066    0.9906
    0.6220    0.5098    0.1802    0.3811    0.4616       NaN       NaN       NaN       NaN    2.6083    0.5599    1.0503
    0.6241    0.5112    0.1796    0.3791    0.4647       NaN       NaN       NaN       NaN    2.5077    0.6120    1.1142
    0.6262    0.5128    0.1790    0.3771    0.4678       NaN       NaN       NaN       NaN    2.3848    0.6604    1.1812
    0.6284    0.5145    0.1784    0.3752    0.4710       NaN       NaN       NaN       NaN    2.2477    0.7025    1.2508
    0.6305    0.5163    0.1777    0.3733    0.4744       NaN       NaN       NaN       NaN    2.1012    0.7351    1.3223
    0.6326    0.5182    0.1771    0.3714    0.4778       NaN       NaN       NaN       NaN    1.9484    0.7544    1.3953
    0.6347    0.5203    0.1765    0.3695    0.4814       NaN       NaN       NaN       NaN    1.7913    0.7554    1.4699
    0.6368    0.5226    0.1759    0.3677    0.4850       NaN       NaN       NaN       NaN    1.6307    0.7316    1.5467
    0.6390    0.5250    0.1753    0.3659    0.4888       NaN       NaN       NaN       NaN    1.4662    0.6733    1.6286
    0.6411    0.5275    0.1747    0.3641    0.4927       NaN       NaN       NaN       NaN    1.3043    0.5847    1.7373
    0.6432    0.5302    0.1741    0.3623    0.4967       NaN       NaN       NaN       NaN    1.1641    0.5008    1.8264
    0.6453    0.5330    0.1735    0.3605    0.5008       NaN       NaN       NaN       NaN    1.0489    0.4327    1.9015
    0.6475    0.5360    0.1729    0.3588    0.5050       NaN       NaN       NaN       NaN    0.9549    0.3790    1.9679
    0.6496    0.5391    0.1723    0.3570    0.5093       NaN       NaN       NaN       NaN    0.8772    0.3366    2.0284
    0.6517    0.5423    0.1717    0.3553    0.5138       NaN       NaN       NaN       NaN    0.8122    0.3025    2.0850
    0.6538    0.5456    0.1712    0.3537    0.5183       NaN       NaN       NaN       NaN    0.7571    0.2748    2.1388
    0.6560    0.5491    0.1706    0.3520    0.5229       NaN       NaN       NaN       NaN    0.7097    0.2517    2.1904
    0.6581    0.5527    0.1700    0.3503    0.5275       NaN       NaN       NaN       NaN    0.6686    0.2324    2.2403
    0.6602    0.5565    0.1695    0.3487    0.5323       NaN       NaN       NaN       NaN    0.6326    0.2159    2.2888
    0.6623    0.5603    0.1689    0.3471    0.5371       NaN       NaN       NaN       NaN    0.6009    0.2018    2.3362
    0.6644    0.5642    0.1683    0.3455    0.5419       NaN       NaN       NaN       NaN    0.5727    0.1895    2.3825
    0.6666    0.5682    0.1678    0.3439    0.5469       NaN       NaN       NaN       NaN    0.5475    0.1787    2.4281
    0.6687    0.5722    0.1672    0.3424    0.5518       NaN       NaN       NaN       NaN    0.5248    0.1692    2.4730
    0.6708    0.5763    0.1667    0.3408    0.5567       NaN       NaN       NaN       NaN    0.5043    0.1608    2.5173
    0.6729    0.5805    0.1661    0.3393    0.5617       NaN       NaN       NaN       NaN    0.4858    0.1532    2.5611
    0.6751    0.5847    0.1656    0.3378    0.5666       NaN       NaN       NaN       NaN    0.4689    0.1465    2.6047
    0.6772    0.5888    0.1651    0.3363    0.5715       NaN       NaN       NaN       NaN    0.4536    0.1404    2.6481
    0.6793    0.5930    0.1645    0.3348    0.5763       NaN       NaN       NaN       NaN    0.4395    0.1349    2.6909
    0.6814    0.5971    0.1640    0.3333    0.5811       NaN       NaN       NaN       NaN    0.4266    0.1299    2.7345
    0.6836    0.6011    0.1635    0.3319    0.5857       NaN       NaN       NaN       NaN    0.4148    0.1253    2.7785
    0.6857    0.6050    0.1629    0.3304    0.5902       NaN       NaN       NaN       NaN    0.4038    0.1211    2.8233
    0.6878    0.6088    0.1624    0.3290    0.5946       NaN       NaN       NaN       NaN    0.3938    0.1173    2.8693
    0.6899    0.6125    0.1619    0.3276    0.5987       NaN       NaN       NaN       NaN    0.3845    0.1138    2.9169
    0.6920    0.6159    0.1614    0.3262    0.6027       NaN       NaN       NaN       NaN    0.3759    0.1106    2.9670
    0.6942    0.6192    0.1609    0.3248    0.6064       NaN       NaN       NaN       NaN    0.3680    0.1076    3.0202
    0.6963    0.6222    0.1604    0.3234    0.6098       NaN       NaN       NaN       NaN    0.3606    0.1049    3.0779
    0.6984    0.6250    0.1599    0.3221    0.6130       NaN       NaN       NaN       NaN    0.3538    0.1023    3.1416
    0.7005    0.6275    0.1594    0.3207    0.6159       NaN       NaN       NaN       NaN    0.3475    0.1000    3.2139
    0.7027    0.6297    0.1589    0.3194    0.6184       NaN       NaN       NaN       NaN    0.3417    0.0978    3.2981
    0.7048    0.6316    0.1584    0.3181    0.6206       NaN       NaN       NaN       NaN    0.3363    0.0958    3.3999
    0.7069    0.6332    0.1579    0.3168    0.6224       NaN       NaN       NaN       NaN    0.3314    0.0939    3.5287
    0.7090    0.6344    0.1574    0.3155    0.6239       NaN       NaN       NaN       NaN    0.3268    0.0922    3.7026
    0.7111    0.6353    0.1569    0.3142    0.6250       NaN       NaN       NaN       NaN    0.3225    0.0905    3.9623
    0.7133    0.6358    0.1564    0.3129    0.6257       NaN       NaN       NaN       NaN    0.3186    0.0890    4.4365
    0.7154    0.6360    0.1560    0.3117    0.6260       NaN       NaN       NaN       NaN    0.3151    0.0877    5.2780
    0.7175    0.6392    0.1555    0.3104    0.6260    0.2427       NaN       NaN       NaN    0.3118    0.0864    3.3047
    0.7196    0.6432    0.1550    0.3092    0.6256    0.3012       NaN       NaN       NaN    0.3088    0.0852    2.5309
    0.7218    0.6464    0.1545    0.3080    0.6249    0.3346       NaN       NaN       NaN    0.3061    0.0841    2.1056
    0.7239    0.6487    0.1541    0.3067    0.6238    0.3572       NaN       NaN       NaN    0.3037    0.0831    1.8268
    0.7260    0.6504    0.1536    0.3055    0.6224    0.3739       NaN       NaN       NaN    0.3015    0.0821    1.6265
    0.7281    0.6513    0.1531    0.3043    0.6207    0.3868       NaN       NaN       NaN    0.2996    0.0813    1.4739
    0.7303    0.6517    0.1527    0.3032    0.6188    0.3969       NaN       NaN       NaN    0.2979    0.0805    1.3531
    0.7324    0.6516    0.1522    0.3020    0.6165    0.4050       NaN       NaN       NaN    0.2965    0.0798    1.2546
    0.7345    0.6510    0.1518    0.3008    0.6141    0.4116       NaN       NaN       NaN    0.2953    0.0792    1.1725
    0.7366    0.6500    0.1513    0.2997    0.6114    0.4169       NaN       NaN       NaN    0.2943    0.0786    1.1030
    0.7387    0.6487    0.1509    0.2985    0.6085    0.4212       NaN       NaN       NaN    0.2935    0.0781    1.0431
    0.7409    0.6470    0.1504    0.2974    0.6055    0.4246       NaN       NaN       NaN    0.2930    0.0776    0.9911
    0.7430    0.6450    0.1500    0.2963    0.6023    0.4273       NaN       NaN       NaN    0.2927    0.0772    0.9454
    0.7451    0.6428    0.1495    0.2952    0.5989    0.4294       NaN       NaN       NaN    0.2925    0.0769    0.9050
    0.7472    0.6404    0.1491    0.2940    0.5955    0.4309       NaN       NaN       NaN    0.2926    0.0766    0.8690
    0.7494    0.6378    0.1486    0.2930    0.5919    0.4321       NaN       NaN       NaN    0.2930    0.0764    0.8367
    0.7515    0.6350    0.1482    0.2919    0.5882    0.4328       NaN       NaN       NaN    0.2935    0.0762    0.8077
    0.7536    0.6321    0.1478    0.2908    0.5845    0.4333       NaN       NaN       NaN    0.2942    0.0761    0.7814
    0.7557    0.6290    0.1474    0.2897    0.5808    0.4334       NaN       NaN       NaN    0.2952    0.0760    0.7575
    0.7579    0.6259    0.1469    0.2887    0.5769    0.4333       NaN       NaN       NaN    0.2964    0.0760    0.7358
    0.7600    0.6227    0.1465    0.2876    0.5731    0.4331       NaN       NaN       NaN    0.2978    0.0761    0.7159
    0.7621    0.6194    0.1461    0.2866    0.5692    0.4326       NaN       NaN       NaN    0.2995    0.0762    0.6978
    0.7642    0.6161    0.1457    0.2855    0.5653    0.4320       NaN       NaN       NaN    0.3014    0.0763    0.6811
    0.7663    0.6127    0.1452    0.2845    0.5614    0.4313       NaN       NaN       NaN    0.3035    0.0765    0.6658
    0.7685    0.6093    0.1448    0.2835    0.5575    0.4305       NaN       NaN       NaN    0.3060    0.0768    0.6517
    0.7706    0.6059    0.1444    0.2825    0.5536    0.4296       NaN       NaN       NaN    0.3086    0.0771    0.6388
    0.7727    0.6026    0.1440    0.2815    0.5497    0.4287       NaN       NaN       NaN    0.3116    0.0775    0.6268
    0.7748    0.5992    0.1436    0.2805    0.5459    0.4277       NaN       NaN       NaN    0.3149    0.0780    0.6159
    0.7770    0.5958    0.1432    0.2795    0.5421    0.4267       NaN       NaN       NaN    0.3184    0.0785    0.6057
    0.7791    0.5924    0.1428    0.2785    0.5383    0.4256       NaN       NaN       NaN    0.3223    0.0790    0.5964
    0.7812    0.5891    0.1424    0.2775    0.5345    0.4245       NaN       NaN       NaN    0.3265    0.0796    0.5879
    0.7833    0.5858    0.1420    0.2766    0.5307    0.4235       NaN       NaN       NaN    0.3311    0.0803    0.5800
    0.7854    0.5825    0.1416    0.2756    0.5270    0.4224       NaN       NaN       NaN    0.3361    0.0811    0.5728
    0.7876    0.5793    0.1412    0.2747    0.5234    0.4213       NaN       NaN       NaN    0.3415    0.0819    0.5663
    0.7897    0.5761    0.1408    0.2737    0.5197    0.4203       NaN       NaN       NaN    0.3473    0.0829    0.5603
    0.7918    0.5730    0.1404    0.2728    0.5161    0.4192       NaN       NaN       NaN    0.3536    0.0839    0.5549
    0.7939    0.5699    0.1400    0.2719    0.5126    0.4182       NaN       NaN       NaN    0.3604    0.0849    0.5494
    0.7961    0.5669    0.1396    0.2709    0.5091    0.4172       NaN       NaN       NaN    0.3678    0.0861    0.5450
    0.7982    0.5639    0.1393    0.2700    0.5056    0.4163       NaN       NaN       NaN    0.3757    0.0874    0.5412
    0.8003    0.5609    0.1389    0.2691    0.5022    0.4154       NaN       NaN       NaN    0.3843    0.0888    0.5378
    0.8024    0.5581    0.1385    0.2682    0.4989    0.4145       NaN       NaN       NaN    0.3936    0.0902    0.5349
    0.8046    0.5553    0.1381    0.2673    0.4955    0.4137       NaN       NaN       NaN    0.4037    0.0918    0.5324
    0.8067    0.5525    0.1377    0.2664    0.4923    0.4129       NaN       NaN       NaN    0.4147    0.0936    0.5304
    0.8088    0.5498    0.1374    0.2656    0.4890    0.4121       NaN       NaN       NaN    0.4266    0.0954    0.5288
    0.8109    0.5472    0.1370    0.2647    0.4858    0.4115       NaN       NaN       NaN    0.4396    0.0975    0.5277
    0.8130    0.5446    0.1366    0.2638    0.4827    0.4108       NaN       NaN       NaN    0.4537    0.0996    0.5270
    0.8152    0.5421    0.1363    0.2630    0.4796    0.4102       NaN       NaN       NaN    0.4692    0.1020    0.5267
    0.8173    0.5396    0.1359    0.2621    0.4765    0.4097       NaN       NaN       NaN    0.4862    0.1045    0.5269
    0.8194    0.5372    0.1355    0.2612    0.4735    0.4092       NaN       NaN       NaN    0.5049    0.1073    0.5275
    0.8215    0.5349    0.1352    0.2604    0.4705    0.4088       NaN       NaN       NaN    0.5255    0.1103    0.5286
    0.8237    0.5326    0.1348    0.2596    0.4676    0.4084       NaN       NaN       NaN    0.5484    0.1135    0.5302
    0.8258    0.5304    0.1344    0.2587    0.4647    0.4081       NaN       NaN       NaN    0.5738    0.1171    0.5322
    0.8279    0.5283    0.1341    0.2579    0.4618    0.4078       NaN       NaN       NaN    0.6023    0.1210    0.5347
    0.8300    0.5262    0.1337    0.2571    0.4590    0.4076       NaN       NaN       NaN    0.6343    0.1252    0.5378
    0.8321    0.5242    0.1334    0.2563    0.4563    0.4074       NaN       NaN       NaN    0.6706    0.1299    0.5414
    0.8343    0.5223    0.1330    0.2555    0.4535    0.4073       NaN       NaN       NaN    0.7120    0.1350    0.5456
    0.8364    0.5204    0.1327    0.2547    0.4508    0.4073       NaN       NaN       NaN    0.7596    0.1407    0.5505
    0.8385    0.5186    0.1323    0.2539    0.4482    0.4073       NaN       NaN       NaN    0.8148    0.1470    0.5561
    0.8406    0.5169    0.1320    0.2531    0.4456    0.4074       NaN       NaN       NaN    0.8796    0.1541    0.5624
    0.8428    0.5152    0.1316    0.2523    0.4430    0.4076       NaN       NaN       NaN    0.9564    0.1620    0.5695
    0.8449    0.5136    0.1313    0.2515    0.4404    0.4078       NaN       NaN       NaN    1.0487    0.1711    0.5776
    0.8470    0.5121    0.1310    0.2507    0.4379    0.4081       NaN       NaN       NaN    1.1611    0.1814    0.5868
    0.8491    0.5106    0.1306    0.2499    0.4354    0.4084       NaN       NaN       NaN    1.2993    0.1932    0.5973
    0.8513    0.5092    0.1303    0.2492    0.4330    0.4088       NaN       NaN       NaN    1.4698    0.2071    0.6092
    0.8534    0.5079    0.1299    0.2484    0.4306    0.4093       NaN       NaN       NaN    1.6762    0.2234    0.6229
    0.8555    0.5066    0.1296    0.2477    0.4282    0.4098       NaN       NaN       NaN    1.9109    0.2427    0.6388
    0.8576    0.5054    0.1293    0.2469    0.4259    0.4104       NaN       NaN       NaN    2.1434    0.2658    0.6574
    0.8597    0.5043    0.1289    0.2462    0.4236    0.4111       NaN       NaN       NaN    2.3284    0.2929    0.6792
    0.8619    0.5033    0.1286    0.2454    0.4213    0.4118       NaN       NaN       NaN    2.4403    0.3241    0.7044
    0.8640    0.5023    0.1283    0.2447    0.4190    0.4126       NaN       NaN       NaN    2.4841    0.3588    0.7331
    0.8661    0.5014    0.1280    0.2439    0.4168    0.4134       NaN       NaN       NaN    2.4778    0.3962    0.7653
    0.8682    0.5006    0.1276    0.2432    0.4146    0.4144       NaN       NaN       NaN    2.4378    0.4352    0.8005
    0.8704    0.4998    0.1273    0.2425    0.4124    0.4154       NaN       NaN       NaN    2.3754    0.4750    0.8383
    0.8725    0.4992    0.1270    0.2418    0.4103    0.4164       NaN       NaN       NaN    2.2980    0.5146    0.8785
    0.8746    0.4986    0.1267    0.2411    0.4082    0.4176       NaN       NaN       NaN    2.2105    0.5533    0.9206
    0.8767    0.4981    0.1264    0.2403    0.4061    0.4188       NaN       NaN       NaN    2.1159    0.5901    0.9642
    0.8789    0.4977    0.1260    0.2396    0.4041    0.4200       NaN       NaN       NaN    2.0164    0.6242    1.0091
    0.8810    0.4973    0.1257    0.2389    0.4020    0.4214       NaN       NaN       NaN    1.9134    0.6544    1.0549
    0.8831    0.4970    0.1254    0.2382    0.4000    0.4228       NaN       NaN       NaN    1.8076    0.6793    1.1016
    0.8852    0.4969    0.1251    0.2376    0.3980    0.4243       NaN       NaN       NaN    1.6995    0.6972    1.1489
    0.8873    0.4968    0.1248    0.2369    0.3961    0.4259       NaN       NaN       NaN    1.5893    0.7057    1.1970
    0.8895    0.4967    0.1245    0.2362    0.3942    0.4275       NaN       NaN       NaN    1.4765    0.7013    1.2460
    0.8916    0.4968    0.1242    0.2355    0.3922    0.4292       NaN       NaN       NaN    1.3603    0.6795    1.2968
    0.8937    0.4970    0.1239    0.2348    0.3904    0.4310       NaN       NaN       NaN    1.2391    0.6332    1.3513
    0.8958    0.4972    0.1236    0.2342    0.3885    0.4329       NaN       NaN       NaN    1.1118    0.5542    1.4318
    0.8980    0.4976    0.1233    0.2335    0.3867    0.4348       NaN       NaN       NaN    0.9977    0.4743    1.5091
    0.9001    0.4980    0.1230    0.2328    0.3848    0.4368       NaN       NaN       NaN    0.9040    0.4094    1.5682
    0.9022    0.4985    0.1227    0.2322    0.3830    0.4389       NaN       NaN       NaN    0.8280    0.3587    1.6183
    0.9043    0.4992    0.1224    0.2315    0.3813    0.4411       NaN       NaN       NaN    0.7654    0.3188    1.6628
    0.9064    0.4999    0.1221    0.2309    0.3795    0.4433       NaN       NaN       NaN    0.7132    0.2870    1.7036
    0.9086    0.5007    0.1218    0.2302    0.3778    0.4457       NaN       NaN       NaN    0.6688    0.2612    1.7422
    0.9107    0.5016    0.1215    0.2296    0.3761    0.4481       NaN       NaN       NaN    0.6307    0.2398    1.7791
    0.9128    0.5026    0.1212    0.2289    0.3744    0.4505       NaN       NaN       NaN    0.5975    0.2218    1.8147
    0.9149    0.5037    0.1209    0.2283    0.3727    0.4531       NaN       NaN       NaN    0.5684    0.2065    1.8494
    0.9171    0.5049    0.1206    0.2277    0.3710    0.4557       NaN       NaN       NaN    0.5426    0.1933    1.8833
    0.9192    0.5062    0.1203    0.2270    0.3694    0.4584       NaN       NaN       NaN    0.5196    0.1818    1.9166
    0.9213    0.5076    0.1200    0.2264    0.3678    0.4612       NaN       NaN       NaN    0.4990    0.1718    1.9494
    0.9234    0.5091    0.1197    0.2258    0.3662    0.4641       NaN       NaN       NaN    0.4805    0.1630    1.9818
    0.9256    0.5107    0.1194    0.2252    0.3646    0.4670       NaN       NaN       NaN    0.4637    0.1551    2.0139
    0.9277    0.5123    0.1192    0.2246    0.3630    0.4700       NaN       NaN       NaN    0.4484    0.1481    2.0457
    0.9298    0.5141    0.1189    0.2239    0.3614    0.4731       NaN       NaN       NaN    0.4345    0.1417    2.0773
    0.9319    0.5160    0.1186    0.2233    0.3599    0.4762       NaN       NaN       NaN    0.4218    0.1361    2.1087
    0.9340    0.5179    0.1183    0.2227    0.3584    0.4794       NaN       NaN       NaN    0.4101    0.1309    2.1400
    0.9362    0.5200    0.1180    0.2221    0.3569    0.4826       NaN       NaN       NaN    0.3994    0.1262    2.1712
    0.9383    0.5221    0.1178    0.2215    0.3554    0.4859       NaN       NaN       NaN    0.3896    0.1219    2.2024
    0.9404    0.5243    0.1175    0.2210    0.3539    0.4893       NaN       NaN       NaN    0.3805    0.1180    2.2336
    0.9425    0.5266    0.1172    0.2204    0.3525    0.4927       NaN       NaN       NaN    0.3721    0.1144    2.2649
    0.9447    0.5290    0.1169    0.2198    0.3510    0.4961       NaN       NaN       NaN    0.3644    0.1111    2.2964
    0.9468    0.5314    0.1167    0.2192    0.3496    0.4996       NaN       NaN       NaN    0.3573    0.1081    2.3280
    0.9489    0.5339    0.1164    0.2186    0.3482    0.5031       NaN       NaN       NaN    0.3507    0.1053    2.3599
    0.9510    0.5364    0.1161    0.2180    0.3468    0.5066       NaN       NaN       NaN    0.3446    0.1027    2.3921
    0.9532    0.5390    0.1158    0.2175    0.3454    0.5102       NaN       NaN       NaN    0.3390    0.1003    2.4248
    0.9553    0.5416    0.1156    0.2169    0.3440    0.5137       NaN       NaN       NaN    0.3338    0.0981    2.4580
    0.9574    0.5442    0.1153    0.2163    0.3426    0.5172       NaN       NaN       NaN    0.3291    0.0960    2.4920
    0.9595    0.5469    0.1150    0.2158    0.3413    0.5207       NaN       NaN       NaN    0.3247    0.0941    2.5267
    0.9616    0.5495    0.1148    0.2152    0.3400    0.5242       NaN       NaN       NaN    0.3206    0.0924    2.5625
    0.9638    0.5521    0.1145    0.2146    0.3386    0.5276       NaN       NaN       NaN    0.3169    0.0908    2.5996
    0.9659    0.5548    0.1143    0.2141    0.3373    0.5310       NaN       NaN       NaN    0.3135    0.0892    2.6382
    0.9680    0.5573    0.1140    0.2135    0.3360    0.5343       NaN       NaN       NaN    0.3105    0.0878    2.6787
    0.9701    0.5599    0.1137    0.2130    0.3348    0.5375       NaN       NaN       NaN    0.3077    0.0866    2.7217
    0.9723    0.5624    0.1135    0.2124    0.3335    0.5406       NaN       NaN       NaN    0.3052    0.0854    2.7677
    0.9744    0.5647    0.1132    0.2119    0.3322    0.5437       NaN       NaN       NaN    0.3029    0.0842    2.8174
    0.9765    0.5671    0.1130    0.2114    0.3310    0.5466       NaN       NaN       NaN    0.3009    0.0832    2.8720
    0.9786    0.5692    0.1127    0.2108    0.3297    0.5493       NaN       NaN       NaN    0.2992    0.0823    2.9329
    0.9807    0.5713    0.1124    0.2103    0.3285    0.5519       NaN       NaN       NaN    0.2977    0.0814    3.0021
    0.9829    0.5733    0.1122    0.2098    0.3273    0.5543       NaN       NaN       NaN    0.2965    0.0807    3.0826
    0.9850    0.5751    0.1119    0.2092    0.3261    0.5566       NaN       NaN       NaN    0.2954    0.0799    3.1791
    0.9871    0.5767    0.1117    0.2087    0.3249    0.5587       NaN       NaN       NaN    0.2946    0.0793    3.2992
    0.9892    0.5781    0.1114    0.2082    0.3237    0.5605       NaN       NaN       NaN    0.2941    0.0787    3.4570
    0.9914    0.5794    0.1112    0.2077    0.3226    0.5622       NaN       NaN       NaN    0.2937    0.0782    3.6823
    0.9935    0.5805    0.1109    0.2071    0.3214    0.5636       NaN       NaN       NaN    0.2936    0.0778    4.0597
    0.9956    0.5814    0.1107    0.2066    0.3203    0.5648       NaN       NaN       NaN    0.2937    0.0774    5.1406
    0.9977    0.5821    0.1105    0.2061    0.3191    0.5658       NaN       NaN       NaN    0.2940    0.0770    3.7257
    0.9999    0.5890    0.1102    0.2056    0.3180    0.5665    0.2680       NaN       NaN    0.2945    0.0768    2.5935
    1.0020    0.5933    0.1100    0.2051    0.3169    0.5671    0.3041       NaN       NaN    0.2953    0.0765    2.0884
    1.0041    0.5969    0.1097    0.2046    0.3157    0.5673    0.3275       NaN       NaN    0.2963    0.0764    1.7831
    1.0062    0.5998    0.1095    0.2041    0.3146    0.5674    0.3446       NaN       NaN    0.2975    0.0763    1.5732
    1.0083    0.6021    0.1092    0.2036    0.3136    0.5672    0.3579       NaN       NaN    0.2990    0.0762    1.4180
    1.0105    0.6039    0.1090    0.2031    0.3125    0.5668    0.3686       NaN       NaN    0.3007    0.0762    1.2976
    1.0126    0.6053    0.1088    0.2026    0.3114    0.5661    0.3775       NaN       NaN    0.3026    0.0763    1.2009
    1.0147    0.6063    0.1085    0.2021    0.3103    0.5653    0.3848       NaN       NaN    0.3048    0.0764    1.1212
    1.0168    0.6069    0.1083    0.2016    0.3093    0.5643    0.3910       NaN       NaN    0.3073    0.0765    1.0544
    1.0190    0.6072    0.1080    0.2011    0.3082    0.5630    0.3963       NaN       NaN    0.3100    0.0768    0.9949
    1.0211    0.6072    0.1078    0.2006    0.3072    0.5616    0.4007       NaN       NaN    0.3130    0.0770    0.9459
    1.0232    0.6069    0.1076    0.2002    0.3062    0.5600    0.4044       NaN       NaN    0.3164    0.0774    0.9032
    1.0253    0.6063    0.1073    0.1997    0.3051    0.5583    0.4076       NaN       NaN    0.3200    0.0778    0.8655
    1.0275    0.6055    0.1071    0.1992    0.3041    0.5564    0.4102       NaN       NaN    0.3240    0.0782    0.8322
    1.0296    0.6044    0.1069    0.1987    0.3031    0.5543    0.4124       NaN       NaN    0.3284    0.0787    0.8024
    1.0317    0.6032    0.1067    0.1983    0.3021    0.5522    0.4142       NaN       NaN    0.3331    0.0793    0.7758
    1.0338    0.6018    0.1064    0.1978    0.3011    0.5499    0.4156       NaN       NaN    0.3382    0.0800    0.7518
    1.0359    0.6002    0.1062    0.1973    0.3001    0.5475    0.4167       NaN       NaN    0.3438    0.0807    0.7301
    1.0381    0.5985    0.1060    0.1969    0.2992    0.5450    0.4176       NaN       NaN    0.3498    0.0815    0.7105
    1.0402    0.5966    0.1057    0.1964    0.2982    0.5424    0.4182       NaN       NaN    0.3563    0.0824    0.6926
    1.0423    0.5946    0.1055    0.1959    0.2972    0.5398    0.4186       NaN       NaN    0.3634    0.0833    0.6764
    1.0444    0.5925    0.1053    0.1955    0.2963    0.5371    0.4188       NaN       NaN    0.3710    0.0843    0.6616
    1.0466    0.5903    0.1051    0.1950    0.2954    0.5343    0.4188       NaN       NaN    0.3793    0.0855    0.6482
    1.0487    0.5881    0.1048    0.1946    0.2944    0.5315    0.4187       NaN       NaN    0.3882    0.0867    0.6359
    1.0508    0.5858    0.1046    0.1941    0.2935    0.5287    0.4185       NaN       NaN    0.3979    0.0880    0.6247
    1.0529    0.5834    0.1044    0.1937    0.2926    0.5258    0.4182       NaN       NaN    0.4085    0.0894    0.6145
    1.0550    0.5810    0.1042    0.1932    0.2916    0.5229    0.4177       NaN       NaN    0.4200    0.0910    0.6053
    1.0572    0.5785    0.1040    0.1928    0.2907    0.5199    0.4172       NaN       NaN    0.4325    0.0926    0.5968
    1.0593    0.5760    0.1037    0.1923    0.2898    0.5170    0.4166       NaN       NaN    0.4462    0.0945    0.5893
    1.0614    0.5735    0.1035    0.1919    0.2889    0.5140    0.4160       NaN       NaN    0.4612    0.0964    0.5824
    1.0635    0.5710    0.1033    0.1914    0.2880    0.5111    0.4153       NaN       NaN    0.4776    0.0985    0.5763
    1.0657    0.5685    0.1031    0.1910    0.2872    0.5081    0.4146       NaN       NaN    0.4957    0.1008    0.5709
    1.0678    0.5660    0.1029    0.1906    0.2863    0.5051    0.4139       NaN       NaN    0.5158    0.1033    0.5661
    1.0699    0.5634    0.1027    0.1901    0.2854    0.5022    0.4131       NaN       NaN    0.5380    0.1060    0.5620
    1.0720    0.5609    0.1025    0.1897    0.2846    0.4992    0.4124       NaN       NaN    0.5628    0.1089    0.5585
    1.0742    0.5584    0.1022    0.1893    0.2837    0.4963    0.4116       NaN       NaN    0.5906    0.1121    0.5556
    1.0763    0.5559    0.1020    0.1889    0.2828    0.4934    0.4108       NaN       NaN    0.6219    0.1156    0.5534
    1.0784    0.5535    0.1018    0.1884    0.2820    0.4905    0.4101       NaN       NaN    0.6575    0.1194    0.5517
    1.0805    0.5510    0.1016    0.1880    0.2812    0.4876    0.4093       NaN       NaN    0.6982    0.1236    0.5506
    1.0826    0.5486    0.1014    0.1876    0.2803    0.4847    0.4085       NaN       NaN    0.7452    0.1282    0.5502
    1.0848    0.5462    0.1012    0.1872    0.2795    0.4819    0.4078       NaN       NaN    0.7999    0.1334    0.5504
    1.0869    0.5439    0.1010    0.1867    0.2787    0.4791    0.4071       NaN       NaN    0.8644    0.1391    0.5513
    1.0890    0.5416    0.1008    0.1863    0.2779    0.4763    0.4064       NaN       NaN    0.9414    0.1455    0.5530
    1.0911    0.5393    0.1006    0.1859    0.2770    0.4735    0.4058       NaN       NaN    1.0344    0.1527    0.5554
    1.0933    0.5370    0.1004    0.1855    0.2762    0.4708    0.4052       NaN       NaN    1.1483    0.1610    0.5586
    1.0954    0.5348    0.1002    0.1851    0.2754    0.4681    0.4046       NaN       NaN    1.2887    0.1704    0.5628
    1.0975    0.5327    0.1000    0.1847    0.2747    0.4654    0.4040       NaN       NaN    1.4608    0.1814    0.5682
    1.0996    0.5305    0.0998    0.1843    0.2739    0.4627    0.4035       NaN       NaN    1.6639    0.1942    0.5749
    1.1017    0.5285    0.0996    0.1839    0.2731    0.4601    0.4030       NaN       NaN    1.8808    0.2094    0.5832
    1.1039    0.5264    0.0994    0.1835    0.2723    0.4575    0.4026       NaN       NaN    2.0753    0.2273    0.5934
    1.1060    0.5244    0.0992    0.1831    0.2715    0.4549    0.4021       NaN       NaN    2.2143    0.2482    0.6058
    1.1081    0.5225    0.0990    0.1827    0.2708    0.4524    0.4018       NaN       NaN    2.2918    0.2718    0.6205
    1.1102    0.5206    0.0988    0.1823    0.2700    0.4499    0.4015       NaN       NaN    2.3203    0.2979    0.6376
    1.1124    0.5188    0.0986    0.1819    0.2692    0.4474    0.4012       NaN       NaN    2.3146    0.3260    0.6570
    1.1145    0.5170    0.0985    0.1815    0.2685    0.4450    0.4010       NaN       NaN    2.2862    0.3554    0.6785
    1.1166    0.5152    0.0983    0.1811    0.2677    0.4425    0.4008       NaN       NaN    2.2425    0.3857    0.7021
    1.1187    0.5135    0.0981    0.1807    0.2670    0.4401    0.4007       NaN       NaN    2.1884    0.4165    0.7273
    1.1209    0.5119    0.0978    0.1803    0.2663    0.4378    0.4006       NaN       NaN    2.1269    0.4474    0.7541
    1.1230    0.5102    0.0889    0.1799    0.2655    0.4354    0.4005       NaN       NaN    2.0602    0.4780    0.7822
    1.1251    0.5086    0.0802    0.1795    0.2648    0.4331    0.4006       NaN       NaN    1.9897    0.5079    0.8115
    1.1272    0.5071    0.0737    0.1792    0.2641    0.4308    0.4006       NaN       NaN    1.9164    0.5368    0.8417
    1.1293    0.5057    0.0685    0.1788    0.2634    0.4286    0.4007       NaN       NaN    1.8411    0.5642    0.8726
    1.1315    0.5043    0.0643    0.1784    0.2626    0.4264    0.4009       NaN       NaN    1.7641    0.5897    0.9042
    1.1336    0.5029    0.0608    0.1780    0.2619    0.4242    0.4011       NaN       NaN    1.6858    0.6128    0.9362
    1.1357    0.5017    0.0578    0.1776    0.2612    0.4220    0.4014       NaN       NaN    1.6063    0.6328    0.9688
    1.1378    0.5004    0.0552    0.1773    0.2605    0.4198    0.4017       NaN       NaN    1.5257    0.6486    1.0017
    1.1400    0.4993    0.0529    0.1769    0.2598    0.4177    0.4021       NaN       NaN    1.4439    0.6592    1.0350
    1.1421    0.4982    0.0509    0.1765    0.2591    0.4156    0.4026       NaN       NaN    1.3603    0.6627    1.0687
    1.1442    0.4971    0.0491    0.1762    0.2584    0.4135    0.4030       NaN       NaN    1.2745    0.6568    1.1031
    1.1463    0.4961    0.0475    0.1758    0.2578    0.4115    0.4036       NaN       NaN    1.1852    0.6379    1.1389
    1.1485    0.4952    0.0460    0.1754    0.2571    0.4095    0.4042       NaN       NaN    1.0905    0.6002    1.1774
    1.1506    0.4943    0.0447    0.1750    0.2564    0.4075    0.4048       NaN       NaN    0.9869    0.5335    1.2324
    1.1527    0.4935    0.0435    0.1747    0.2557    0.4055    0.4056       NaN       NaN    0.8877    0.4565    1.3081
    1.1548    0.4927    0.0423    0.1743    0.2551    0.4036    0.4063       NaN       NaN    0.8061    0.3931    1.3595
    1.1569    0.4920    0.0413    0.1740    0.2544    0.4016    0.4072       NaN       NaN    0.7405    0.3439    1.4008
    1.1591    0.4914    0.0403    0.1736    0.2537    0.3997    0.4080       NaN       NaN    0.6871    0.3056    1.4363
    1.1612    0.4908    0.0394    0.1732    0.2531    0.3979    0.4090       NaN       NaN    0.6427    0.2752    1.4682
    1.1633    0.4903    0.0386    0.1729    0.2524    0.3960    0.4100       NaN       NaN    0.6052    0.2506    1.4978
    1.1654    0.4899    0.0378    0.1725    0.2518    0.3942    0.4110       NaN       NaN    0.5731    0.2303    1.5258
    1.1676    0.4895    0.0370    0.1722    0.2511    0.3924    0.4122       NaN       NaN    0.5451    0.2133    1.5527
    1.1697    0.4892    0.0363    0.1718    0.2505    0.3906    0.4133       NaN       NaN    0.5206    0.1988    1.5787
    1.1718    0.4889    0.0356    0.1715    0.2499    0.3888    0.4146       NaN       NaN    0.4989    0.1864    1.6042
    1.1739    0.4888    0.0350    0.1711    0.2492    0.3870    0.4159       NaN       NaN    0.4796    0.1756    1.6291
    1.1760    0.4887    0.0344    0.1708    0.2486    0.3853    0.4172       NaN       NaN    0.4622    0.1661    1.6537
    1.1782    0.4886    0.0338    0.1704    0.2480    0.3836    0.4187       NaN       NaN    0.4466    0.1578    1.6779
    1.1803    0.4886    0.0333    0.1701    0.2473    0.3819    0.4201       NaN       NaN    0.4324    0.1503    1.7018
    1.1824    0.4887    0.0328    0.1697    0.2467    0.3802    0.4217       NaN       NaN    0.4195    0.1437    1.7256
    1.1845    0.4889    0.0323    0.1694    0.2461    0.3786    0.4233       NaN       NaN    0.4078    0.1377    1.7491
    1.1867    0.4891    0.0318    0.1691    0.2455    0.3769    0.4249       NaN       NaN    0.3970    0.1324    1.7726
    1.1888    0.4894    0.0313    0.1687    0.2449    0.3753    0.4266       NaN       NaN    0.3872    0.1275    1.7959
    1.1909    0.4898    0.0309    0.1684    0.2443    0.3737    0.4284       NaN       NaN    0.3782    0.1230    1.8192
    1.1930    0.4903    0.0305    0.1680    0.2437    0.3721    0.4303       NaN       NaN    0.3699    0.1190    1.8424
    1.1952    0.4908    0.0301    0.1677    0.2431    0.3706    0.4321       NaN       NaN    0.3622    0.1153    1.8655
    1.1973    0.4914    0.0297    0.1674    0.2425    0.3690    0.4341       NaN       NaN    0.3552    0.1119    1.8887
    1.1994    0.4920    0.0293    0.1670    0.2419    0.3675    0.4361       NaN       NaN    0.3487    0.1088    1.9119
    1.2015    0.4928    0.0290    0.1667    0.2413    0.3660    0.4382       NaN       NaN    0.3428    0.1059    1.9351
    1.2036    0.4936    0.0286    0.1664    0.2407    0.3645    0.4403       NaN       NaN    0.3373    0.1033    1.9584
    1.2058    0.4945    0.0283    0.1661    0.2402    0.3630    0.4425       NaN       NaN    0.3322    0.1008    1.9818
    1.2079    0.4954    0.0280    0.1657    0.2396    0.3615    0.4447       NaN       NaN    0.3276    0.0986    2.0053
    1.2100    0.4964    0.0277    0.1654    0.2390    0.3601    0.4470       NaN       NaN    0.3234    0.0965    2.0289
    1.2121    0.4975    0.0274    0.1651    0.2384    0.3586    0.4493       NaN       NaN    0.3195    0.0945    2.0527
    1.2143    0.4987    0.0271    0.1648    0.2379    0.3572    0.4517       NaN       NaN    0.3159    0.0927    2.0767
    1.2164    0.4999    0.0268    0.1644    0.2373    0.3558    0.4541       NaN       NaN    0.3127    0.0911    2.1009
    1.2185    0.5012    0.0265    0.1641    0.2367    0.3544    0.4566       NaN       NaN    0.3098    0.0896    2.1254
    1.2206    0.5026    0.0263    0.1638    0.2362    0.3530    0.4591       NaN       NaN    0.3072    0.0881    2.1502
    1.2228    0.5040    0.0260    0.1635    0.2356    0.3516    0.4617       NaN       NaN    0.3048    0.0868    2.1754
    1.2249    0.5055    0.0258    0.1632    0.2351    0.3503    0.4643       NaN       NaN    0.3027    0.0856    2.2010
    1.2270    0.5070    0.0255    0.1628    0.2345    0.3490    0.4669       NaN       NaN    0.3009    0.0845    2.2270
    1.2291    0.5086    0.0253    0.1625    0.2340    0.3476    0.4695       NaN       NaN    0.2994    0.0835    2.2536
    1.2312    0.5102    0.0251    0.1622    0.2334    0.3463    0.4722       NaN       NaN    0.2980    0.0825    2.2808
    1.2334    0.5118    0.0248    0.1619    0.2329    0.3450    0.4748       NaN       NaN    0.2970    0.0816    2.3087
    1.2355    0.5136    0.0246    0.1616    0.2323    0.3437    0.4775       NaN       NaN    0.2961    0.0809    2.3374
    1.2376    0.5153    0.0244    0.1613    0.2318    0.3425    0.4802       NaN       NaN    0.2955    0.0801    2.3671
    1.2397    0.5171    0.0242    0.1610    0.2313    0.3412    0.4829       NaN       NaN    0.2951    0.0795    2.3979
    1.2419    0.5188    0.0240    0.1607    0.2307    0.3399    0.4856       NaN       NaN    0.2950    0.0789    2.4299
    1.2440    0.5206    0.0238    0.1604    0.2302    0.3387    0.4883       NaN       NaN    0.2951    0.0784    2.4635
    1.2461    0.5224    0.0236    0.1600    0.2297    0.3375    0.4909       NaN       NaN    0.2954    0.0780    2.4985
    1.2482    0.5242    0.0234    0.1597    0.2292    0.3363    0.4936       NaN       NaN    0.2959    0.0776    2.5359
    1.2503    0.5260    0.0232    0.1594    0.2286    0.3350    0.4961       NaN       NaN    0.2967    0.0772    2.5759
    1.2525    0.5278    0.0231    0.1591    0.2281    0.3339    0.4987       NaN       NaN    0.2977    0.0770    2.6190
    1.2546    0.5296    0.0229    0.1588    0.2276    0.3327    0.5012       NaN       NaN    0.2989    0.0767    2.6660
    1.2567    0.5313    0.0227    0.1585    0.2271    0.3315    0.5036       NaN       NaN    0.3004    0.0766    2.7176
    1.2588    0.5330    0.0226    0.1582    0.2266    0.3303    0.5060       NaN       NaN    0.3022    0.0765    2.7752
    1.2610    0.5346    0.0224    0.1579    0.2261    0.3292    0.5082       NaN       NaN    0.3042    0.0764    2.8406
    1.2631    0.5362    0.0222    0.1576    0.2256    0.3280    0.5104       NaN       NaN    0.3064    0.0764    2.9162
    1.2652    0.5377    0.0221    0.1574    0.2251    0.3269    0.5125       NaN       NaN    0.3090    0.0765    3.0058
    1.2673    0.5391    0.0219    0.1571    0.2246    0.3258    0.5145       NaN       NaN    0.3118    0.0766    3.1156
    1.2695    0.5404    0.0218    0.1568    0.2241    0.3247    0.5163       NaN       NaN    0.3149    0.0768    3.2560
    1.2716    0.5417    0.0216    0.1565    0.2236    0.3236    0.5181       NaN       NaN    0.3184    0.0770    3.4480
    1.2737    0.5428    0.0215    0.1562    0.2231    0.3225    0.5197       NaN       NaN    0.3221    0.0773    3.7417
    1.2758    0.5439    0.0213    0.1559    0.2226    0.3214    0.5211       NaN       NaN    0.3263    0.0777    4.3195
    1.2779    0.5448    0.0212    0.1556    0.2221    0.3203    0.5225       NaN       NaN    0.3308    0.0781    4.4199
    1.2801    0.5502    0.0211    0.1553    0.2216    0.3193    0.5236    0.2359       NaN    0.3356    0.0786    2.7113
    1.2822    0.5550    0.0209    0.1550    0.2212    0.3182    0.5247    0.2772       NaN    0.3410    0.0791    2.1009
    1.2843    0.5590    0.0208    0.1547    0.2207    0.3172    0.5255    0.3018       NaN    0.3467    0.0797    1.7621
    1.2864    0.5623    0.0207    0.1545    0.2202    0.3161    0.5262    0.3193       NaN    0.3530    0.0804    1.5398
    1.2886    0.5652    0.0205    0.1542    0.2197    0.3151    0.5268    0.3328       NaN    0.3598    0.0812    1.3801
    1.2907    0.5676    0.0204    0.1539    0.2192    0.3141    0.5271    0.3437       NaN    0.3672    0.0820    1.2586
    1.2928    0.5696    0.0203    0.1536    0.2188    0.3131    0.5274    0.3529       NaN    0.3752    0.0829    1.1625
    1.2949    0.5714    0.0202    0.1533    0.2183    0.3121    0.5274    0.3606       NaN    0.3838    0.0839    1.0844
    1.2970    0.5728    0.0201    0.1531    0.2178    0.3111    0.5273    0.3673       NaN    0.3932    0.0850    1.0194
    1.2992    0.5740    0.0200    0.1528    0.2174    0.3101    0.5270    0.3731       NaN    0.4035    0.0861    0.9644
    1.3013    0.5749    0.0198    0.1525    0.2169    0.3091    0.5266    0.3781       NaN    0.4147    0.0874    0.9174
    1.3034    0.5755    0.0197    0.1522    0.2165    0.3082    0.5260    0.3826       NaN    0.4269    0.0888    0.8766
    1.3055    0.5759    0.0196    0.1519    0.2160    0.3072    0.5252    0.3864       NaN    0.4402    0.0903    0.8410
    1.3077    0.5761    0.0195    0.1517    0.2155    0.3062    0.5244    0.3898       NaN    0.4548    0.0919    0.8096
    1.3098    0.5761    0.0194    0.1514    0.2151    0.3053    0.5233    0.3927       NaN    0.4709    0.0937    0.7819
    1.3119    0.5759    0.0193    0.1511    0.2146    0.3044    0.5222    0.3953       NaN    0.4886    0.0956    0.7572
    1.3140    0.5755    0.0192    0.1509    0.2142    0.3034    0.5209    0.3975       NaN    0.5083    0.0977    0.7352
    1.3162    0.5749    0.0191    0.1506    0.2137    0.3025    0.5195    0.3994       NaN    0.5301    0.1000    0.7155
    1.3183    0.5742    0.0190    0.1503    0.2133    0.3016    0.5180    0.4010       NaN    0.5546    0.1024    0.6979
    1.3204    0.5734    0.0189    0.1500    0.2129    0.3007    0.5164    0.4023       NaN    0.5820    0.1051    0.6820
    1.3225    0.5724    0.0188    0.1498    0.2124    0.2998    0.5147    0.4035       NaN    0.6131    0.1080    0.6678
    1.3246    0.5712    0.0187    0.1495    0.2120    0.2989    0.5129    0.4044       NaN    0.6484    0.1112    0.6551
    1.3268    0.5700    0.0186    0.1492    0.2115    0.2980    0.5111    0.4051       NaN    0.6890    0.1147    0.6438
    1.3289    0.5686    0.0185    0.1490    0.2111    0.2971    0.5091    0.4057       NaN    0.7361    0.1186    0.6337
    1.3310    0.5672    0.0184    0.1487    0.2107    0.2962    0.5071    0.4061       NaN    0.7913    0.1228    0.6249
    1.3331    0.5656    0.0183    0.1485    0.2102    0.2953    0.5050    0.4063       NaN    0.8567    0.1276    0.6171
    1.3353    0.5640    0.0183    0.1482    0.2098    0.2945    0.5029    0.4065       NaN    0.9353    0.1329    0.6105
    1.3374    0.5623    0.0182    0.1479    0.2094    0.2936    0.5007    0.4065       NaN    1.0309    0.1389    0.6049
    1.3395    0.5606    0.0181    0.1477    0.2090    0.2928    0.4985    0.4064       NaN    1.1485    0.1456    0.6004
    1.3416    0.5588    0.0180    0.1474    0.2085    0.2919    0.4962    0.4063       NaN    1.2934    0.1534    0.5971
    1.3438    0.5570    0.0179    0.1472    0.2081    0.2911    0.4939    0.4060       NaN    1.4683    0.1625    0.5949
    1.3459    0.5551    0.0178    0.1469    0.2077    0.2903    0.4916    0.4057       NaN    1.6656    0.1731    0.5940
    1.3480    0.5532    0.0178    0.1466    0.2073    0.2894    0.4892    0.4054       NaN    1.8591    0.1855    0.5946
    1.3501    0.5512    0.0177    0.1464    0.2069    0.2886    0.4869    0.4050       NaN    2.0145    0.2001    0.5968
    1.3522    0.5492    0.0176    0.1461    0.2064    0.2878    0.4845    0.4046       NaN    2.1155    0.2168    0.6008
    1.3544    0.5473    0.0175    0.1459    0.2060    0.2870    0.4821    0.4041       NaN    2.1681    0.2356    0.6066
    1.3565    0.5453    0.0174    0.1456    0.2056    0.2862    0.4797    0.4036       NaN    2.1855    0.2562    0.6143
    1.3586    0.5433    0.0174    0.1454    0.2052    0.2854    0.4772    0.4031       NaN    2.1791    0.2781    0.6236
    1.3607    0.5413    0.0173    0.1451    0.2048    0.2846    0.4748    0.4025       NaN    2.1568    0.3012    0.6349
    1.3629    0.5393    0.0172    0.1449    0.2044    0.2838    0.4724    0.4020       NaN    2.1236    0.3250    0.6478
    1.3650    0.5373    0.0171    0.1446    0.2040    0.2830    0.4700    0.4014       NaN    2.0828    0.3494    0.6623
    1.3671    0.5353    0.0171    0.1444    0.2036    0.2822    0.4676    0.4009       NaN    2.0366    0.3741    0.6781
    1.3692    0.5333    0.0170    0.1441    0.2032    0.2815    0.4652    0.4003       NaN    1.9864    0.3988    0.6952
    1.3713    0.5314    0.0169    0.1439    0.2028    0.2807    0.4628    0.3998       NaN    1.9333    0.4234    0.7134
    1.3735    0.5294    0.0169    0.1436    0.2024    0.2799    0.4604    0.3993       NaN    1.8780    0.4478    0.7327
    1.3756    0.5275    0.0168    0.1434    0.2020    0.2792    0.4580    0.3987       NaN    1.8211    0.4717    0.7528
    1.3777    0.5256    0.0167    0.1431    0.2016    0.2784    0.4556    0.3982       NaN    1.7628    0.4950    0.7736
    1.3798    0.5237    0.0167    0.1429    0.2012    0.2777    0.4533    0.3977       NaN    1.7035    0.5175    0.7952
    1.3820    0.5219    0.0166    0.1427    0.2008    0.2770    0.4509    0.3973       NaN    1.6434    0.5389    0.8173
    1.3841    0.5201    0.0165    0.1424    0.2004    0.2762    0.4486    0.3969       NaN    1.5825    0.5590    0.8398
    1.3862    0.5183    0.0165    0.1422    0.2001    0.2755    0.4463    0.3964       NaN    1.5211    0.5774    0.8628
    1.3883    0.5165    0.0164    0.1419    0.1997    0.2748    0.4440    0.3961       NaN    1.4590    0.5939    0.8860
    1.3905    0.5148    0.0163    0.1417    0.1993    0.2740    0.4418    0.3957       NaN    1.3963    0.6078    0.9096
    1.3926    0.5131    0.0163    0.1415    0.1989    0.2733    0.4395    0.3954       NaN    1.3327    0.6184    0.9334
    1.3947    0.5114    0.0162    0.1412    0.1985    0.2726    0.4373    0.3951       NaN    1.2681    0.6250    0.9575
    1.3968    0.5098    0.0162    0.1410    0.1982    0.2719    0.4351    0.3948       NaN    1.2020    0.6261    0.9818
    1.3989    0.5082    0.0161    0.1407    0.1978    0.2712    0.4329    0.3946       NaN    1.1337    0.6199    1.0068
    1.4011    0.5066    0.0160    0.1405    0.1974    0.2705    0.4307    0.3944       NaN    1.0622    0.6037    1.0327
    1.4032    0.5051    0.0160    0.1403    0.1970    0.2698    0.4286    0.3943       NaN    0.9857    0.5730    1.0607
    1.4053    0.5036    0.0159    0.1400    0.1967    0.2691    0.4265    0.3942       NaN    0.9003    0.5189    1.0964
    1.4074    0.5022    0.0159    0.1398    0.1963    0.2684    0.4244    0.3941       NaN    0.8109    0.4449    1.1687
    1.4096    0.5008    0.0158    0.1396    0.1959    0.2678    0.4223    0.3941       NaN    0.7364    0.3818    1.2168
    1.4117    0.4995    0.0157    0.1393    0.1956    0.2671    0.4202    0.3941       NaN    0.6771    0.3331    1.2533
    1.4138    0.4981    0.0157    0.1391    0.1952    0.2664    0.4182    0.3941       NaN    0.6295    0.2955    1.2833
    1.4159    0.4969    0.0156    0.1389    0.1948    0.2657    0.4162    0.3942       NaN    0.5903    0.2660    1.3095
    1.4181    0.4956    0.0156    0.1386    0.1945    0.2651    0.4142    0.3944       NaN    0.5575    0.2422    1.3333
    1.4202    0.4945    0.0155    0.1384    0.1941    0.2644    0.4122    0.3946       NaN    0.5294    0.2227    1.3556
    1.4223    0.4933    0.0155    0.1382    0.1937    0.2638    0.4102    0.3948       NaN    0.5052    0.2064    1.3768
    1.4244    0.4922    0.0154    0.1380    0.1934    0.2631    0.4083    0.3951       NaN    0.4839    0.1925    1.3973
    1.4265    0.4912    0.0154    0.1377    0.1930    0.2625    0.4064    0.3954       NaN    0.4652    0.1806    1.4171
    1.4287    0.4902    0.0153    0.1375    0.1927    0.2618    0.4045    0.3958       NaN    0.4485    0.1703    1.4366
    1.4308    0.4892    0.0153    0.1373    0.1923    0.2612    0.4026    0.3962       NaN    0.4335    0.1613    1.4556
    1.4329    0.4883    0.0152    0.1371    0.1920    0.2605    0.4008    0.3966       NaN    0.4200    0.1534    1.4744
    1.4350    0.4875    0.0152    0.1368    0.1916    0.2599    0.3989    0.3971       NaN    0.4078    0.1463    1.4930
    1.4372    0.4866    0.0151    0.1366    0.1913    0.2593    0.3971    0.3977       NaN    0.3968    0.1400    1.5114
    1.4393    0.4859    0.0151    0.1364    0.1909    0.2587    0.3953    0.3983       NaN    0.3867    0.1343    1.5296
    1.4414    0.4852    0.0150    0.1362    0.1906    0.2580    0.3936    0.3989       NaN    0.3775    0.1292    1.5477
    1.4435    0.4845    0.0150    0.1359    0.1902    0.2574    0.3918    0.3996       NaN    0.3691    0.1245    1.5657
    1.4456    0.4839    0.0149    0.1357    0.1899    0.2568    0.3901    0.4004       NaN    0.3614    0.1203    1.5837
    1.4478    0.4833    0.0149    0.1355    0.1895    0.2562    0.3884    0.4011       NaN    0.3544    0.1165    1.6016
    1.4499    0.4828    0.0148    0.1353    0.1892    0.2556    0.3867    0.4020       NaN    0.3479    0.1129    1.6194
    1.4520    0.4824    0.0148    0.1351    0.1888    0.2550    0.3850    0.4029       NaN    0.3420    0.1097    1.6373
    1.4541    0.4820    0.0147    0.1349    0.1885    0.2544    0.3833    0.4038       NaN    0.3366    0.1067    1.6551
    1.4563    0.4816    0.0147    0.1346    0.1882    0.2538    0.3817    0.4048       NaN    0.3316    0.1040    1.6729
    1.4584    0.4813    0.0146    0.1344    0.1878    0.2532    0.3801    0.4058       NaN    0.3270    0.1015    1.6907
    1.4605    0.4811    0.0146    0.1342    0.1875    0.2526    0.3784    0.4069       NaN    0.3229    0.0992    1.7086
    1.4626    0.4809    0.0145    0.1340    0.1872    0.2520    0.3769    0.4081       NaN    0.3191    0.0970    1.7265
    1.4648    0.4808    0.0145    0.1338    0.1868    0.2514    0.3753    0.4092       NaN    0.3156    0.0951    1.7444
    1.4669    0.4807    0.0145    0.1336    0.1865    0.2509    0.3737    0.4105       NaN    0.3125    0.0932    1.7624
    1.4690    0.4807    0.0144    0.1333    0.1862    0.2503    0.3722    0.4118       NaN    0.3097    0.0915    1.7805
    1.4711    0.4807    0.0144    0.1331    0.1858    0.2497    0.3707    0.4131       NaN    0.3072    0.0900    1.7987
    1.4732    0.4808    0.0143    0.1329    0.1855    0.2491    0.3692    0.4145       NaN    0.3050    0.0885    1.8169
    1.4754    0.4809    0.0143    0.1327    0.1852    0.2486    0.3677    0.4159       NaN    0.3031    0.0872    1.8353
    1.4775    0.4811    0.0142    0.1325    0.1849    0.2480    0.3662    0.4174       NaN    0.3015    0.0859    1.8538
    1.4796    0.4814    0.0142    0.1323    0.1845    0.2474    0.3647    0.4189       NaN    0.3001    0.0848    1.8725
    1.4817    0.4817    0.0142    0.1321    0.1842    0.2469    0.3633    0.4205       NaN    0.2989    0.0838    1.8914
    1.4839    0.4821    0.0141    0.1319    0.1839    0.2463    0.3619    0.4221       NaN    0.2980    0.0828    1.9104
    1.4860    0.4825    0.0141    0.1317    0.1836    0.2458    0.3604    0.4237       NaN    0.2974    0.0819    1.9297
    1.4881    0.4830    0.0140    0.1315    0.1833    0.2452    0.3590    0.4254       NaN    0.2969    0.0811    1.9492
    1.4902    0.4835    0.0140    0.1313    0.1829    0.2447    0.3576    0.4272       NaN    0.2968    0.0804    1.9689
    1.4923    0.4841    0.0140    0.1310    0.1826    0.2441    0.3563    0.4290       NaN    0.2968    0.0798    1.9889
    1.4945    0.4848    0.0139    0.1308    0.1823    0.2436    0.3549    0.4308       NaN    0.2971    0.0792    2.0092
    1.4966    0.4855    0.0139    0.1306    0.1820    0.2431    0.3536    0.4327       NaN    0.2977    0.0787    2.0299
    1.4987    0.4862    0.0139    0.1304    0.1817    0.2425    0.3522    0.4346       NaN    0.2984    0.0782    2.0509
    1.5008    0.4870    0.0138    0.1302    0.1814    0.2420    0.3509    0.4365       NaN    0.2995    0.0778    2.0724
    1.5030    0.4879    0.0138    0.1300    0.1810    0.2415    0.3496    0.4385       NaN    0.3007    0.0775    2.0943
    1.5051    0.4888    0.0137    0.1298    0.1807    0.2409    0.3483    0.4405       NaN    0.3023    0.0772    2.1168
    1.5072    0.4897    0.0137    0.1296    0.1804    0.2404    0.3470    0.4425       NaN    0.3040    0.0770    2.1398
    1.5093    0.4907    0.0137    0.1294    0.1801    0.2399    0.3458    0.4445       NaN    0.3061    0.0768    2.1634
    1.5115    0.4918    0.0136    0.1292    0.1798    0.2394    0.3445    0.4466       NaN    0.3084    0.0768    2.1878
    1.5136    0.4928    0.0136    0.1290    0.1795    0.2389    0.3433    0.4487       NaN    0.3110    0.0767    2.2129
    1.5157    0.4940    0.0136    0.1288    0.1792    0.2384    0.3420    0.4508       NaN    0.3140    0.0767    2.2390
    1.5178    0.4951    0.0135    0.1286    0.1789    0.2378    0.3408    0.4530       NaN    0.3172    0.0768    2.2662
    1.5199    0.4963    0.0135    0.1284    0.1786    0.2373    0.3396    0.4551       NaN    0.3208    0.0770    2.2945
    1.5221    0.4975    0.0134    0.1282    0.1783    0.2368    0.3384    0.4572       NaN    0.3247    0.0771    2.3242
    1.5242    0.4987    0.0134    0.1280    0.1780    0.2363    0.3372    0.4594       NaN    0.3290    0.0774    2.3554
    1.5263    0.5000    0.0134    0.1278    0.1777    0.2358    0.3360    0.4615       NaN    0.3337    0.0777    2.3885
    1.5284    0.5012    0.0133    0.1276    0.1774    0.2353    0.3349    0.4637       NaN    0.3388    0.0781    2.4238
    1.5306    0.5025    0.0133    0.1274    0.1771    0.2348    0.3337    0.4658       NaN    0.3444    0.0785    2.4616
    1.5327    0.5038    0.0133    0.1272    0.1768    0.2343    0.3326    0.4679       NaN    0.3504    0.0790    2.5024
    1.5348    0.5050    0.0132    0.1271    0.1765    0.2339    0.3314    0.4699       NaN    0.3570    0.0796    2.5469
    1.5369    0.5063    0.0132    0.1269    0.1762    0.2334    0.3303    0.4720       NaN    0.3641    0.0802    2.5959
    1.5391    0.5076    0.0132    0.1267    0.1759    0.2329    0.3292    0.4740       NaN    0.3719    0.0810    2.6505
    1.5412    0.5088    0.0131    0.1265    0.1756    0.2324    0.3281    0.4759       NaN    0.3803    0.0817    2.7123
    1.5433    0.5100    0.0131    0.1263    0.1753    0.2319    0.3270    0.4778       NaN    0.3895    0.0826    2.7834
    1.5454    0.5112    0.0131    0.1261    0.1751    0.2314    0.3259    0.4797       NaN    0.3995    0.0836    2.8670
    1.5475    0.5124    0.0130    0.1259    0.1748    0.2310    0.3249    0.4815       NaN    0.4104    0.0846    2.9681
    1.5497    0.5135    0.0130    0.1257    0.1745    0.2305    0.3238    0.4832       NaN    0.4224    0.0858    3.0951
    1.5518    0.5145    0.0130    0.1255    0.1742    0.2300    0.3227    0.4849       NaN    0.4355    0.0870    3.2640
    1.5539    0.5156    0.0130    0.1253    0.1739    0.2296    0.3217    0.4865       NaN    0.4498    0.0884    3.5100
    1.5560    0.5165    0.0129    0.1251    0.1736    0.2291    0.3207    0.4880       NaN    0.4657    0.0899    3.9420
    1.5582    0.5174    0.0129    0.1250    0.1733    0.2286    0.3196    0.4894       NaN    0.4832    0.0915    5.7283
    1.5603    0.5182    0.0129    0.1248    0.1731    0.2282    0.3186    0.4907       NaN    0.5026    0.0932    3.0283
    1.5624    0.5261    0.0128    0.1246    0.1728    0.2277    0.3176    0.4919    0.2524    0.5243    0.0951    2.1834
    1.5645    0.5301    0.0128    0.1244    0.1725    0.2272    0.3166    0.4930    0.2790    0.5486    0.0972    1.7823
    1.5666    0.5335    0.0128    0.1242    0.1722    0.2268    0.3156    0.4940    0.2970    0.5760    0.0994    1.5362
    1.5688    0.5364    0.0127    0.1240    0.1719    0.2263    0.3146    0.4948    0.3107    0.6071    0.1019    1.3663
    1.5709    0.5390    0.0127    0.1238    0.1717    0.2259    0.3136    0.4956    0.3218    0.6427    0.1046    1.2405
    1.5730    0.5413    0.0127    0.1237    0.1714    0.2254    0.3127    0.4962    0.3310    0.6837    0.1076    1.1430
    1.5751    0.5434    0.0127    0.1235    0.1711    0.2250    0.3117    0.4968    0.3389    0.7316    0.1108    1.0649
    1.5773    0.5452    0.0126    0.1233    0.1708    0.2245    0.3107    0.4972    0.3458    0.7880    0.1145    1.0009
    1.5794    0.5468    0.0126    0.1231    0.1706    0.2241    0.3098    0.4974    0.3519    0.8553    0.1185    0.9475
    1.5815    0.5481    0.0126    0.1229    0.1703    0.2237    0.3089    0.4976    0.3573    0.9368    0.1230    0.9023
    1.5836    0.5493    0.0125    0.1227    0.1700    0.2232    0.3079    0.4976    0.3621    1.0366    0.1280    0.8638
    1.5858    0.5503    0.0125    0.1226    0.1697    0.2228    0.3070    0.4975    0.3664    1.1596    0.1338    0.8294
    1.5879    0.5511    0.0125    0.1224    0.1695    0.2224    0.3061    0.4973    0.3703    1.3102    0.1404    0.8010
    1.5900         0       NaN       NaN       NaN       NaN       NaN       NaN       NaN    1.4865    0.1481    0.7765];

% proceed by phase
switch lower(ph)
    case {'r' 'rayl' 'rayleigh'}
        if(nargout>1)
            c=interp1(xc(:,1),xc(:,3:9),x,'spline',0);
        else
            c=interp1(xc(:,1),xc(:,2),x,'spline',0);
        end
    case 'p'
        c=interp1(xc(:,1),xc(:,10),x,'spline',0);
    case 's'
        c=interp1(xc(:,1),xc(:,11),x,'spline',0);
    case {'sa' 'sall'}
        c=interp1(xc(:,1),xc(:,12),x,'spline',0);
    otherwise
        error('seizmo:bathy_micro_excite:badInput',...
            'PH not a recognized type!');
end

% set nans and negative values to 0
c(c<0 | isnan(c))=0;

% nargout
if(nargout>1)
    varargout={c(:,1) c(:,2) c(:,3) c(:,4) c(:,5) c(:,6) c(:,7)};
else
    varargout={c};
end

end
