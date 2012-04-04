function [mout]=prem(varargin)
%PREM    Returns the isotropic PREM Earth Model at a reference period of 1s
%
%    Usage:    model=prem()
%              model=prem(...,'depths',depths,...)
%              model=prem(...,'dcbelow',false,...)
%              model=prem(...,'range',[top bottom],...)
%              model=prem(...,'crust',true|false,...)
%              model=prem(...,'ocean',true|false,...)
%
%    Description:
%     MODEL=PREM() returns a struct containing the isotropic version of the
%     1D radial Earth model PREM.  Note that by default the ocean has been
%     replaced with an extended upper crust.  The struct has the following
%     fields:
%      MODEL.name      - model name ('PREM')
%           .ocean     - true if OCEAN & CRUST are both TRUE
%           .crust     - true if CRUST is TRUE
%           .isotropic - always true here
%           .refperiod - always 1sec here
%           .flattened - always false here (see FLATTEN_1DMODEL)
%           .depth     - km depths from 0 to 6371
%           .vp        - isotropic p-wave velocity at 1s periods (km/s)
%           .vs        - isotropic s-wave velocity at 1s periods (km/s)
%           .rho       - density (g/cm^3)
%           .qk        - bulk moduli quality factor
%           .qu        - shear moduli quality factor
%     Note that the model includes repeated depths at discontinuities.
%
%     MODEL=PREM(...,'DEPTHS',DEPTHS,...) returns the model parameters
%     only at the depths in DEPTHS.  DEPTHS is assumed to be in km.  The
%     model parameters are found by linear interpolation between known
%     values.  DEPTHS at discontinuities return values from the deeper
%     (bottom) side of the discontinuity for the first time and from the
%     top side for the second time.  Depths can not be repeated more than
%     twice and must be monotonically non-decreasing.
%
%     MODEL=PREM(...,'DCBELOW',FALSE,...) returns values from the
%     shallow (top) side of the discontinuity the first time a depth is
%     given at one (using the DEPTHS option) if DCBELOW is FALSE.  The
%     default is TRUE (returns value from bottom-side the first time).  The
%     second time a depth is used, the opposite side is given.
%
%     MODEL=PREM(...,'RANGE',[TOP BOTTOM],...) specifies the range of
%     depths that known model parameters are returned.  [TOP BOTTOM] must
%     be a 2 element array in km.  Note this does not block depths given by
%     the DEPTHS option.
%
%     MODEL=PREM(...,'CRUST',TRUE|FALSE,...) indicates if the crust of
%     PREM is to be removed or not.  Setting CRUST to FALSE will return a
%     crustless model (the mantle is extended to the surface using linear
%     interpolation).  This will also remove the ocean.
%
%     MODEL=PREM(...,'OCEAN',TRUE|FALSE,...) indicates if the ocean of
%     PREM is to be removed or not.  Setting CRUST to FALSE will return a
%     oceanless model (the crust (or mantle when CRUST is FALSE) is
%     extended to the surface using linear interpolation).  The default is
%     FALSE (no ocean).
%
%    Notes:
%     - Note that this is just a linear interpolation of PREM isotropic at
%       1 second reference period
%     - PREM reference:
%        Dziewonski & Anderson 1981, Preliminary reference earth model,
%        Phys. Earth planet. Inter. 25, pp. 297-356
%
%    Examples:
%     % Plot parameters for the CMB region:
%     model=prem('r',[2600 3400]);
%     figure;
%     plot(model.depth,model.vp,'r',...
%          model.depth,model.vs,'g',...
%          model.depth,model.rho,'b','linewidth',2);
%     title('PREM')
%     legend({'Vp' 'Vs' '\rho'})
%
%    See also: AK135, IASP91

%     Version History:
%        May  19, 2010 - initial version
%        May  20, 2010 - discon on edge handling, quicker
%        May  23, 2010 - use output from my prem_perfect routine rather
%                        than using values from taup (bugged!), qk/qu
%                        rather than qp/qs
%        May  24, 2010 - added several struct fields for info
%        Aug.  8, 2010 - minor doc touch, dcbelow option
%        Aug.  9, 2010 - fix Qu to be 1e10 in the ocean & outercore
%        Aug. 17, 2010 - added reference
%        Sep. 19, 2010 - back to Inf Qu in ocean & outercore, doc update,
%                        better discontinuity support
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 14:45 GMT

% todo:

% check nargin
if(mod(nargin,2))
    error('seizmo:prem:badNumInputs',...
        'Unpaired Option/Value!');
end

% option defaults
varargin=[{'d' [] 'b' true 'c' true 'o' false 'r' [0 6371]} varargin];

% check options
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:prem:badOption',...
        'All Options must be specified with a string!');
end
for i=1:2:numel(varargin)
    % skip empty
    skip=false;
    if(isempty(varargin{i+1})); skip=true; end

    % check option is available
    switch lower(varargin{i})
        case {'o' 'ocean'}
            if(skip); continue; end
            if(~islogical(varargin{i+1}) || ~isscalar(varargin{i+1}))
                error('seizmo:prem:badOCEAN',...
                    'OCEAN must be a TRUE or FALSE!');
            end
            ocean=varargin{i+1};
        case {'d' 'dep' 'depth' 'depths'}
            if(~isempty(varargin{i+1}))
                if(~isreal(varargin{i+1}) || any(varargin{i+1}<0 ...
                        | varargin{i+1}>6371) || any(isnan(varargin{i+1})))
                    error('seizmo:prem:badDEPTHS',...
                        ['DEPTHS must be real-valued km depths within ' ...
                        'the range [0 6371] in km!']);
                elseif(any(diff(varargin{i+1})<0))
                    error('seizmo:prem:badDEPTHS',...
                        'DEPTHS must be monotonically non-increasing!');
                elseif(any(histc(varargin{i+1},...
                        varargin{i+1}([find(diff(varargin{i+1}));end]))>3))
                    error('seizmo:prem:badDEPTHS',...
                        'DEPTHS has values repeated 3+ times!');
                end
            end
            depths=varargin{i+1}(:);
        case {'dcb' 'dc' 'below' 'b' 'dcbelow'}
            if(skip); continue; end
            if(~islogical(varargin{i+1}) || ~isscalar(varargin{i+1}))
                error('seizmo:prem:badDCBELOW',...
                    'DCBELOW must be a TRUE or FALSE!');
            end
            dcbelow=varargin{i+1};
        case {'c' 'cru' 'crust'}
            if(skip); continue; end
            if(~islogical(varargin{i+1}) || ~isscalar(varargin{i+1}))
                error('seizmo:prem:badCRUST',...
                    'CRUST must be a TRUE or FALSE!');
            end
            crust=varargin{i+1};
        case {'r' 'rng' 'range'}
            if(skip); continue; end
            if(~isreal(varargin{i+1}) || numel(varargin{i+1})~=2)
                error('seizmo:prem:badRANGE',...
                    ['RANGE must be a 2 element vector specifying ' ...
                    '[TOP BOTTOM] in km!']);
            end
            range=sort(varargin{i+1});
        otherwise
            error('seizmo:prem:badOption',...
                'Unknown Option: %s',varargin{i});
    end
end

% the prem model
% - extracted from my perfect prem routine :P
model=[
      0      1.45         0      1.02     57823       inf
      3      1.45         0      1.02     57823       inf
      3       5.8       3.2       2.6     57823       600
     15       5.8       3.2       2.6     57823       600
     15       6.8       3.9       2.9     57823       600
   24.4       6.8       3.9       2.9     57823       600
   24.4   8.11062   4.49101   3.38075     57823       600
 46.855   8.09674   4.48273   3.37831     57823       600
69.2687   8.08288   4.47447   3.37587     57823       600
     80   8.07625   4.47052   3.37471     57823       600
     80   8.07625   4.47052   3.37471     57823        80
102.353   8.06243   4.46228   3.37228     57823        80
124.664   8.04864   4.45405   3.36985     57823        80
146.934   8.03487   4.44585   3.36743     57823        80
169.163   8.02113   4.43765   3.36502     57823        80
191.352   8.00742   4.42948    3.3626     57823        80
213.499   7.99373   4.42131    3.3602     57823        80
    220   7.98971   4.41892   3.35949     57823        80
    220   8.55895    4.6439   3.43577     57823       143
 243.22   8.60362   4.66015   3.44964     57823       143
 266.52   8.64845   4.67646   3.46355     57823       143
289.903   8.69343   4.69283   3.47752     57823       143
313.367   8.73857   4.70926   3.49153     57823       143
336.913   8.78387   4.72574   3.50559     57823       143
360.542   8.82933   4.74228    3.5197     57823       143
384.253   8.87495   4.75888   3.53386     57823       143
    400   8.90524    4.7699   3.54326     57823       143
    400   9.13392   4.93249   3.72375     57823       143
424.662   9.26018   5.00443   3.75483     57823       143
449.685   9.38828   5.07743   3.78637     57823       143
475.072   9.51825   5.15149   3.81836     57823       143
500.829   9.65012   5.22663   3.85083     57823       143
526.962   9.78391   5.30286   3.88377     57823       143
553.477   9.91965   5.38021   3.91718     57823       143
580.378   10.0574   5.45869   3.95109     57823       143
    600   10.1578   5.51593   3.97582     57823       143
    600   10.1578   5.51602   3.97582     57823       143
 627.58   10.2005   5.53737   3.98224     57823       143
655.267   10.2434    5.5588   3.98869     57823       143
    670   10.2662   5.57021   3.99212     57823       143
    670   10.7513   5.94513   4.38074     57823       312
699.726   10.8439   6.03204   4.39927     57823       312
729.886   10.9377   6.12025   4.41794     57823       312
760.487   11.0329   6.20978   4.43676     57823       312
    771   11.0656   6.24054    4.4432     57823       312
    771   11.0656   6.24039    4.4432     57823       312
802.202   11.1225   6.26274   4.46223     57823       312
833.516   11.1787   6.28483   4.48121     57823       312
 864.94   11.2343   6.30667   4.50013     57823       312
896.473   11.2892   6.32827   4.51899     57823       312
928.114   11.3434   6.34962    4.5378     57823       312
959.863    11.397   6.37073   4.55656     57823       312
991.716     11.45    6.3916   4.57526     57823       312
1023.67   11.5023   6.41224   4.59391     57823       312
1055.74    11.554   6.43264   4.61251     57823       312
 1087.9    11.605   6.45282   4.63106     57823       312
1120.16   11.6555   6.47277   4.64955     57823       312
1152.53   11.7054    6.4925   4.66799     57823       312
1184.99   11.7548   6.51202   4.68639     57823       312
1217.55   11.8035   6.53132   4.70474     57823       312
1250.21   11.8518   6.55041   4.72303     57823       312
1282.96   11.8994    6.5693   4.74128     57823       312
 1315.8   11.9466     6.588   4.75949     57823       312
1348.74   11.9933   6.60649   4.77765     57823       312
1381.78   12.0395    6.6248   4.79577     57823       312
 1414.9   12.0852   6.64293   4.81384     57823       312
1448.12   12.1304   6.66087   4.83188     57823       312
1481.42   12.1752   6.67865   4.84987     57823       312
1514.81   12.2196   6.69626   4.86783     57823       312
1548.29   12.2636   6.71371   4.88575     57823       312
1581.86   12.3072     6.731   4.90364     57823       312
1615.52   12.3505   6.74814   4.92149     57823       312
1649.26   12.3934   6.76514   4.93932     57823       312
1683.08   12.4359   6.78201   4.95711     57823       312
1716.99   12.4782   6.79875   4.97487     57823       312
1750.99   12.5201   6.81536   4.99261     57823       312
1785.06   12.5618   6.83186   5.01033     57823       312
1819.22   12.6033   6.84825   5.02802     57823       312
1853.47   12.6445   6.86454    5.0457     57823       312
1887.79   12.6855   6.88074   5.06335     57823       312
1922.19   12.7264   6.89685     5.081     57823       312
1956.68   12.7671   6.91289   5.09863     57823       312
1991.24   12.8076   6.92885   5.11624     57823       312
2025.88   12.8481   6.94475   5.13385     57823       312
2060.61   12.8884    6.9606   5.15146     57823       312
2095.41   12.9287    6.9764   5.16906     57823       312
2130.29    12.969   6.99216   5.18666     57823       312
2165.25   13.0092    7.0079   5.20427     57823       312
2200.29   13.0495   7.02362   5.22188     57823       312
2235.41   13.0898   7.03932   5.23949     57823       312
2270.61   13.1301   7.05502   5.25712     57823       312
2305.88   13.1706   7.07073   5.27476     57823       312
2341.24   13.2111   7.08646   5.29242     57823       312
2376.67   13.2519   7.10221    5.3101     57823       312
2412.18   13.2928     7.118    5.3278     57823       312
2447.77   13.3339   7.13383   5.34552     57823       312
2483.44   13.3752   7.14972   5.36328     57823       312
2519.19   13.4168   7.16567   5.38106     57823       312
2555.02   13.4587    7.1817   5.39889     57823       312
2590.93   13.5009   7.19781   5.41675     57823       312
2626.91   13.5435   7.21403   5.43466     57823       312
2662.98   13.5864   7.23034   5.45261     57823       312
2699.14   13.6298   7.24678   5.47061     57823       312
2735.37   13.6736   7.26334   5.48867     57823       312
   2741   13.6804   7.26593   5.49148     57823       312
   2741   13.6804   7.26593   5.49148     57823       312
2777.33   13.6891   7.26565    5.5096     57823       312
2813.66   13.6978   7.26535   5.52774     57823       312
2849.98   13.7066   7.26502   5.54591     57823       312
2886.31   13.7155   7.26465    5.5641     57823       312
   2891   13.7166    7.2646   5.56646     57823       312
   2891   8.06479         0   9.90045     57823       inf
2931.32   8.13333         0   9.96452     57823       inf
2971.99   8.20106         0   10.0281     57823       inf
   3013   8.26796         0   10.0912     57823       inf
3054.34     8.334         0   10.1538     57823       inf
3096.01   8.39918         0   10.2158     57823       inf
   3138   8.46347         0   10.2772     57823       inf
3180.32   8.52686         0   10.3381     57823       inf
3222.95   8.58933         0   10.3983     57823       inf
 3265.9   8.65087         0   10.4579     57823       inf
3309.15   8.71148         0   10.5169     57823       inf
3352.71   8.77113         0   10.5752     57823       inf
3396.57   8.82982         0   10.6328     57823       inf
3440.72   8.88754         0   10.6897     57823       inf
3485.15   8.94429         0   10.7458     57823       inf
3529.88   9.00006         0   10.8013     57823       inf
3574.88   9.05485         0    10.856     57823       inf
3620.15   9.10866         0   10.9099     57823       inf
3665.69   9.16149         0    10.963     57823       inf
 3711.5   9.21334         0   11.0153     57823       inf
3757.57    9.2642         0   11.0669     57823       inf
3803.89    9.3141         0   11.1176     57823       inf
3850.46   9.36302         0   11.1675     57823       inf
3897.27   9.41099         0   11.2165     57823       inf
3944.33     9.458         0   11.2647     57823       inf
3991.62   9.50407         0   11.3121     57823       inf
4039.14   9.54921         0   11.3586     57823       inf
4086.89   9.59343         0   11.4042     57823       inf
4134.85   9.63675         0   11.4489     57823       inf
4183.04   9.67919         0   11.4928     57823       inf
4231.43   9.72075         0   11.5358     57823       inf
4280.04   9.76146         0   11.5778     57823       inf
4328.84   9.80133         0   11.6191     57823       inf
4377.85    9.8404         0   11.6594     57823       inf
4427.05   9.87868         0   11.6988     57823       inf
4476.45   9.91619         0   11.7373     57823       inf
4526.03   9.95296         0   11.7749     57823       inf
4575.79   9.98901         0   11.8117     57823       inf
4625.74   10.0244         0   11.8476     57823       inf
4675.86   10.0591         0   11.8825     57823       inf
4726.15   10.0932         0   11.9166     57823       inf
4776.62   10.1266         0   11.9498     57823       inf
4827.25   10.1596         0   11.9822     57823       inf
4878.05   10.1919         0   12.0136     57823       inf
4929.01   10.2238         0   12.0442     57823       inf
4980.13   10.2552         0    12.074     57823       inf
5031.41   10.2862         0   12.1029     57823       inf
5082.84   10.3167         0    12.131     57823       inf
5134.42    10.347         0   12.1582     57823       inf
 5149.5   10.3557         0    12.166     57823       inf
 5149.5   11.0283   3.50431   12.7636    1327.7      84.6
5167.02   11.0349   3.50897   12.7729    1327.7      84.6
5184.57   11.0415   3.51356    12.782    1327.7      84.6
5202.13    11.048    3.5181    12.791    1327.7      84.6
5219.72   11.0544   3.52257   12.7999    1327.7      84.6
5237.34   11.0607   3.52698   12.8087    1327.7      84.6
5254.97   11.0669   3.53133   12.8173    1327.7      84.6
5272.63    11.073   3.53561   12.8258    1327.7      84.6
5290.31   11.0791   3.53983   12.8342    1327.7      84.6
5308.01    11.085   3.54399   12.8425    1327.7      84.6
5325.73   11.0909   3.54808   12.8506    1327.7      84.6
5343.47   11.0967   3.55211   12.8586    1327.7      84.6
5361.23   11.1023   3.55608   12.8665    1327.7      84.6
5379.01   11.1079   3.55998   12.8742    1327.7      84.6
5396.81   11.1134   3.56381   12.8819    1327.7      84.6
5414.63   11.1188   3.56758   12.8893    1327.7      84.6
5432.46   11.1241   3.57128   12.8967    1327.7      84.6
5450.32   11.1293   3.57492   12.9039    1327.7      84.6
 5468.2   11.1344   3.57849    12.911    1327.7      84.6
5486.09   11.1394     3.582    12.918    1327.7      84.6
   5504   11.1443   3.58544   12.9248    1327.7      84.6
5521.93   11.1492   3.58881   12.9315    1327.7      84.6
5539.87   11.1539   3.59211   12.9381    1327.7      84.6
5557.83   11.1585   3.59535   12.9445    1327.7      84.6
5575.81   11.1631   3.59851   12.9508    1327.7      84.6
 5593.8   11.1675   3.60161    12.957    1327.7      84.6
5611.81   11.1718   3.60465    12.963    1327.7      84.6
5629.83   11.1761   3.60761   12.9689    1327.7      84.6
5647.87   11.1802    3.6105   12.9746    1327.7      84.6
5665.92   11.1843   3.61333   12.9803    1327.7      84.6
5683.99   11.1882   3.61608   12.9857    1327.7      84.6
5702.07    11.192   3.61877   12.9911    1327.7      84.6
5720.16   11.1958   3.62139   12.9963    1327.7      84.6
5738.27   11.1994   3.62393   13.0013    1327.7      84.6
5756.39    11.203   3.62641   13.0062    1327.7      84.6
5774.52   11.2064   3.62882    13.011    1327.7      84.6
5792.66   11.2098   3.63115   13.0157    1327.7      84.6
5810.82    11.213   3.63342   13.0202    1327.7      84.6
5828.99   11.2161   3.63561   13.0245    1327.7      84.6
5847.17   11.2192   3.63773   13.0288    1327.7      84.6
5865.35   11.2221   3.63978   13.0328    1327.7      84.6
5883.55   11.2249   3.64177   13.0368    1327.7      84.6
5901.76   11.2277   3.64367   13.0406    1327.7      84.6
5919.98   11.2303   3.64551   13.0442    1327.7      84.6
5938.21   11.2328   3.64728   13.0477    1327.7      84.6
5956.44   11.2353   3.64897   13.0511    1327.7      84.6
5974.69   11.2376   3.65059   13.0543    1327.7      84.6
5992.94   11.2398   3.65214   13.0574    1327.7      84.6
 6011.2   11.2419   3.65362   13.0603    1327.7      84.6
6029.47   11.2439   3.65502   13.0631    1327.7      84.6
6047.75   11.2458   3.65635   13.0657    1327.7      84.6
6066.03   11.2476   3.65761   13.0682    1327.7      84.6
6084.32   11.2493   3.65879   13.0706    1327.7      84.6
6102.61   11.2509   3.65991   13.0728    1327.7      84.6
6120.91   11.2524   3.66095   13.0749    1327.7      84.6
6139.21   11.2538   3.66191   13.0768    1327.7      84.6
6157.52   11.2551   3.66281   13.0786    1327.7      84.6
6175.84   11.2562   3.66363   13.0802    1327.7      84.6
6194.16   11.2573   3.66437   13.0817    1327.7      84.6
6212.48   11.2583   3.66505    13.083    1327.7      84.6
 6230.8   11.2591   3.66565   13.0842    1327.7      84.6
6249.13   11.2599   3.66617   13.0853    1327.7      84.6
6267.46   11.2605   3.66663   13.0862    1327.7      84.6
6285.79   11.2611     3.667   13.0869    1327.7      84.6
6304.13   11.2615   3.66731   13.0875    1327.7      84.6
6322.47   11.2618   3.66754    13.088    1327.7      84.6
 6340.8   11.2621    3.6677   13.0883    1327.7      84.6
6359.14   11.2622   3.66778   13.0885    1327.7      84.6
   6371   11.2622    3.6678   13.0885    1327.7      84.6];

% remove ocean if desired
if(~ocean)
    % linearly extrapolated to the surface
    model(1,:)=[      0       5.8       3.2       2.6     57823       600];
    model(2:3,:)=[];
end

% remove crust & ocean if desired
if(~crust)
    % linearly extrapolated to the surface
    model(1,:)=[      0    8.1257       4.5    3.3834     57823       600];
    model(2:(4+2*ocean),:)=[];
end

% replace infinity with a high value
% - this is mainly for so interpolation doesn't produce NaNs
%model(isinf(model))=1e10;

% interpolate depths if desired
if(~isempty(depths))
    %depths=depths(depths>=range(1) & depths<=range(2));
    [bot,top]=interpdc1(model(:,1),model(:,2:end),depths);
    if(dcbelow)
        [tidx,tidx]=unique(depths);
        top(tidx,:)=bot(tidx,:);
        model=[depths top];
    else
        [tidx,tidx]=unique(depths,'first');
        bot(tidx,:)=top(tidx,:);
        model=[depths bot];
    end
else
    % get index range (assumes depths are always non-decreasing in model)
    idx1=find(model(:,1)>range(1),1);
    idx2=find(model(:,1)<range(2),1,'last');
    
    % are range points amongst the knots?
    tf=ismember(range,model(:,1));
    
    % if they are, just use the knot, otherwise interpolate
    if(tf(1))
        idx1=idx1-1;
    else
        vtop=interp1q(...
            model(idx1-1:idx1,1),model(idx1-1:idx1,2:end),range(1));
    end
    if(tf(2))
        idx2=idx2+1;
    else
        vbot=interp1q(...
            model(idx2:idx2+1,1),model(idx2:idx2+1,2:end),range(2));
    end
    
    % clip model
    model=model(idx1:idx2,:);
    
    % pad range knots if not there
    if(~tf(1)); model=[range(1) vtop; model]; end
    if(~tf(2)); model=[model; range(2) vbot]; end
end

% fix interpolation of infinity
model(isnan(model))=inf;

% array to struct
mout.name='PREM';
mout.ocean=ocean & crust;
mout.crust=crust;
mout.isotropic=true;
mout.refperiod=1;
mout.flattened=false;
mout.depth=model(:,1);
mout.vp=model(:,2);
mout.vs=model(:,3);
mout.rho=model(:,4);
mout.qk=model(:,5);
mout.qu=model(:,6);

end
