%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TABLE A, PDS Octet 6                             %
% Orig Center #7 (NCEP) Generating Process/Model   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gen_proc=tableA(ival)
switch ival
   case 02, gen_proc='Ultra Violet Index Model';
   case 05, gen_proc='Satellite Derived Precipitation and temperatures, from IR';
   case 10, gen_proc='Global Wind-Wave Forecast Model';
   case 19, gen_proc='Limited-area Fine Mesh (LFM) analysis';
   case 25, gen_proc='Snow Cover Analysis';
   case 30, gen_proc='Forecaster generated field';
   case 31, gen_proc='Value added post processed field';
   case 39, gen_proc='Nested Grid forecast Model (NGM)';
   case 42, gen_proc='Global Optimum Interpolation Analysis (GOI) from AVN run';
   case 43, gen_proc='Global Optimum Interpolation Analysis (GOI) from "Final" run';
   case 44, gen_proc='Sea Surface Temperature Analysis';
   case 45, gen_proc='Coastal Ocean Circulation Model';
   case 49, gen_proc='Ozone Analysis from TIROS Observations';
   case 52, gen_proc='Ozone Analysis from Nimbus 7 Observations';
   case 53, gen_proc='LFM-Fourth Order Forecast Model';
   case 64, gen_proc='Regional Optimum Interpolation Analysis (ROI)';
   case 68, gen_proc='80 wave tri., 18-layer Spec. model from AVN run';
   case 69, gen_proc='80 wave tri., 18 layer Spec. model from MRF run';
   case 70, gen_proc='Quasi-Lagrangian Hurricane Model (QLM)';
   case 73, gen_proc='Fog Forecast model - Ocean Prod. Center';
   case 74, gen_proc='Gulf of Mexico Wind/Wave';
   case 75, gen_proc='Gulf of Alaska Wind/Wave';
   case 76, gen_proc='Bias corrected MRF';
   case 77, gen_proc='126 wave tri., 28 layer spec. model from AVN run';
   case 78, gen_proc='126 wave tri., 28 layer spec. model from MRF run';
   case 79, gen_proc='Backup from the previous run';
   case 80, gen_proc='62 wave triangular, 28 layer Spectral model from MRF run';
   case 81, gen_proc='Spectral Statistical Interpolation (SSI) anal from AVN run';
   case 82, gen_proc='Spectral Statistical Interpolation (SSI) anal from "Final" run';
   case 83, gen_proc='MESO ETA Model - Backup Version (currently 80 km)';
   case 84, gen_proc='MESO ETA Model (currently 32 km)';
   case 85, gen_proc='No longer used';
   case 86, gen_proc='RUC Model, from FSL (isentropic; scale: 60km at 40N)';
   case 87, gen_proc='CAC Ensemble Forecasts from Spectral (ENSMB)';
   case 88, gen_proc='Ocean Wave model with additional physics (PWAV)';
   case 90, gen_proc='62 wave tri., 28 layer spec. model extension of MRF run';
   case 91, gen_proc='62 wave tri., 28 layer spec. model extension of AVN run';
   case 92, gen_proc='62 wave tri., 28 layer spec. model run from MRF final analysis';
   case 93, gen_proc='62 wave tri., 28 layer spec. model run from T62 GDAS anal of MRF run';
   case 94, gen_proc='T170/L42 Global Spectral Model from MRF Run';
   case 95, gen_proc='T126/L42 Global Spectral Model from MRF Run';
   case 96, gen_proc='Aviation Model (currently T170/L42 Global Spectral Model)';
   case 100, gen_proc='RUC Surface Analysis (scale: 60km at 40N)';
   case 101, gen_proc='RUC Surface Analysis (scale: 40km at 40N)';
   case 105, gen_proc='RUC Model from FSL (isentropic; scale: 40km at 40N)';
   case 110, gen_proc='ETA Model - 15km version';
   case 150, gen_proc='NWS River Forecast System (NWSRFS)';
   case 151, gen_proc=' NWS Flash Flood Guidance System (NWSFFGS)';
   case 152, gen_proc='WSR-88D Stage II Precipi Anal';
   case 153, gen_proc='WSR-88D Stage III Precip Anal';
   otherwise, gen_proc='N/A';
end
