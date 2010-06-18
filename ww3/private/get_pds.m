%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% GET THE PDS STRUCTURE %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pds_struct=get_pds(fid,lenpds)

pds=fread(fid,lenpds);%read the whole pds into memory

pds_struct.len=lenpds;              % the length of pds
pds_struct.PTV=pds(4);              % parameter table version number 
pds_struct.IC=table0(pds(5));       % identification of center(table0)
pds_struct.GenProcID=pds(6);        % generating process ID mumber(allocated by the originating center),table A
pds_struct.GridID=pds(7);           % grid identification (geographical location and area, defined by the originating center,table B
temp=table1(pds(8));                % flag specifying the presence or abscence of GDS or a BMS; table 1
pds_struct.HASGDS=temp(1);
pds_struct.HASBMS=temp(2);
pds_struct.parameter=Get_Parameter(pds(9),1);   % parameter name
pds_struct.units=Get_Parameter(pds(9),3);       % parameter units
pds_struct.description=Get_Parameter(pds(9),2); % parameter description
[layerstr,layerabbrev]=table3a(pds(10));        %indicator of type of level of layer, table 3 & 3a

if ~iscell(layerstr)
   pds_struct.layer=sprintf('%d - %s [%s]',int2str(pds(10)),layerstr,layerabbrev);
else
   pds_struct.layer=layerstr;
end

%pds_struct.layer=[int2str(pds(10)) ' - ' layerstr ' [' layerabbrev ']'];
pds_struct.h_p_etc=pds(11:12);    % height, pressure,etc. of the level or layer,table 3
pds_struct.year=pds(13)+100*(pds(25)-1);
pds_struct.month=pds(14);
pds_struct.day=pds(15);
pds_struct.hour=pds(16);
pds_struct.min=pds(17);
pds_struct.fcast_time_unit=table4(pds(18));
pds_struct.P1=pds(19);
pds_struct.P2=pds(20);
pds_struct.TRI=table5(pds(21));
pds_struct.TRI_N=bitshift2(pds(22),pds(23));
pds_struct.TRI_Nmissing=pds(24);
pds_struct.century=pds(25);
pds_struct.SubCID=pds(26);
pds_struct.DecSF=int2(pds(27),pds(28));
pds_struct.pdsvals=pds;
