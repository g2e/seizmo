%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% GET THE BMS STRUCTURE %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bms_struct=get_bms(fid,lenbms)
bms=fread(fid,lenbms);
bms_struct.len=lenbms;               % Length in octets of Bit Map Section
bms_struct.NUnused_Bits=bms(4);      % Number of unused bits at end of Section 3.
bms_struct.oct5=bms(5);
bms_struct.oct6=bms(6);
bms_struct.bitmap=uint8(bms(7:lenbms)); % Bit map, zero filled to an even number of octets
bms_struct.bmsvals=bms;
