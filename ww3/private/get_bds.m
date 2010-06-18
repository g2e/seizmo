%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% GET THE BDS STRUCTURE %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bds_struct=get_bds(fid,lenbds)

bds=fread(fid,11);
bds_struct.len=lenbds;
bds_struct.oct4=bds(4);
bds_struct.bsfE=int2(bds(5),bds(6));
%bds_struct.RefVal=ibm2fltmex5(bds(7:10));
bds_struct.RefVal=ibm2flt(bds(7:10));
bds_struct.nbits=bds(11);

bds=fread(fid,lenbds-11);
bds_struct.bindata=uint8(bds);
