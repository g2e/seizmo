function grib_struct=read_grib(gribname,irec,varargin)
%READ_GRIB    Reads certain WMO GRiB-formatted binary data files
%
%  Inputs: The first 2 arguments are required:
%          gribname - filename containing GRiB records.
%          irec - specifies which GRiB records to read.
%                 If irec is a vector, it specifies which GRiB records to return.
%                 If irec is a scalar, is specifies how far to read into the GRiB file.
%                 If irec==-1, READ_GRIB reads all records(default). 
%                 Irec can be a CELL ARRAY of parameter names to extract.  
%                 Type read_grib('paramtable') for a list of parameter names. 
%                 Irec can also be the string 'inv{entory}', so that READ_GRIB prints 
%                 a GRiB contents list.
%         
%          READ_GRIB accepts the following property/value pairs:
%
%          HeaderFlag - (0|1) report only headers if==1 (default=1) no data 
%                       structures are returned unless DataFlag==1.
%          DataFlag   - (0|1) return decoded BDS if==1 (default=1).  The data for 
%                       the parameter is stored in the .fltarray field of the structure.
%          ScreenDiag - (0|1) control diagnostics to screen (default=1)
%          ParamTable - ('NCEPOPER'|'NCEPREAN'|'ECMWF128'|'ECMWF160') selects the parameter 
%                       table to use for matching kpds6 number to the correct parameter 
%                       name. (default='NCEPOPER')
%
%          Additional help comments can be obtained by passing one argument to READ_GRIB.  
%          Try: 'mainhelp','paramtableshort','paramtablelong','output_struct'.
%
% Outputs: grib_struct - a MATLAB structure containing the GRiB contents for the requested records.
%                        Call read_grib('output_struct') for a better description.
%
% Call as:grib_struct=read_grib(gribname,irec,p1,v1,p2,v2,...)
%

% Written by Brian Blanton
%            Renaissance Computing Institute
%            Uni. of North Carolina at Chapel Hill
%            Chapel Hill, NC, USA
%            brian_blanton@renci.org
%
%            The decoding of the GRiB BDS is done with mex files that use code
%            from wgrib by Wesley Ebisuzaki at NCEP (http://wesley.wwb.noaa.gov).
%            BDSunpk -->> BDS_unpack_mex5.c
%
% HISTORY: Version 1.2 (Sep 2002):
%             BOB:    pre-2002: first writing
%             BOB:    Jun 2002: added "find_grib_marker" routine to catch info between  
%                               "end_of_grib" marker (7777) and "beginning_of_grib" marker (GRIB).
%             BOB: 06 Sep 2002: fixed bug in HEX2DEC call in INT2 and INT3 functions, used 
%                               in GET_PDS function.  
%                               See comments in INT2. 
%                               Also fixed "Help" code.
%          Version 1.3 (Jun 2003)
%             BOB: 01 Jun 2003: added varargins, generalized parameter table code, 
%                               added ECMWF 128,160, and NCEPOPER tables.
%          Version 1.3.1 (11 Nov 2003)
%             BOB: 11 Nov 2003: fixed bug in HeaderFlag varargin
%          Version 1.4.0 (1 Sep 2005)
%             BOB: 01 Sep 2005: Added gds tables for more Data Representation Types (DRT)
%          Version 1.4.3 (6 Nov 2008)
%             BOB: 06 Nov 2008: Added gds for Polar Stereographic grids (DRT)
%          Version 1.4.4 (20 Nov 2008)
%             BOB: 20 Nov 2008: Fixed handling of garbage in between grib records
%                               Fixed handling of year/century        
%
%          Mid-Jun 2009 : moved code to SVN repo.  1.4.4 corresponds to revision 1.  
%          Revision 2: added comments
%          Revision 3: Replaced undef value of 10e20 with NaN, within BDS_unpack_mex5.c
% 

global rgversion rgversiondate
rgversion='r3';
rgversiondate=datenum(2010,2,5);

if nargin==0,gribhelp('mainhelp');return,end
if nargin==1,gribhelp(gribname);return,end

%if nargin~=5,error('READ_GRIB must have 0|4 input arguments'),end

if isempty(gribname)
   [fname,fpath]=uigetfile('*','Which GRiB File?');
   if fname==0,return,end
   gribname=[fpath fname];
else
   [fpath,fname,fext,fver]=fileparts(gribname);
   fname=[fname fext];
end

inventory_flag=0;
select_by_param=0;
if isempty(irec)
   irec=-1;
elseif iscell(irec)
   select_by_param=1;
   params_to_get=irec;
   irec=-1;
elseif isstr(irec)
   if strncmp(lower(irec),'inv',3)
      inventory_flag=1;
      ScreenDiag=0;
      dataskip=1;
      irec=-1;
   else
      error('String value for irec MUST be "inventory"')
   end
elseif irec==-1
   % This means get entire GRiB file
else
  [m,n]=size(irec);
  if m~=1&n~=1
     error('Size of irec must be 1xn or nx1')
  end
   irect=irec(:);
   [sirect,irecsort]=sort(irect);
   irect=sirect;
end

% Default propertyname values
HeaderFlag=1;
DataFlag=1;
ScreenDiag=1;
global ParamTable
ParamTable='NCEPOPER';

% Process propertyname/value pairs 
k=1;
while k<length(varargin),
  switch lower(varargin{k}),
    case 'headerflag',
      HeaderFlag=varargin{k+1};
      if ~(HeaderFlag==1 | HeaderFlag==0)
         error('Invalid HeaderFlag to READ_GRIB.')
      end
      varargin([k k+1])=[];
    case 'dataflag',
      DataFlag=varargin{k+1};
      if ~(DataFlag==1 | DataFlag==0)
         error('Invalid DataFlag to READ_GRIB.')
      end
      varargin([k k+1])=[];
    case 'screendiag',
      ScreenDiag=varargin{k+1};
      if ~(ScreenDiag==1 | ScreenDiag==0)
         error('Invalid ScreenDiag to READ_GRIB.')
      end
      varargin([k k+1])=[];
    case 'paramtable',
      ParamTable=varargin{k+1};
      if ~(any(strcmp(ParamTable,{'NCEPOPER','ECMWF160','NCEPREAN','ECMWF128'})))
         error('Invalid Parameter Table to READ_GRIB.')
      end
      varargin([k k+1])=[];
    otherwise
      k=k+2;
  end
end

% Initialize parameter table 
global Parameter_Table  
Set_Parameter_Table(ParamTable,ScreenDiag);


% Open GRiB file and Find first 'GRIB' string
fid=fopen(gribname,'r');
if fid<0
   error(['GRiB file named ' gribname ' not found.'])
end

if ~is_grib_file(fid,ScreenDiag)
   % No first GRiB marker found. Abort.
   disp('ERROR: READ_GRIB:')
   disp([gribname ' is not a GRiB file.'])
   disp(['No GRiB marker "GRIB" found in file.'])
   grib_struct=-1;
   return
else
   mark1=ftell(fid);
   fseek(fid,0,-1);   % rewind the GRiB file
   % The -4 is because of the first GRIB marker found previously.
   if mark1>4
      gribfileheader=fread(fid,mark1-4);
   else
      gribfileheader=[];      
   end
   % Report pre-firstmarker header if length <= 90
   if(mark1<91)
      if ScreenDiag, disp(sprintf('GRiB header = %s',char(gribfileheader'))),end
   else
      % Reposition fid to point just before "G" char
      fseek(fid,-4,0);
   end
end

first4octets=fread(fid,4);
if ScreenDiag,disp(char(first4octets'));,end

if inventory_flag
   inventory_str=sprintf('###################################################\n');
   inventory_str=[inventory_str sprintf('Inventory for GRiB file %s \n',fname)];
   inventory_str=[inventory_str sprintf('###################################################\n')];
   inventory_str=[inventory_str sprintf('\n')];
   inventory_str=[inventory_str sprintf(...
   '           Parameter\n')];
   inventory_str=[inventory_str sprintf(...
   '  Rec #  #    Name  Units              Date    Hr Mn   P1    P2  Quantity                Level           IC        Grid     Description\n')];
   inventory_str=[inventory_str sprintf(...
   '--------------------------------------------------------------------------------------------------------------------------\n')];
end

% initialization
grec=0;crec=0;
fpos=ftell(fid);

% Each record in a GRiB file contains sections 0-5, possibly with 
% sections 2,3 (GDS,BMS) omitted.  Basically, each GRIB record 
% contains a parameter (weather model variable like U_FLX)
% and definitions of the model grid, etc.  So the size of the 
% grib_struct will increase with each new record. 

grib_struct=[];

while ~isempty(first4octets)

   oct1to4=char(first4octets'); %'

   % Make sure first 4 octets are "GRIB"
   if ~strcmp(oct1to4,'GRIB')
      % we've already read in 4 bytes at the end of the loop, looking for
      % "GRIB", but we could have read in "  GR", for example. So, rewind 4
      % bytes, and start looking for "GRIB" from there. 
      fseek(fid,-4,0); 
      % scanning for next grib marker
      iret=find_grib_marker(fid,ScreenDiag);
      if iret<1 && feof(fid)
         if ScreenDiag
	    str=sprintf('\n%s\n',['End of file reached.']);
            disp(str)
	 end
         return
      end
   end      

   grec=grec+1;
   this_one_extracted=0;
   
   if ScreenDiag,fprintf(1,'GRec=%4d : FPos1=%8d',grec,fpos),end
   
   % read section 0, the IS (Indicator Section);
   % this is always extracted
   oct5to7=fread(fid,3);
   lengrib=bitshift3(oct5to7(1),oct5to7(2),oct5to7(3));
   edition=fread(fid,1);

   % Edition check
   if edition~=1
      str=sprintf('\n%s\n',['READ_GRIB cannot read edition ' int2str(edition) ' GRiB records.']);
      str=[str sprintf('%s\n',['Edition ' int2str(edition) ' detected at record ' int2str(grec) ' in ' fname]);];
      disp(str)
      grib_struct=[];
      return
   end
   
   % Determine if this is a record to skip or keep
   headerskip=1;dataskip=1;
   if select_by_param
      % Look forward into the PDS of the current GRiB 
      % record to see what parameter this is
      current_fpos=ftell(fid);
      oct1to3=fread(fid,3);
      lenpds=bitshift3(oct1to3(1),oct1to3(2),oct1to3(3));
      fseek(fid,-3,0);
      pds=fread(fid,lenpds);
      param=Get_Parameter(pds(9),1);
      if any(strcmp(param,params_to_get))
          %"any(strcmp(param,params_to_get))" equivalent to "ismember(param,params_to_get)"
          headerskip=0;
          if DataFlag,dataskip=0;end
          crec=crec+1;     
      end
      fseek(fid,current_fpos,-1);
  elseif irec(1)==-1
      headerskip=0;
      if DataFlag,dataskip=0;end
      crec=crec+1;     
  else
      idx=find(irect==grec);
      if ~isempty(idx)
          headerskip=0;
          if DataFlag,dataskip=0;end
          crec=crec+1;
          irect(find(irect==sirect(crec)))=[];
      end
   end
    
  if headerskip & dataskip
      % Skip GRiB record, except for '7777' delimiter at end;
      % We''ve already read in 8 octets!!
      fseek(fid,lengrib-12,0);
   else
      [pds_struct,gds_struct,bms_struct,bds_struct,dataarray]=...
           extract_grib(fid,grec,fpos,headerskip,dataskip);
      % Package current grib records into returning grib_struct
      grib_struct(crec).sec1_1=oct1to4;
      grib_struct(crec).lengrib=lengrib;
      grib_struct(crec).edition=edition;
      grib_struct(crec).file=gribname;
      grib_struct(crec).record=grec;
      grib_struct(crec).description=pds_struct.description;
      temp1=pds_struct.parameter;
      [temp2,temp3]=table3a(pds_struct.pdsvals(10));
      grib_struct(crec).parameter=pds_struct.parameter;
      grib_struct(crec).layer=temp3;
      grib_struct(crec).units=pds_struct.units;
      dtime=datenum(pds_struct.year,pds_struct.month,pds_struct.day,pds_struct.hour,pds_struct.min,00);
      grib_struct(crec).stime=datestr(dtime);
      % 14 Jun 2002: BOB: added level to output struct
      grib_struct(crec).level=levels(pds_struct.pdsvals(10),pds_struct.pdsvals(11),pds_struct.pdsvals(12));
      grib_struct(crec).gridtype=gds_struct.DRT;
      grib_struct(crec).pds=pds_struct;
      grib_struct(crec).gds=gds_struct;
      grib_struct(crec).bms=bms_struct;
      grib_struct(crec).bds=bds_struct;
      grib_struct(crec).fltarray=dataarray;
      this_one_extracted=1;
      if inventory_flag
         clear grib_struct;
         crec=0;  % Reset structure logging
	 minute=['0' int2str(pds_struct.min)];
	 if pds_struct.min>9,minute=int2str(pds_struct.minute);,end
	 hour=['0' int2str(pds_struct.hour)];
	 if pds_struct.hour>9,hour=int2str(pds_struct.hour);,end
	 day=['0' int2str(pds_struct.day)];
	 if pds_struct.day>9,day=int2str(pds_struct.day);,end
	 month=['0' int2str(pds_struct.month)];
	 if pds_struct.month>9,month=int2str(pds_struct.month);,end
         level=levels(pds_struct.pdsvals(10),pds_struct.pdsvals(11),pds_struct.pdsvals(12));
	 l=length(level);lm=floor(l/2);ll=1:lm;lr=lm+1:l;
	 level_left=level(ll);level_right=level(lr);
	 inventory_str=[inventory_str sprintf(...
	            '%4d    %3d %6s  %-14s  %2d%2d/%2s/%2s %2s %2s %5.1f %5.1f %-16s %9s%-9s  %10s %10s %-s\n',...
	            grec,...
		    pds_struct.pdsvals(9),...
		    pds_struct.parameter,...
		    pds_struct.units,...
	            pds_struct.century-1,...
		    pds_struct.year,...
		    month,day,hour,minute,...
		    pds_struct.P1,pds_struct.P2,pds_struct.TRI,level_left,level_right,...
		    pds_struct.IC,gds_struct.DRT(1:8),...
		    pds_struct.description)];
      end
   end   
   
   % The last 4 octets (bytes) should contain the ascii chars 55;
   % this is the string '7777' end-of-grib record delimiter.
   end_grib_delim=char(fread(fid,4)');

   % Get file pointer position 
   fpos=ftell(fid);
   if ScreenDiag
      fprintf(1,' : GLen=%8d : FPos2=%8d : EOGRiB=%4s',lengrib,fpos,end_grib_delim)
      if ~this_one_extracted
         fprintf(1,'\n')
      else
         fprintf(1,' Extracted\n')
      end

   end
   if ~strcmp(end_grib_delim,'7777')
      disp(['Alignment problem reading GRiB message ' gribname])
      disp(['at GRiB record number ' int2str(grec)])
      disp(['Should be at the end of GRiB record, and we''re not.'])
      disp(['Returning from READ_GRIB with GRiB up to this point.'])
      break
%BOB      error('read_grid alignment problem a')   
   end
   
   % Read the next 4 octets, which SHOULD be the beginning of the next GRIB record
   first4octets=fread(fid,4);
   temp=char(first4octets');
   if ~strcmp(temp,'GRIB')
      if ScreenDiag
         disp(sprintf('   Next 4 octets is NOT end of grib.  Searching for next beginning ...',temp))
      end
   end
   
end

if exist('irecsort'),grib_struct=grib_struct(irecsort);,end
if inventory_flag
      helpwin(inventory_str,'Inventory',['Inventory for ' gribname]);
      grib_struct=[];
end

% Permute if needed
if exist('irecsort')
   grib_struct=grib_struct(irecsort);
end


% close grib stream
fclose(fid);
if(ScreenDiag); disp(' '); end
return




