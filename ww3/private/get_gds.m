%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET THE GDS STRUCTURE                                          %
% BOB 28 Oct 2005 split details for Gaussian and Equidist grids  %
%                 added {'Rotated latitude/longitude grid'}      %
% BOB 08 Feb 2007 Added code for thinned/reduced Gaussian Grid   %
%                 handling, as in wgrib                          %
% BOB 06 Nov 1008 Added code for Polar Stereographic grids       %
%                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gds_struct=get_gds(fid,lengds)

gds=fread(fid,lengds);
gds_struct.len=lengds;
gds_struct.NV=gds(4);
gds_struct.PV=gds(5);
gds_struct.DRT=table6(gds(6));

% 01 Apr 2004: BOB GDS table D (Sundry Grid Definitions)

switch gds_struct.DRT
   
   case {'Lambert Conf.'}        % DRT=3
      % For Lambert Conformal Grids
      % GDS Octets 7-42
      gds_struct.Nx=bitshift2(gds(7),gds(8));
      gds_struct.Ny=bitshift2(gds(9),gds(10));
      gds_struct.La1=int3(gds(11),gds(12),gds(13))/1000;
      gds_struct.Lo1=int3(gds(14),gds(15),gds(16))/1000;
      gds_struct.rcf=gds(17);
      gds_struct.LOV=int3(gds(18),gds(19),gds(20))/1000;
      gds_struct.Dx=int3(gds(21),gds(22),gds(23))/1000;
      gds_struct.Dy=int3(gds(24),gds(25),gds(26))/1000;
      gds_struct.Latin1=int3(gds(29),gds(30),gds(31))/1000;
      gds_struct.Latin2=int3(gds(32),gds(33),gds(34))/1000;
      gds_struct.Lat_of_SP=int3(gds(35),gds(36),gds(37))/1000;
      gds_struct.Lon_of_SP=int3(gds(38),gds(39),gds(40))/1000;

   case {'Gaussian Lat/Lon'}     % DRT=4
      % For Gaussian Grids and some others.
      % GDS Octets 7-32
      gds_struct.Ni=bitshift2(gds(7),gds(8));
      gds_struct.Nj=bitshift2(gds(9),gds(10));
      % 8 Feb, 2007 Added code for thinned/reduced Gaussian Grid handling, as in wgrib
      if gds_struct.Ni==65535
         fprintf(1,' ### This is a reduced Gaussian Grid ### ' );
         gds_struct.Ni=-1;
         pl = gds(4) * 4 + gds(5)-1;
         %pl
         %gds_struct.Nj

         isum=0;
         for i = 0:gds_struct.Nj-1
            idx1=(pl+(i)*2);
            idx2=(pl+(i)*2+1);
            iadd1=gds(idx1+1)*256;
            iadd2=gds(idx2+1);
            isum = isum + iadd1+iadd2;
            %sprintf('%d %d %d %d %d\n',i,idx1,idx2,iadd1,iadd2)
            %if i==10, return,end
         end
         gds_struct.NxNy = isum;
      else
         gds_struct.NxNy=gds_struct.Ni*gds_struct.Nj;
      end
      
      gds_struct.La1=int3(gds(11),gds(12),gds(13))/1000;
      gds_struct.Lo1=int3(gds(14),gds(15),gds(16))/1000;
      gds_struct.rcf=gds(17);
      gds_struct.La2=int3(gds(18),gds(19),gds(20))/1000;
      gds_struct.Lo2=int3(gds(21),gds(22),gds(23))/1000;
      gds_struct.Di=bitshift2(gds(24),gds(25))/1000;
      gds_struct.N=bitshift2(gds(26),gds(27));
      gds_struct.smf=gds(28);
      gds_struct.oct29to32=gds(29:32);                                          

   case {'Polar Stereogrphic'}    % DRT=5 
      % For Polar Stereographic grids
      % GDS Octets 7-32
      gds_struct.Ni=bitshift2(gds(7),gds(8));
      gds_struct.Nj=bitshift2(gds(9),gds(10));
      gds_struct.La1=int3(gds(11),gds(12),gds(13))/1000;
      gds_struct.Lo1=int3(gds(14),gds(15),gds(16))/1000;
      gds_struct.rcf=gds(17);
      gds_struct.LOV=int3(gds(18),gds(19),gds(20))/1000;
      gds_struct.Dx=int3(gds(21),gds(22),gds(23))/1000;
      gds_struct.Dy=int3(gds(24),gds(25),gds(26))/1000;
      pds_struct.pcf=gds(27);
      gds_struct.smf=gds(28);
      gds_struct.oct29to32=gds(29:32);      

   case {'Equidis. Cyl. Lat/Lon' , 'Rotated latitude/longitude grid'} % DRT=0,10      
      % For Gaussian Grids and some others.
      % GDS Octets 7-32
      gds_struct.Ni=bitshift2(gds(7),gds(8));
      gds_struct.Nj=bitshift2(gds(9),gds(10));
      gds_struct.La1=int3(gds(11),gds(12),gds(13))/1000;
      gds_struct.Lo1=int3(gds(14),gds(15),gds(16))/1000;
      gds_struct.rcf=gds(17);
      gds_struct.La2=int3(gds(18),gds(19),gds(20))/1000;
      gds_struct.Lo2=int3(gds(21),gds(22),gds(23))/1000;
      gds_struct.Di=bitshift2(gds(24),gds(25))/1000;
      gds_struct.Dj=bitshift2(gds(26),gds(27))/1000;
      gds_struct.smf=gds(28);
      gds_struct.oct29to32=gds(29:32);                                          

   case {'Arakawa Semi-Staggered E-Grid'} % DRT=201
      % Arakawa Semi-Staggered E-Grid on Rotated Latitude/Longitude Grid
      % GDS Octets 7-32
      gds_struct.Ni=bitshift2(gds(7),gds(8));              % Total number of actual points in grid
      gds_struct.Nj=bitshift2(gds(9),gds(10));             % Dummy second dimensionl; set to 1
      gds_struct.La1=int3(gds(11),gds(12),gds(13))/1000;   % Lat of first grid point, millidegrees 
      gds_struct.Lo1=int3(gds(14),gds(15),gds(16))/1000;   % Lon of first grid point, millidegrees
      gds_struct.rcf=gds(17);                              % Resolution and component flag
      gds_struct.La2=int3(gds(18),gds(19),gds(20))/1000;   % Number of mass points along southernmost grid row
      gds_struct.Lo2=int3(gds(21),gds(22),gds(23))/1000;   % Number of rows in each column
      gds_struct.Di=bitshift2(gds(24),gds(25))/1000;       % Longitudinal Direction Increment
      gds_struct.Dj=bitshift2(gds(26),gds(27))/1000;       % Latitudinal Direction Increment
      gds_struct.smf=gds(28);                              % Scanning mode flags
      gds_struct.oct29to32=gds(29:32);                     % Reserved (set to zero)         

   case {'Arakawa Filled E-Grid'} % DRT=202
      % Arakawa Filled E-Grid on Rotated Latitude/Longitude Grid
      % GDS Octets 7-32
      gds_struct.Ni=bitshift2(gds(7),gds(8));              % Total number of actual points in grid
      gds_struct.Nj=bitshift2(gds(9),gds(10));             % Dummy second dimensionl; set to 1
      gds_struct.La1=int3(gds(11),gds(12),gds(13))/1000;   % Lat of first grid point, millidegrees 
      gds_struct.Lo1=int3(gds(14),gds(15),gds(16))/1000;   % Lon of first grid point, millidegrees
      gds_struct.rcf=gds(17);                              % Resolution and component flag
      gds_struct.La2=int3(gds(18),gds(19),gds(20))/1000;   % Number of (zonal) points is each row
      gds_struct.Lo2=int3(gds(21),gds(22),gds(23))/1000;   % Number of (meridional) points in each column
      gds_struct.Di=bitshift2(gds(24),gds(25))/1000;       % Longitudinal Direction Increment
      gds_struct.Dj=bitshift2(gds(26),gds(27))/1000;       % Latitudinal Direction Increment
      gds_struct.smf=gds(28);                              % Scanning mode flags
      gds_struct.oct29to32=gds(29:32);                     % Reserved (set to zero)         
       
   otherwise
      error(['DRT ' int2str(gds(6)) ' not yet coded.'])
end

gds_struct.gdsvals=gds;
