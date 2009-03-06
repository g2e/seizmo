function m_gshhs(resolution,varargin);
% M_GSHHS Add a coastline to a given map using 
%           the Global Self-consistant Hierarchical High-resolution 
%           Shorelines.
%
%         M_GSHHS(RES, (standard line option,...,...) ) draws the coastline
%         as a simple line.
%         M_GSHHS(RES,'patch' ( ,standard patch options,...,...) ) draws the 
%         coastline as a number of patches. 
%
%         M_GSHHS(RES,'save',FILENAME) saves the extracted coastline data
%         for the current projection in a file FILENAME. This allows 
%         speedier replotting using M_USERCOAST(FILENAME). 
%    
%         RES - selections resolution
%                  1  or 'crude'	
%                  2  or 'low'  	
%                  3  or 'intermediate'  
%                  4  or 'high' 	
%                  5  or 'full  	
%
%
%
%         See also M_PROJ, M_GRID, M_COAST, M_GSHHS_L, M_GSHHS_H, M_GSHHS_C 
%         M_USERCOAST    

% Rich Pawlowicz (rich@ocgy.ubc.ca) 15/June/98
%
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.
%
%  16/Dec/2005
%*********************************************************************
%  Modified after code provided by Bruce Lipphardt (brucel@udel.edu) to 
%  reduce the hierarchy of M_GSHHS_* routines to a single routine with a
%  variable resolution input:
%

% Root of directories where gshhs_X.b files live
FILNAME='private/';


res_list = char('c','l','i','h','f') ;

if isstr(resolution),
 resolution = strmatch(lower(resolution(1)),res_list);
end;

if isempty(resolution) | resolution<1 | resolution> length(res_list),
  error('**Don''t recognize the specified resolution');
end;
  
res_char = res_list(resolution) ;
file     = [FILNAME,sprintf('gshhs_%s.b',res_char)] ;
tag_name = sprintf('m_gshhs_%s',res_char) ;


% Set current projection to geographic
Currentmap=m_coord('set');
m_coord('geographic');


if length(varargin)>1 & strcmp(varargin{1},'save'),
  [ncst,Area,k]=mu_coast(res_char,file);
  eval(['save ' varargin{2} ' ncst k Area']);
else
  mu_coast(res_char,file,varargin{:},'tag',tag_name);
end;

m_coord(Currentmap.name);

