function clr=mkreadclr(pfad,silent);
% mkreadclr........read continuous layer representation file
%
% call: clr=mkreadclr(pfad);
%       clr=mkreadclr(pfad,silent);
%
%           pfad: full pth to the .clr file to be read
%           silent: if present, nothing is written to the screen.
%                   otherwise, the programm produces some output on read progress.
%
% result: clr: a CLR structure representing the file content.
%              This structure has the following fields:
%
%                .name: model name
%                .year: publication year
%              .planet: name of planet for which model is valid
%                  .rp: radius of planet [km]
%                       MUST be defined in file!
%              .lyrcnt: number of leyers defined by model
%              .layers: a structure array which contains the layer specifications.
%                       .layers(i) is the sub-structure for the i'th layer, which
%                       is not necessarily the i-th layer from the surface!!
%                       This sub-structure has the folowing fields:
%                       .depth: depth extent of the layer, measured in km from the
%                               surface. Depths defined as radii in the file are
%                               transformed into depth.
%                               min(depth) is upper boundary, max(depth) is lower boundary.
%                        .name: name of layer
%                          .vp: P wave velocity polynomial coefficients
%                               the i-th coefficient .vp(i) is associated with the
%                               (i-1)-th power of the normalized radius in the polynomial.
%                          .vs: S wave velocity polynomial coefficients
%                         .rho: density polynomial coefficients
%                          .qp: P wave Q factor polynomial coefficients
%                          .qs: S wave Q factor polynomial coefficients
%                       If any quantity is not defined in the file, NaN is returned.
%                .conr: depth of conrad discontinuity [km]
%                       All discontinuity depths are NaN if not defined.
%                .moho: depth of mohorovicic discontinuity
%                .d410: depth of olivine alpha - beta phase transition discontinuity
%                       (these fields are named after standard depths on earth - just
%                        to have short field names)
%                .d520: depth of olivine beta - gamma phase transition discontinuity
%                .d660: depth of olivine gamma - perovskite boundary (lower mantle)
%                 .cmb: depth of core mantle boundary
%                 .icb: depth of outer core - inner core boundary
%                  .dz: list of depths of other named discontinuities
%                       will be empty if no non-standard doscontinuities are defined.
%               .dname: list of names of other named discontinuities
%                       .dname(i) is the name of discontinuity at depth .dz(i)
%                       will be empty if no non-standard doscontinuities are defined.
%
%
% Example for evaluation:
% to comput vp in depth z, take the .vp coefficient of the correct layer
% and compute
%
% x=(rp-z)/rp; % normalized radius
% velocity=vp(1)+vp(2)*x+vp(3)*x^2 ... ; 
%
% Martin Knapmeyer, 01.09.2003


%%% init result
clr=mkmakeclr('clr');

%%% init helper structure for radius/depth conversion
coord.conr=[];
coord.moho=[];
coord.d410=[];
coord.d520=[];
coord.d660=[];
coord.cmb=[];
coord.icb=[];
coord.dz=[];
coord.layers=[];


%%% run silent?
if nargin==2
   silent=1;
else
   silent=0;
end; % if nargin==2


%%% open file
[fid,msg]=fopen(pfad,'r');
if fid==-1
   error(['MKREADCLR: ' msg]);
end; % if fid

%%% read & interpret lines
done=0;
linecnt=1; % line counter
while done==0
   if (feof(fid)) %|(linecnt>25)
      %%% end of file reached. quit file loop.
      done=1;
   else
      %%% we gonna get another line from the file!
      %%% read the next line
      line=fgetl(fid);
      
            
      %%% strip off comments
      line=mkstripcomments(line);
      %disp(['stripped:' line]);
      
      %%% split line into tokens
      [tokcnt,tokenlist]=mksplitline(line);
      
      %%% identify keyword and modifier
      [keyword,modifier,parmlist]=mksplittokenlist(tokenlist,linecnt,silent);
      
      %%% consolidation of parmlist
      [numparms,stringparms]=mkparmconsolid(keyword,modifier,parmlist,silent);
      
      %%% build identified parameters into field of CLR structure
      [clr,coord]=mkparm2clr(clr,coord,keyword,modifier,numparms,stringparms);
      
      %%% increment line counter
      linecnt=linecnt+1;
   end; % if feof(fid)
end; % while done

%%% close file
fclose(fid);


%%% transform depth representations into radius
clr=mkradius2depth(clr,coord);



%%% leave MKREADCLR
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 HELPER FUNCTIONS                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newline=mkstripcomments(oldline)
% mkstripcomments.......strip off .clr file comments from line
%
% call: newline=mkstripcomments(oldline)
%
%       oldline: string containing a line from a .clr file
%
% result: newline: the same als OLDLINE but all comments stripped off.
%
% comments begin with '#', '//', or '/*' and end at the end of the line.
%
% Martin Knapmeyer, 01.09.2003

%%% find all comment start sequences
indies=[findstr(oldline,'#') findstr(oldline,'//') findstr(oldline,'/*')];

%%% remove all comments
if ~isempty(indies)
   %%% a comment exists
   indies=sort(indies);
   commentstart=indies(1);
   newline=oldline(1:(commentstart-1));
else
   %%% no comment exists
   newline=oldline;
end; % if ~isempty

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [tokcnt,tokenlist]=mksplitline(line);
% mksplitline.......split text line into tokens used in .clr files
%
% call: [tokcnt,tokenlist]=mksplitline(line);
%
%       line: string containing a line from a .clr file
%
% result: tokcnt: number of tokens found
%         tokenlist: strin matrix, containing on token from LINE per row
%
% A token is defined as a sequence of non-whitespace characters, where
% whitespace are ASCII characters 9 and 32.
%
% Martin Knapmeyer, 01.09.2003

%%% init result
tokenlist=[];

%%% split line
done=0;
while done==0
    [tok,rem]=strtok(line);
    if ~isempty(tok)
       %%% add token to tokenlist
       tokenlist=strvcat(tokenlist,tok);
       line=rem;
    else
       %%% no more tokens
       done=1;
    end; % if ~isempty
end; % while done

tokcnt=size(tokenlist,1);

%disp(['MKSPLITLINE: ' int2str(tokcnt) ' tokens found.']);


return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [keyword,modifier,parmlist]=mksplittokenlist(tokenlist,linecnt,silent);
% mksplittokenlist.........divide token list into keyword, modifier, and parameter list
%
% call: [keyword,modifier,parmlist]=mksplittokenlist(tokenlist,linecnt,silent);
%
%       tokenlist: string matrix containing token list, as returned by MKSPLITLINE
%       linecnt: line counter, needed in error messages
%       silent: verbosity flag:
%               1: run silent
%               0: run verbose
%
% result: keyword: string containing the keyword (without "!")
%                   empty if no keyword present
%         modifier: string containing the modifier to the keyword (without "!"),
%                   empty if no modifier present
%         parmlist: string matrix, containing all parameters
%                   empty if no parameters present
%
% Keywords and Modifiers are recognized by 
% a) position in TOKENLIST and
% b) preceeding exclamation mark
%
% Martin Knapmeyer, 01.09.2003


%%% init output
keyword=[];
modifier=[];
parmlist=[];

tokcnt=size(tokenlist,1);
if tokcnt==0
   %%% token list is empty
else
   %%% token list is not empty
   keyword=deblank(tokenlist(1,:));
   if strcmp(keyword(1),'!')~=1
      error(['MKSPLITTOKENLIST: Exclamation mark missing at begin of line '...
             int2str(linecnt) ' of file.']);
   else
      keyword=keyword(2:end); % strip off exclamation mark
      if tokcnt>1
         %%% other tokens might be modifier and parameters
         modifier=deblank(tokenlist(2,:));
         if strcmp(modifier(1),'!')~=1
            %%% not identified as modifier ("!" missing)
            modifier=[];
            parmlist=tokenlist(2:tokcnt,:);
         else
            modifier=modifier(2:end); % strip off exclamation mark
            if tokcnt>2
               %%% remaining tokens are pramaters
               parmlist=tokenlist(3:tokcnt,:);
            end; % if tokcnt>2
         end; % strcmp(modifier(1),'!')~=1 else
      else
         %%% no more tokens
         modifier=[];
         parmlist=[];
      end; % if tokcnt>1 else
   end; % if strcmp((keyword(1),'!')~=1
end; % if tokcnt

%%% tell them what we've found
if ~silent
   if ~isempty(keyword)
      disp(['MKSPLITTOKENLIST: keyword ' upper(keyword) ' identified.']);
   end; % if ~isempty(keyword)
   if ~isempty(modifier)
      disp(['MKSPLITTOKENLIST: modifier ' upper(modifier) ' identified.']);
   end; % if ~isempty(modifier)
end; % if ~silent

%%% case insensitive output of keywords and modifiers
keyword=lower(keyword);
modifier=lower(modifier);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [numparms,stringparms]=mkparmconsolid(keyword,modifier,parmlist,silent);
% mkparmconsolid..........analysis & consolidation of parameter list
%
% call: [numparms,stringparms]=mkparmconsolid(keyword,modifier,parmlist,silent);
%
%       keyword: the identified keyword (without then "!")
%       modifier: the identified modifier (without the "!")
%       parmlist: the raw parameter list
%       silent: verbosity flag:
%               1: run silent
%               0: run verbose
%
%       Input parameters are as returned by MKSPLITTOKENLIST
%
% result: numparms: numeric column vector containing all numeric parameters
%         stringprams: string row vector containing the string parameter (only one possible)
%
% The input matrix PARMLIST is a very raw version of a parameter list. It is split
% up into the parameter values as intended by the file.
% To do so, names of keywords and modifiers are hard-wired within this routine.
%
% What this routine does is essentially the following: look at which keyword and
% modifier are given and thencut and paste the lines of PARMLIST into the suitable
% formats.
% Additionally, it checks many syntactical and semantical issues.
%
% Martin Knapmeyer, 02.09.2003


%%% init result
numparms=[];
stringparms=[];


if ~isempty(keyword)
   %%% there is a keyword, so there might be a modifier and parameters
   
   %%% make things not case sensitive
   keyword=lower(keyword);
   modifier=lower(modifier);

   switch keyword
      case {'name'}
           if ~isempty(modifier)
               error(['MKPARMCONSOLID: no modifier allowed for keyword '...
                       upper(keyword) '.']);
           else
               stringparms=mkstrmat2line(parmlist);
           end; % if ~isepmty(modifier)
      case {'year'}
           if ~isempty(modifier)
               error(['MKPARMCONSOLID: no modifier allowed for keyword '...
                       upper(keyword) '.']);
           else
               buffy=mkstrmat2line(parmlist);
               numparms=str2num(buffy);
           end; % if ~isepmty(modifier)

      case {'planet'}
           if ~isempty(modifier)
              switch modifier
                 case {'name'}
                      stringparms=mkstrmat2line(parmlist);
                 case {'radius'}
                      buffy=mkstrmat2line(parmlist);
                      numparms=str2num(buffy);
                 otherwise
                     error(['MKPARMCONSOLID: modifier ' upper(modifier)...
                             'illegal for keyword' upper(keyword)]);
              end; % switch modifier
           else
              error(['MKPARMCONSOLID: modifier required for keyword '...
                     upper(keyword) '.']);
           end; % ~isempty(modifier)
      case {'layer'}
           if ~isempty(modifier)
              switch modifier
                 case {'start'}
                      stringparms=mkstrmat2line(parmlist);
                 case {'depth','radius','vp','vs','rho','qp','qs'}
                      buffy=mkstrmat2line(parmlist);
                      numparms=str2num(buffy);
                 case {'end'}
                      %%% no parameters allowed
                      if ~isempty(parmlist)
                         error(['MKPARMCONSOLID: no parameter allowed for modifier '...
                                upper(modifier) ' of keyword ' upper(keyword) '.']);
                      end; % if ~isempty(parmlist)
                 otherwise
                     error(['MKPARMCONSOLID: modifier ' upper(modifier)...
                             ' illegal for keyword ' upper(keyword)]);
              end; % switch modifier
           else
              error(['MKPARMCONSOLID: modifier required for keyword '...
                     upper(keyword) '.']);
           end; % ~isempty(modifier)

      case {'discon'}
            if ~isempty(modifier)
              switch modifier
                 case {'depth','radius'}
                      tokcnt=size(parmlist,1);
                      if tokcnt<2
                         error(['MKPARMCONSOLID: keyword ' upper(keyword),...
                                ' requires two parameters!']);
                      else
                         numparms=str2num(parmlist(1,:));
                         stringparms=mkstrmat2line(parmlist(2:tokcnt,:));
                      end; % if tokcnt<2 else
                 otherwise
                     error(['MKPARMCONSOLID: modifier ' upper(modifier)...
                             'illegal for keyword' upper(keyword)]);
              end; % switch modifier
           else
              error(['MKPARMCONSOLID: modifier required for keyword '...
                     upper(keyword) '.']);
           end; % ~isempty(modifier)
      case {'usertag'}
           if ~isempty(modifier)
               error(['MKPARMCONSOLID: no modifier allowed for keyword '...
                       upper(keyword) '.']);
           else
               stringparms=mkstrmat2line(parmlist);
           end; % if ~isepmty(modifier)
      otherwise
          error(['MKPARMCONSOLID: unknown keyword ' upper(keyword)]);
   end; % switch keyword
   
   %%% format output
   numparms=numparms(:);
   stringparms=stringparms(:)';

   %%% tell them what we've found
   if ~silent
      if ~isempty(numparms)
         disp(['MKPARMCONSOLID: numeric parameters [' num2str(numparms') ']']);
      end; % if ~isempty
      if ~isempty(stringparms)
         disp(['MKPARMCONSOLID: string parameters "' stringparms '"']);
      end; % if ~isempty
   end; % if ~silent
   
else
   %%% there is no keyword   
end; % if ~isempty(keyword)




return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function line=mkstrmat2line(strmat);
% mkstrmat2line.......concatenates all lines of a string matrix into a single line
%
% call: line=mkstrmat2line(strmat);
%
%       strmat: string matrix as generated by strvcat()
%
% result: line: a single-line string in which all rows of STRMAT are concatenated
%
% Rows of STRMAT are deblanked and then separated by a single ASCII 32 character.
%
% Martin Knapmeyer 03.09.2003

%%% init result
line=[];

%%% how many lines?
linecnt=size(strmat,1);

%%% loop over all lines
for indy=1:linecnt
    if indy==1
       line=[deblank(strmat(indy,:))];
    else
       line=[line char(32) deblank(strmat(indy,:))];
    end; % if indy==1
end; % for indy

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [newclr,newcoord]=mkparm2clr(oldclr,oldcoord,keyword,modifier,numparms,stringparms);
% mkparm2clr.........builds CLR structure field from identified keywords & parameters
%
% call: [newclr,newcoord]=mkparm2clr(oldclr,oldcoord,keyword,modifier,numparms,stringparms);
%
%       oldclr: existing CLR structure
%       oldcoord: COORD structure which stores the coordinate type
%                 ('depth' or 'radius') for each layer and discontinuity depth
%                 fields are named as in CLR-structures, but contain the modifier words
%                 given with each declaration.
%       keyword: keyword identified by MKSPLITTOKENLIST
%       modifer: modifier identified by MKSPLITTOKENLIST
%       numparms: numeric paramater list returned by MKPARMCONSOLID
%       stringparms: string parameter returned by MKPARMCONSOLID
%
% result: newclr: new CLR structure with additional fields as defined by input.
%         newcoord: updated version of OLDCOORD.
%
% For definition of the CLR structure see help lines of MKREADCLR.
%
%
% Martin Knapmeyer 03.09.2003

% MKPARMCONSOLID has already checked syntax. Here we can assume aeverything is OK!
% we assume that
% - KEYWORD and MODIFIER are lowercase
% - NUMPARMS and STRINGPARMS are not empty when numeric or string parms are required


%%% make numparms a row vector
numparms=numparms(:)';

%%% init output
newclr=oldclr;
newcoord=oldcoord;

%%% init empty dummy layer
newlayer=mkmakeclr('layer');


if ~isempty(keyword)

   %%% dispatch input to respective fields
   switch keyword
       case {'name'}
            newclr.name=stringparms;
       case {'year'}
            newclr.year=numparms;
       case {'planet'}
            switch modifier
                case {'name'}
                     newclr.planet=stringparms;
                case {'radius'}
                     newclr.rp=numparms;
                otherwise
                     error(['MKPARM2CLR: unknown modifier ' upper(modifier)]);
            end; % switch modifier
       case {'layer'}
             switch modifier
                case {'start'}
                     newclr.layers=[newclr.layers newlayer];
                     newclr.lyrcnt=newclr.lyrcnt+1;
                     newclr.layers(newclr.lyrcnt).name=stringparms;
                case {'depth','radius'}
                     newclr.layers(newclr.lyrcnt).depth=numparms;
                     newcoord.layers=strvcat(newcoord.layers,modifier);
                case {'vp'}
                     newclr.layers(newclr.lyrcnt).vp=numparms;
                case {'vs'}
                     newclr.layers(newclr.lyrcnt).vs=numparms;
                case {'rho'}
                     newclr.layers(newclr.lyrcnt).rho=numparms;
                case {'qp'}
                     newclr.layers(newclr.lyrcnt).qp=numparms;
                case {'qs'}
                     newclr.layers(newclr.lyrcnt).qs=numparms;
                case {'end'}
                     % nothing special happens
                otherwise
                     error(['MKPARM2CLR: unknown modifier ' upper(modifier)]);
             end; % switch modifier
       case {'discon'}
            %%% store standard discontinuities within respective fields
            %%% and non-standard discontinuities in a list
            stringparms=lower(stringparms);
            switch stringparms
                case {'conrad'}
                     newclr.conr=numparms;
                     newcoord.conr=modifier;
                case {'moho','mantle'}
                     newclr.moho=numparms;
                     newcoord.moho=modifier;
                case {'olivine alpha beta','transition zone'}
                     newclr.d410=numparms;
                     newcoord.d410=modifier;
                case {'olivine beta gamma'}
                     newclr.d520=numparms;
                     newcoord.d520=modifier;
                case {'olivine gamma perovskite','lower mantle'}
                     newclr.d660=numparms;
                     newcoord.d660=modifier;
                case {'cmb','outer core','outer-core'}
                     newclr.cmb=numparms;
                     newcoord.cmb=modifier;
                case {'icb','inner core','inner-core'}
                     newclr.icb=numparms;
                     newcoord.icb=modifier;
                otherwise
                     %%% non-standard-discontinuity
                     newclr.dz=[newclr.dz numparms];
                     newclr.dname=strvcat(newclr.dname,stringparms);
                     newcoord.dz=strvcat(newcoord.dz,modifier);
            end; % switch stringparms
       case {'usertag'}
            newclr.tag=stringparms;
       otherwise
           error(['MKPARM2CLR: unknown keyword ' upper(keyword)]);
   end; % switch keyword

end; % if ~isempty

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newclr=mkradius2depth(oldclr,coord);
% mkradius2depth.........transform radius coordinate into depth coordinate
%
% call: newclr=mkradius2depth(oldclr,coord);
%
%       oldclr: complete CLR structure
%       coord: COORD structure as defined in MKREADCLR.
%              For each definition of a radius/deopth coordinate in OLDCLR, this
%              structure contains a string flag that determines if it is a radius
%              or a depth.
%
% result: newclr: the same as OLDCLR, but all radius coordinates are transformed
%                 into depth coordinates.
%
% Algorithm: if COORD indicates a radius, the value is replaced byb oldclr.rp minus
%            this radius.
%
% Martin Knapmeyer, 02.09.2003

%%% prepare result
newclr=oldclr;



%%% is a radius defined?
if isempty(newclr.rp)
   error('MKRADIUS2DEPTH: No planetary radius defined!');
end; % if isempty(newclr.rp)

%%% transform layers
for indy=1:newclr.lyrcnt
    if strcmp(coord.layers(indy,:),'radius')==1
       newclr.layers(indy).depth=newclr.rp-newclr.layers(indy).depth;
    end; % if strcmp(coord.layers(indy),'radius')==1
end; % for indy

%%% transform standard discontinuities
if strcmp(coord.conr,'radius')==1
   newclr.conr=newclr.rp-newclr.conr;
end; % if strcmp
if strcmp(coord.moho,'radius')==1
   newclr.moho=newclr.rp-newclr.moho;
end; % if strcmp
if strcmp(coord.d410,'radius')==1
   newclr.d410=newclr.rp-newclr.d410;
end; % if strcmp
if strcmp(coord.d520,'radius')==1
   newclr.d520=newclr.rp-newclr.d520;
end; % if strcmp
if strcmp(coord.d660,'radius')==1
   newclr.d660=newclr.rp-newclr.d660;
end; % if strcmp
if strcmp(coord.cmb,'radius')==1
   newclr.cmb=newclr.rp-newclr.cmb;
end; % if strcmp
if strcmp(coord.icb,'radius')==1
   newclr.icb=newclr.rp-newclr.icb;
end; % if strcmp

%%% transform non-standard discontinuities
dzcnt=length(newclr.dz);
for indy=1:dzcnt
    if strcmp(coord.dz(indy,:),'radius')==1
       newclr.dz(indy)=newclr.rp-newclr.dz(indy);
    end; % if strcmp(coord.layers(indy),'radius')==1
end; % for indy