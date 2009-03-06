function [model,dz,dname]=mkreadnd(pfad,silent);
% mkreadnd............read "named discontinuity"-file as used by TauP Toolkit
%
% call: [vmod,dz,dname]=mkreadnd(pfad);
%       model=mkreadnd(pfad);
%         ...=mkreadnd(pfad,silent);
%
%       pfad: path to the .nd/file
%      silent: if present, MKREADND does not write any comment to the screen
%              the value of silent is ignored.
%
% results: vmod: n-by-6 array containing the model parameters
%                  vmod[:,1]: depth [km below surface]
%                  vmod[:,2]: Vp [km/s]
%                  vmod[:,3]: Vs [km/s]
%                  vmod[:,4]: rho [g/ccm]
%                  vmod[:,5]: Qp (attenuation factor for P waves)
%                  vmod[:,6]: Qs (attenuation factor for S waves)
%          dz: array of discontinuity depths
%          dname: array of discontinuity names dname(i) corresponds to dz(i)
%          model: if only one output parameter is specified, this will be a structure with the
%                 following fields:
%                 model.z: depth [km below surface]
%                 model.vp: Vp
%                 model.vs: Vs
%                 model.rho: rho [g/ccm]
%                 model.qp: Qp
%                 model.qs: Qs
%                 model.conr: depth of conrad discontinuity
%                 model.moho: depth of moho
%                 model.d410: depth of Mantle-Transition Zone-discontinuity (the "410" on earth, where
%                             the olivine-alpha-beta-transition occurs)
%                 model.d520: depth of the olivine-beta-gamma transition (520km on earth)
%                 model.d660: depth of lower mantle discontinuity (the "660" on earth, where the 
%                             olivine-gamma-perovskite transition occurs)
%                 model.cmb: depth of core mantle boundary
%                 model.icb: depth of inner core boundary
%                 model.dz: depths of additional discontinuities
%                 model.dname: names of additional discvontinuities
%                 model.rp: planetary radius
%                           Defined usgin the !radius-Keyword. If this is nbot given,
%                           the largest value of z is assumed to be the planetary radius!
%                 model.name: name of model, determined by !name keyword
%                 model.year: year of publication of model, defined through the !year-Keyword.
%                             If this is not given, the year will be empty
%
%                 if a discontinuity is not specified in the .nd-file, its depth
%                 will be NaN in MODEL.
%
%                 in order to identify the standard discontinuities within the .nd-file,
%                 use the following discontinuity names:
%                 "conrad" -> model.conr
%                 "moho" or "mantle" -> model.moho
%                 "transition zone" or "olivine alpha beta" -> model.d410
%                 "olivine beta gamma" -> model.d520
%                 "lower mantle" or "olivine gamma perovskite" -> model.d660
%                 "outer core" or "outer-core" -> core mantle boundary
%                 "inner core" or "inner-core" -> inner core boundary
%                 names are not case sensitive.
%
%                 The routine recognizes several keywords to define special fields of the output
%                 These Keyewords are preceeded by an exclamation mark "!". Between the Keyword and its value,
%                 a Space (' ') is expected.
%                 Keywords are:
%
%                 !radius: defines the planets radius
%                 !year: year of publication of model
%                 !name: model name
%
%                 if a keyword occurs more than once, the last value is returned.
%                 Keywords may occur at any place of the file, but have to be at the beginning of a line.
%                 Lines with keywords are allowed to contain comments.
%
% In case of error, vmod==-1 and discont==[].
%
% The TauP Manual (version 1.1) states that the following discontinuity names are recognized by
% TauP Tooklit: "mantle" for Moho, "outer core" for core-mantle-boundary, "inner core" for outer
% core-inner-core boundary.
% However, this routine accepts any string as discontinuity name.
%
% a simple sample model file is
%
% /* below is a simple named discontinuities model. */
% !year 2002
% !radius 6371
% !name SimpleSample
% 0.0 5.0 3.0 2.7
% 20 5.0 3.0 2.7
% 20 6.5 3.7 2.9
% 33 6.5 3.7 2.9
% mantle # the word "mantle" designates that this is the moho
% 33 7.8 4.4 3.3
% 410 8.9 4.7 3.5 // unnamed discontinuity!
% 410 9.1 4.9 3.7
% 670 10.2 5.5 4.0
% 670 10.7 5.9 4.4
% 2891 13.7 7.2 5.6
% outer-core # "outer-core" designates that this is the core mantle boundary
% 2891 8.0 0.0 9.9
% 5149.5 10.3 0.0 12.2
% inner-core # "inner-core" makes this the inner-outer core boundary
% 5149.5 11 3.5 12.7
% 5160 // no values are defined at this depth.
% 5200 -1 -1 -1  # another way to express undefined values. will be returned as NaN
% 6371 11.3 3.7 13
% // end of simple sample
%
% values are separated by whitespace.
% As you see,the files may contain comments begining with "/*, "//" or "#".
% Not all possible parameters are necessarily given.
% Values are assumed to be NaN if not specified. NaN is physically nonsense for all parameters and may therefore
% serve as identifier for invalid or undefined values. (The routine returns what is given in the file
% and does not make assumptions on what is not given)
% if any parameter in the file is -1, NaN is returned.
% discontinuity names have to contain at least one letter! (otherwise undistinguishable
% from an incomplete parameter list!)
%
% Martin Knapmeyer, 04.04.2002, 10.04.2002, 19.06.2002, 05.07.2002,
% 25.04.2006
%
% BUGS: cannot handle comments longer than one line (those beginning with "/*' in one line
%       and ending with '*/' in another). reason: line oriented reading of files.


% struct-to-double assingment warning fixed, MK25042006


%%%%% prepare output vars
vmod=-1;     % returns this in case of error!
dz=[];       % returns this in case of error!
dname=[];    % returns this in case of error!
model=-1;    % returns this in case of error!

%%%% declarations of variables used for keywords
pradius=[];
pubyear=[];

%%%%% check input parameters
if nargin==0
   error('MKREADND: Please specify path to .nd file!');
end; % if nargin==0
if exist(pfad)==0
   error(['MKREADND: File ' pfad ' does not exist!']);
end; % if exist==0

%%%%% open and read file
[fid,msg]=fopen(pfad,'r');
if fid==-1
   %%%%% error during file opening
   warning(msg);
   disp(['MKREADND: problems reading velocity model file'])
else
   %%%%% vmod default
   vmodini=[-1 -1 -1 -1 -1 -1]*NaN;
   vmod=vmodini;
   %%%%% read and interpret
   parmcnt=1;
   namecnt=1;
   linecnt=1;
   done=0;
   while done==0
      %% read current line
      line=fgetl(fid);
      if isempty(line)
         line=' ';
      end; % if isempty
      if line==-1
         done=1;
      else
         %% strip /**/ comments
         indy=findstr(line,'/*');
         if ~isempty(indy)
            line=line(1:indy);
         end; % if indy
         %% strip // and # comments
         [token,remainder1]=strtok([' ' line],'/');
         if ~isempty(remainder1)
            %disp(['MKREADND: stripping comment:' remainder1]);
         end; % if ~isempty
         token=deblank(token);
         [token,remainder2]=strtok([' ' token],'#'); % token
         if ~isempty(remainder2)
            %disp(['MKREADND: stripping comment:' remainder2]);
         end; % if isemtpy(remainder)
         token=deblank(token);
         %% evaluate the token
         %% at this place, the remaining string may have one of the following forms:
         %% - empty
         %% - contain up to six numbers separated by whitespace
         %% - contain a discontinuity name
         %% - contain a keyword.
         %%disp(['MKREADND: remaining parameter text: ' token]);
         if ~isempty(token)
            %% if token contains non-number, non-whitespace characters, it is a name
            %% otherwise, it is a parameter list.
            if sum(mkisletter(mkdelwhitespace(token)))~=0
               %% name or keyword
               keyword=mkdelwhitespace(token);
               if strcmp(keyword(1),'!')
                  % keyword
                  [keyword,parm]=strtok(token,' ');
                  switch lower(keyword)
                     case {'!radius'}
                        pradius=str2num(parm);
                     case {'!year'}
                        pubyear=str2num(parm);
                     case {'!name'}
                        mname=parm;
                     otherwise
                        disp(['MKREADND: Keyword ' keyword ' not recognized. Ignored.']);
                  end; % switch
               else
                  % name
                  dz(namecnt)=vmod(parmcnt-1,1);
                  dname=strvcat(dname,token);
                  if nargin==1
                     disp(['MKREADND: discontinuity: ' token ' at ' num2str(dz(namecnt)) ' km']);
                  end; % if nargin
                  namecnt=namecnt+1;
               end; % if strcmp
            else
               %% parm list
               pcnt=1;
               done2=0;
               vmod(parmcnt,:)=vmodini; % initalize new row with defaults (inf for Qp, Qs!)
               while done2==0
                  [token,rem]=strtok(token);
                  if isempty(rem)
                     done2=1; % all existing parameters found
                  end; % if isempty(rem) end
                  vmod(parmcnt,pcnt)=str2num(token);
                  pcnt=pcnt+1;
                  token=rem;
               end; % while done2
               parmcnt=parmcnt+1;
            end; %% if sum(mkisletter()) else
         end; % if ~isempty(token)
         linecnt=linecnt+1;
      end; % if line==-1 else
   end; % while done==0
   
   
   %%% close file
   fclose(fid);
end; % if fid==-1


%%%%% replace all -1 by NaN
indies=find(vmod==-1);
vmod(indies)=vmod(indies)*NaN;
indies=find(dz==-1);
dz(indies)=dz(indies)*0-1;


if isempty(pradius)
   pradius=max(vmod(:,1));
end; % if isempty

%%%%% return results
model=vmod;

%%% except when a model-structure is to be returned
if nargout==1
   model=struct; % to avoid the struct-to-double assingment warning
   %%% th trivial ones...
   model.z=vmod(:,1);
   model.vp=vmod(:,2);
   model.vs=vmod(:,3);
   model.rho=vmod(:,4);
   model.qp=vmod(:,5);
   model.qs=vmod(:,6);
   model.name=mname;
   model.rp=pradius;
   model.year=pubyear;
   %%% the less trivial ones...
   model.conr=NaN;
   model.moho=NaN;
   model.d410=NaN;
   model.d520=NaN;
   model.d660=NaN;
   model.cmb=NaN;
   model.icb=NaN;
   model.dz=[];
   model.dname=[];
   anz=length(dz); % so many discontinuities in dz/dname
   cnt=1; % counter for nonstandard discontinuities
   for indy=1:anz
      currentname=mkdeblank(lower(dname(indy,:)));
      switch currentname
         case {'conrad'}
            model.conr=dz(indy);
         case {'moho','mantle'}
            model.moho=dz(indy);
         case {'olivine alpha beta','transition zone'}
            model.d410=dz(indy);
         case {'olivine beta gamma'}
            model.d520=dz(indy);
         case {'olivine gamma perovskite','lower mantle'}
            model.d660=dz(indy);
         case {'cmb','outer core','outer-core'}
            model.cmb=dz(indy);
         case {'icb','inner core','inner-core'}
            model.icb=dz(indy);
         otherwise
            if nargin==1
               disp(['MKREADND: non-standard discontinuity: ' currentname]);
            end; % if nargin
            model.dz=[model.dz dz(indy)];
            model.dname=strvcat(model.dname,currentname);
         end; % switch
   end; % for indy
      

end; % if nargout


if nargin==1
      disp(['MKREADND: planet radius is ' num2str(pradius)]);
      disp(['MKREADND: ' num2str(namecnt-1) ' discontinuity names.']);
      disp(['MKREADND: ' num2str(parmcnt-1) ' parameter sets.']);
      disp(['MKREADND: ' num2str(linecnt-1) ' lines total.']);
   end; % if nargin
