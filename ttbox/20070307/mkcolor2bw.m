function handles=mkcolor2bw(axhandle);
% mkcolor2bw.......convert TTBOX color plot into b/w plot
%
% call: handles=mkcolor2bw(axhandle);
%
%       axhandle: current axis handle as returned by gca
%
% result: handles: handles to all modified line objects
%
% This routine converts color plots into b/w plots that are cheaper to
% publish. It is pecialized for TTBOX plots: the "tag" propertiy of
% objects is evaluated in a certain predefined way.
%
% Martin Knapmeyer, 08.12.2003


%%% init result
reshandles=[];


%%% identify all line objects
handles=findobj(axhandle,'type','line');


%%% loop through handle list and redefine color and style
%%% according to tag property
if ~isempty(handles)
   for i=1:length(handles);
   
       h=handles(i);
       
       %%% get teag property and decompose into keyword and modifier
       tag=get(h,'tag');
       [keyword,remainder]=strtok(tag);
       modifier=strtok(remainder);
       if isempty(keyword)
          keyword='unknown';
          modifier='unknown';
       else
          reshandles=[reshandles h];
          %disp(['MKCOLOR2BW: identified object: ' keyword ' ' modifier]);
       end; % if isempty keyword
       
       %%% set all colors to black
       set(h,'Color','k');
       set(h,'MarkerEdgeColor','k');
       set(h,'MarkerFaceColor','k');
       
       %%% set line style and width according to keyword and modifier
       keyword=lower(keyword);
       switch keyword
          case {'phase'}
              switch modifier
                  case {'P'}
                      set(h,'LineStyle','-');
                      set(h,'LineWidth',2);
                  case {'PcP','PKiKP'}
                      set(h,'LineStyle',':');
                      set(h,'LineWidth',1)
                  case {'PKP','PKIKP'}
                      set(h,'LineStyle','-');
                      set(h,'LineWidth',1);
                  case {'PP','PPP','PKKP','PcPPcP','PKIKPPKIKP','PKPPKP'}
                      set(h,'LineStyle','-');
                      set(h,'LineWidth',0.5);
                  case {'S'}
                      set(h,'LineStyle','--');
                      set(h,'LineWidth',2);
                  case {'ScS','SKiKS'}
                      set(h,'LineStyle','.-');
                      set(h,'LineWidth',1)
                  case {'SKS','SKIKS'}
                      set(h,'LineStyle','--');
                      set(h,'LineWidth',1);
                  case {'SS','SSS','SKKS','ScSScS','SKIKSSKIKS','SKSSKS'}
                      set(h,'LineStyle','--');
                      set(h,'LineWidth',0.5);
                  case {'SP','PS'}
                      set(h,'LineStyle','-.');
                      set(h,'LineWidth',2);
                  case {'ScP','PcS'}
                      set(h,'LineStyle','-.');
                      set(h,'LineWidth',1);
                  case {'ScPScP','PcSPcS'}
                      set(h,'LineStyle','-.');
                      set(h,'LineWidth',0.5);  
                  case {'unknown'}
                      %%% don't do anything
                  otherwise
                      %error(['MKCOLOR2BW: tag modifier "' modifier '" not recognized']);
              end; % switch modifier
          case {'ray'}
              modifier=lower(modifier);
              switch modifier
                  case {'p'}
                      set(h,'LineStyle','-');
                      set(h,'LineWidth',1);
                  case {'s'}
                      set(h,'LineStyle','-');
                      set(h,'LineWidth',0.25);
                  case {'unknown'}
                      %%% don't do anything
                  otherwise
                      error(['MKCOLOR2BW: tag modifier "' modifier '" not recognized']);
              end; % switch modifier
          case {'unknown'}
             %%% don't do anything
          otherwise
             error(['MKCOLOR2BW: tag keyword "' keyword '" not recognized']);
       end; % switch tag
       

   end; % for i
end; % if ~isempty


%%% return result
handles=reshandles;
