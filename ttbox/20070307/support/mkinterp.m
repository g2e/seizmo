function yneu=mkinterp(xalt,yalt,xneu,mode);
% mkinterp..............one dimensional interpolation of nonmonotonic functions
%
% call: yneu=mkinterp(xalt,yalt,xneu,mode);
%
%            xalt: x values at which function f(x) is given
%            yalt: function values yalt=f(xalt)
%            xneu: x values at which function f(x) is requested
%                  all values of XALT have to be between the smallest XALT and
%                  the largest XALT (or equal to these) - the routine does not
%                  do EXTRApolation, but INTERpolation.
%            mode: interpolation mode
%                  This determines how to handle multiple y for identical x
%                  the following options are supperoted:
%
%                  'data': the mean of all y1(x1), y2(x1), y3(x1),...
%                          is returned. This makes sense if the y are
%                          measured values of some parameter x and the y
%                          being different because of measurement error
%                  'jump': the first of all y being in the list is
%                          returned. I introduced this option specifically
%                          for use with TTBOX, where this case occurs when
%                          interpolating discontinuous velocity models
%
%                  DEFAULT: 'data'
%
% result: yneu:  yneu=f(xneu)
%                three cases will be distinguished:
%                1) xneu is identical to one of the xalt values
%                   yneu will be the corresponding yalt
%                2) xneu is no identical to one of the xalt values
%                   yneu will be determined by linear interpolation
%                3) xneu is identical to several xalt values
%                   yneu will be determined as mean of f(xalt1),
%                   f(xalt2),..., or as first of all yalt in the list.
%                4) xneu is out of the range defined by xalt.
%                   yneu will be the ylat at the respective border of
%                   the range. A warning will be issued.
%
%
% MatLab's interp1-function requires x to be monotonic for purely technical
% reasons: vectorization as done there is impossible for nonmonotonic x series.
% That is ridiculous.
%
% This routine is not vectorized and therefore slower :-( 
% This routine does not make any assumptions on x and y. f(x) does not even
% have to be unique!
% Only linear interpolation is possible.
%
% Martin Knapmeyer 28.01.2000, 25.08.2004, 26.08.2004


%%% nderstand input
nin=nargin;
switch nin
   case {3}
      mode='data';
   case {4}
      mode=lower(mode);
   otherwise
      error(['MKINTERP: unknown mode ' upper(mode)]);
end; % switch nin


xlen=length(xneu);
[minxalt,minindy]=min(xalt);
[maxxalt,maxindy]=max(xalt);
yneu=zeros(size(xneu));

for i=1:xlen
    if (xneu(i)>=minxalt)&(xneu(i)<=maxxalt)
       indies=find(xalt==xneu(i));
       inlen=length(indies);
       switch inlen
          case 0,
               % i-th xneu is truely between the given xalt values
               % search the limits to interpolate between
               xdiff=xalt-xneu(i);
               [xdiffsort,sortindy]=sort(xdiff);
               xalt=xalt(sortindy); % sort x by distance to xneu
               yalt=yalt(sortindy); % sort y by distance of x to xneu
               indy=find(diff(sign(xdiffsort))==2);
               x0=xalt(indy);   % largest x xmaller than xneu(i)
               x1=xalt(indy+1); % smallest x larger than xneu(i)
               y0=yalt(indy);   % y1=f(x1)
               y1=yalt(indy+1); % y2=f(x2)
               % interpolate
               h=x1-x0;
               delta=y1-y0;
               yneu(i)=y0+delta*(xneu(i)-x0)/h;
          case 1,
               % i-th xneu is identical to one of the xalt values
               yneu(i)=yalt(indies);
          otherwise,
               % i-th xneu is identical to several xalt values
               % mode has to be considered
               switch mode
                  case {'data'}
                      % compute mean of all relevant yalt values
                      yneu(i)=mean(yalt(indies));
                  case {'jump'}
                      % return first of all relevant yalt values
                      yneu(i)=yalt(min(indies));
               end; % switch mode
       end; % switch
    else
       warning(['MKINTERP: XNEU(' int2str(i)...
             ') is outside x range defined by XALT!']);
       if xneu(i)<=minxalt
          yneu=yalt(minindy);
       end; % if xneu(i)<=minxalt
       if xneu(i)>=maxxalt
          yneu=yalt(maxindy);
       end; % if nxeu(i)>=maxxalt
    end; % if (xneu(i)>=minxalt)&(xneu(i)=<maxxalt)
end; % for i
