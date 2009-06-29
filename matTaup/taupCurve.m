function tt_curve=taupCurve(model,depth,phase)

% TAUPCURVE calculate travel time curve using TauP toolkit
%
% taupTime(model,depth,phase)
%
% Input arguments:
%   Model:      Global velocity model. Default is "iasp91".
%   Depth:      Event depth in km
%   Phase:      Phase list separated by comma
% 
% Output argumet:
%   tt is a structure array with fields:
%   tt(index).phaseName
%            .srcDepth
%            .distance (in degree)
%            .time
%   If no output argument specified, travel timve curves will be plotted.
%
% Example:
%   taupCurve([],50,'P,sS')
%   taupCurve('prem',0,'P,PKP,PKIKP,PKiKP')
%
% This program calls TauP toolkit for calculation, which is 
% developed by:
%   H. Philip Crotwell, Thomas J. Owens, Jeroen Ritsema
%   Department of Geological Sciences
%   University of South Carolina
%   http://www.seis.sc.edu
%   crotwell@seis.sc.edu
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   Nov, 2002
%
% Modified:
%  June 29, 2009 (gge) added defaults for depth,phase; minor code cleaning
%                      changed sourceDepth to srcDepth


import edu.sc.seis.TauP.*;
import java.io.*;
import java.lang.*;
import java.util.*;
import java.util.zip.*;

if nargin<1 || isempty(model)
    model='iasp91';
end;
if nargin<2 || isempty(depth)
    depth=0;
end;
if nargin<3 || isempty(phase)
    phase='ttall';
end;

inArgs{1}='-mod';
inArgs{2}=model;
inArgs{3}='-h';
inArgs{4}=num2str(depth);
inArgs{5}='-ph';
inArgs{6}=phase;

try
    matCurve=MatTauP_Curve.run_curve(inArgs);
catch
    fprintf('Java exception occurred! Please check input arguments. \n\n');
    return;
end;

tt(1:matCurve.length)=struct('time',[],'distance',[],'srcDepth',[],...
    'phaseName',[],'rayParam',[]);
for ii=1:matCurve.length
    tt(ii).phaseName=char(matCurve(ii).phaseName);
    tt(ii).srcDepth=matCurve(ii).sourceDepth;
    tt(ii).time=matCurve(ii).time;
    tt(ii).distance=matCurve(ii).dist;
    tt(ii).rayParam=matCurve(ii).rayParam;
end;

c={'b','r','g','m','c','y', ...
   'b--','r--','g--','m--','c--','y--', ... 
   'b-.','r-.','g-.','m-.','c-.','y-.', ...
   'b:','r:','g:','m:','c:','y:', ...
   'b','r','g','m','c','y', ...
   'b--','r--','g--','m--','c--','y--', ... 
   'b-.','r-.','g-.','m-.','c-.','y-.', ...
   'b:','r:','g:','m:','c:','y:', ...
   'b','r','g','m','c','y', ...
   'b--','r--','g--','m--','c--','y--', ... 
   'b-.','r-.','g-.','m-.','c-.','y-.', ...
   'b:','r:','g:','m:','c:','y:'};
if nargout==0
    clf;hold on;box on
    n=0; p=cell(1,1);
    for ii=1:numel(tt)
        if numel(tt(ii).distance)>1
            n=n+1;
            plot(tt(ii).distance,tt(ii).time,c{ii});
            p{n}=tt(ii).phaseName;
        end;
    end;
    legend(p,2,'fontsize',6);
    xlabel('Distance (deg)');
    ylabel('Travel Time (s)');
    title(upper(model));
    return;
end;

tt_curve=tt;
