function tt=taupTime(model,depth,phase,varargin)

% TAUPTIME calculate travel time using TauP toolkit
%
% tt=taupTime(model,depth,phase,'option',value)
%
% Input arguments:
% The first three arguments are fixed:
%   Model:      Global velocity model. Default is "iasp91".
%   Depth:      Event depth in km
%   Phase:      Phase list separated by comma
% The other arguments are variable:
%   'deg' value:   Epicentral distance in degree
%   'km'  value:   Epicentral distance in kilometer
%   'sta'/'station' value:  Station location [Lat Lon]
%   'evt'/'event' value:    Event location [Lat Lon]
% 
% Output argument:
%   tt is a structure array with fields:
%   tt(index).time
%            .distance
%            .srcDepth
%            .rayParam
%            .phase
%
% Example:
%   taupTime([],50,'P,S','deg',45.6)
%   taupTime('prem',50,'Pdiff,PKP','sta',[-40 -100],'evt',[30,50])
%
% This program is only a wrapping program for TauP toolkit, 
% which is developed by:
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

import edu.sc.seis.TauP.*;
import java.io.*;
import java.lang.*;
import java.util.*;
import java.util.zip.*;

if nargin<5
    error('At least 5 input arguments required');
end;

if isempty(model)
    model='iasp91';
end;

inArgs{1}='-mod';
inArgs{2}=model;
inArgs{3}='-h';
inArgs{4}=num2str(depth);
inArgs{5}='-ph';
inArgs{6}=phase;
n_inArgs=6;

dist=0;
sta=0;
evt=0;
ii=1;
while (ii<=length(varargin))
    switch lower(varargin{ii})
    case {'deg','km'}
        if ~(isa(varargin{ii+1},'double') & length(varargin{ii+1})==1)
            error('  Incompatible value for option %s !',varargin{ii});
        end;
        inArgs{n_inArgs+1}=['-' varargin{ii}];
        inArgs{n_inArgs+2}=num2str(varargin{ii+1});
        n_inArgs=n_inArgs+2;
        ii=ii+2;
        dist=1;
    case {'sta','station'}
        if ~(isa(varargin{ii+1},'double') & length(varargin{ii+1})==2)
            error('  Incompatible value for option %s !',varargin{ii});
        end;
        inArgs{n_inArgs+1}='-sta';
        temp=varargin{ii+1};
        inArgs{n_inArgs+2}=num2str(temp(1));
        inArgs{n_inArgs+3}=num2str(temp(2));
        n_inArgs=n_inArgs+3;
        ii=ii+2;
        sta=1;
    case {'evt','event'}
        if ~(isa(varargin{ii+1},'double') & length(varargin{ii+1})==2)
            error('  Incompatible value for option %s !',varargin{ii});
        end;
        inArgs{n_inArgs+1}='-evt';
        temp=varargin{ii+1};
        inArgs{n_inArgs+2}=num2str(temp(1));
        inArgs{n_inArgs+3}=num2str(temp(2));
        n_inArgs=n_inArgs+3;
        ii=ii+2;
        evt=1;
    otherwise
        error('  Unknown option %s \n',varargin{ii});
    end; %switch
end; %for
if ~(dist | (sta & evt))
    error('  Event/source locations or distance not specified !');
end;

%disp(inArgs);

try
    arrivals=MatTauP_Time.run_time(inArgs);
catch
    error('Java exception occurred! Please check model name.');
end;

if arrivals.length==0
    fprintf('   Phases do not exist at specified distance!\n');
end;

if nargout==0
    for ii=1:arrivals.length
        fprintf('  Phase: %-10s  Time: %.3f(s) \n', ...
            char(arrivals(ii).getName),arrivals(ii).getTime);
    end;
    return;
end;

tt = [];
for ii=1:arrivals.length
    tt(ii).time=arrivals(ii).getTime;
    tt(ii).distance=arrivals(ii).getDistDeg;
    tt(ii).srcDepth=arrivals(ii).getSourceDepth;
    tt(ii).phaseName=char(arrivals(ii).getName);
    tt(ii).rayParam=arrivals(ii).getRayParam;
end;
