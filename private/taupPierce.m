function tt_pierce=taupPierce(model,depth,phase,varargin)

% TAUPPIERCE calculate points where specified rays pierce discontinuities
%      using TauP toolkit
%
% tt_pierce=taupPierce(model,depth,phase,'option',value,...) 
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
%   'az' value:    Azimuth (event to station) in degree. This is used to  
%                  calculate lat and lon of pierce points if the event lat lon 
%                  and distance are also given.
%   'baz' value:   Back azimuth (station to event) in degree. This is used to  
%                  calculate lat and lon of pierce points if the event lat lon 
%                  and distance are also given.
%   'dev':         only calculate underside and bottom turn points
%   'turn':        only calculate bottom turning points
%   'under':       only calculate underside reflection points
%   'pierce' depth:  adds depth for calculating pierce points
%   'nodiscon':    only prints pierce points for the depths added
%                     with 'pierce' (have to be after 'pierce' option)
%
% Output argument:
%   tt_pierce is a structure array with fields:
%   tt_pierce(index).time
%            .distance
%            .srcDepth
%            .rayParam
%            .phaseName
%            .pierce.p
%            .pierce.time
%            .pierce.distance
%            .pierce.depth
%            .pierce.latitude
%            .pierce.longitude
%
% Example:
%   taupTime([],50,'P,S',45.6,[],'turn')
%   taupTime('prem',50,'P,PKP',[-40 -100],[30,50],'pierce',100)
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
    case {'dev','under','turn','nodiscon'}
        inArgs{n_inArgs+1}=['-' varargin{ii}];;
        n_inArgs=n_inArgs+1;
        ii=ii+1;
    case {'pierce','az','baz'}
        if ~(isa(varargin{ii+1},'double') & length(varargin{ii+1})==1)
            error('  Incompatible value for option %s !',varargin{ii});
        end;
        inArgs{n_inArgs+1}=['-' varargin{ii}];;
        inArgs{n_inArgs+2}=num2str(varargin{ii+1});
        n_inArgs=n_inArgs+2;
        ii=ii+2;
    otherwise
        error('  Unknown option %s \n',varargin{ii});
    end; %switch
end; %for
if ~(dist | (sta & evt))
    error('  Event/source locations or distance not specified !');
end;

%disp(inArgs);return;

try
    matArrivals=MatTauP_Pierce.run_pierce(inArgs);
catch
    fprintf('Java exception occurred! Please check input arguments. \n\n');
    return;
end;

if nargout==0
    for ii=1:matArrivals.length
        fprintf('  Phase: %-10s  Time: %.3f(s) \n', ...
            char(matArrivals(ii).getName),matArrivals(ii).getTime);
    end;
    return;
end;

tt_pierce = [];
for ii=1:matArrivals.length
    tt(ii).time=matArrivals(ii).getTime;
    tt(ii).distance=matArrivals(ii).getDistDeg;
    tt(ii).srcDepth=matArrivals(ii).getSourceDepth;
    tt(ii).phaseName=char(matArrivals(ii).getName);
    tt(ii).rayParam=matArrivals(ii).getRayParam;
    tt(ii).pierce.p=matArrivals(ii).getMatPath.p;
    tt(ii).pierce.time=matArrivals(ii).getMatPath.time;
    tt(ii).pierce.distance=matArrivals(ii).getMatPath.dist;
    tt(ii).pierce.depth=matArrivals(ii).getMatPath.depth;
    tt(ii).pierce.latitude=matArrivals(ii).getMatPath.lat;
    tt(ii).pierce.longitude=matArrivals(ii).getMatPath.lon;
end;

tt_pierce=tt;
