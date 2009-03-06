function datapath=mkttboxdata;
% mkttboxdata....determine path to TTBOX data directory
%
% call: datapath=mkttboxdata;
%
% result: datapath: a string containing the path to the TTBOX data
%                   directory, in which several velocity model files
%                   and reference data for TTBOX validation are stored.
%
%                   DATAPATH ends with a filesep character (a slash on
%                   Unix/Linux, a backslash on Windows systems) so that
%                   [s filename] is a valid path.
%
% This routine determines the path under which it is residing and assumes
% that the data is stored in a directory 'data' below that directory.
%
% By using this, you can avoid to type long paths every time you need
% one of teh TTBOX standard files.
%
% Martin Knapmeyer, 11.02.2004


%%% which is the name of this routine?
myname='mkttboxdata.m';

%%% where does this routine reside?
s=which(myname);

%%% strip off my own name
s=s(1:(length(s)-length(myname)));

%%% add subdirectory DATA to path
s=[s 'data' filesep];

%%%return result
datapath=s;