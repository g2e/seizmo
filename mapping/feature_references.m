function []=feature_references()
%FEATURE_REFERENCES    References for features available in MAPFEATURE
%
%Volcano locations are from the Smithsonian Institution, Global Volcanism Program.
%Data is available here:
%http://www.volcano.si.edu/world/globallists.cfm
%
%Some plate boundary data was obtained from the Plates Project UTIG.
%Data is available here:
%http://www.ig.utexas.edu/research/projects/plates/data.htm
%
%Large Igneous Province data obtained from the Plates Project of UTIG.
%Data is available here (but some is conspicuously missing compared to the figure):
%ftp://ftp.ig.utexas.edu/pub/LIPS/Data/
%
%Isochrons, extinct ridges, cont-ocean boundaries data:
%http://www.earthbyte.org/Research/Current/digit_isochrons.html
%ftp://ftp.es.usyd.edu.au/pub/agegrid/1997/
%
%Updated sea floor age model:
%http://www.earthbyte.org/Research/Current/agegrid2008.html
%ftp://ftp.earthbyte.org/earthbyte/agegrid/2008/Grids/
%
%Peter Bird plate boundaries (much better than the UTIG set):
%http://peterbird.name/publications/2003_PB2002/2003_PB2002.htm
%http://peterbird.name/oldFTP/PB2002/
%
%Fracture zone & magnetic lineation data:
%http://www.aist.go.jp/GSJ/dMG/dMGold/free/plates/Intro.html
%http://www.soest.hawaii.edu/PT/
%
%Impacts (mine is from an earlier excel of the earth impact database):
%http://keith.aa.washington.edu/craterdata/
%http://impacts.rajmon.cz/
%http://bbs.keyhole.com/ubb/ubbthreads.php?ubb=showflat&Number=1211741
%
%Seamounts:
%ftp://ftp.soest.hawaii.edu/pwessel/
%
%
% I have edited the following datasets:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ridges - removed 1-pt ridges
% cvl - adjusted as my digitization was poor
% jos - adjusted as my digitization was poor
% lips - tons of complex & broken polygons (only 2/5ths retained for now)
% pb2002_plates - altered polygon order and added hacks for polar plates
% seamounts - removed duplicates from 1st release as 2nd is botched
% 
% These need help:
%%%%%%%%%%%%%%%%%%%
% cob - europe is incomplete
% lips - 3/5ths of the data needs fixing
%
help feature_refereences
end
