function snew=mkdeblank(s);
% mkdeblank.......remove heading and trailing blanks from a string
%
% call: snew=mkdeblank(s);
%
%       s: string
%
% result: snew: string with heading and trailing blanks removed.
%
% a blank in the sense of this routine is any character with ASCII<=32
%
% Martin Knapmeyer, 08.05.2002


indies=find(abs(s)>32);

if isempty(indies)
   snew=s;
else
   snew=s(min(indies):max(indies));
end;