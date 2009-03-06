function c = mkdelwhitespace(c)
% MKDELWHITESPACE            delete whitespace from string
%
% this function removes any whitespace character from a given
% string. whitespace are ' ' and characters with ascii <32
%
% the following charatcer code table is assumed:
%
% ' ' = 32
% Tab = 9
%
% the result of MKDELWHITESPACE will be zero for the characters in that small
% table only.
%
% Martin Knapmeyer, 04.04.2002

t =  c>32; %~ (((c>=48)&(c<=57)) | (c==46) | (c==43) | (c==45));
c=c(t);

