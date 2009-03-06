function t = mkisletter(c)
%MKISLETTER               is true for non-ciphers
% 
%  this function returns 1 for other characters than 0123456789-+.
%  and is valid only on SUN4-computers. t will be of the same size as c.
%  This function is VERY similar to the original ISLETTER, but has a
%  smaller set of characters that yield zero. The original form is
%  useless.
% 
% 
%  the following charatcer code table is assumed:
% 
%  '0' = 48
%  '1' = 49
%  '2' = 50
%  '3' = 51
%  '4' = 52
%  '5' = 53
%  '6' = 54
%  '7' = 55
%  '8' = 56
%  '9' = 57
%  '.' = 46
%  '+' = 43
%  '-' = 45
% 
%  the result of MKISLETTER will be zero for the characters in that small
%  table only.
% 
%  Martin Knapmeyer, 01.07.1995


t =  ~ (((c>=48)&(c<=57)) | (c==46) | (c==43) | (c==45));


