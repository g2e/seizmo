function str_pb = progress_bar(percentage)

% Made by Nicolas Le Roux on this beautiful day of July, 20th, 2005
% Distribute this code as much as you want.
% You can send any comment at lerouxni@iro.umontreal.ca

% lots of edits by Garrett Euler
% - changed o to #  and  . to -
% - reduced bar size (each char is 2%)
% - added progress tip change on the odd percents (uses =)
% - added % to percent

% progress_bar(percentage) draws a progress bar at the position percentage

str_perc = [' ' num2str(percentage) '%%'];

% To draw, we only need to consider the closest integer
percentage = floor(percentage);

if percentage < 51
	str_hash = char(ones(1,fix(percentage/2))*35);
    str_equal = char(ones(1,mod(percentage,2))*61);
	str_dash_beg = char(ones(1, max(0, 25 - numel(str_hash) - numel(str_equal)))*45);
	str_dash_end = char([32 ones(1, 25)*45]);
	str_pb = strcat('[', str_hash, str_equal, str_dash_beg, str_perc, str_dash_end, ']');
else
	str_hash_beg = char(ones(1,25)*35);
	str_hash_end = char(ones(1, max(0, fix((percentage-50)/2)))*35);
    str_equal = char(ones(1,mod(percentage,2))*61);
	str_dash = char(ones(1, 25 - numel(str_hash_end) - numel(str_equal))*45);
	str_pb = char(strcat('[', str_hash_beg, str_perc, {' '}, str_hash_end, str_equal, str_dash, ']'));
end
