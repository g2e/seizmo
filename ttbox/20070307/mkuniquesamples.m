function [newx,newy,newz]=mkuniquesamples(x,y,z);
% mkuniquesamples........identify non-repeating samples in three lists
%
% call: keepthese=mkuniquesamples(x,y,z);
%
%       x: first list of samples
%       y: second list of samples
%       z: thrid list of samples
%
%       all three lists must have the same number of elements.
%
% result: keepthese: list of indices that makes the three input lists
%                    unique
%
%
% x, y, z are lists of three coordinates in some parameter space that
% may contain repeated entries, where all three values are identical to an
% earlier occurence. This routine constructs an index list KEEPTHESE such
% that x(keepthese), y(keepthese), z(keepthese) contains no repetitions.
%
% It is assumed that the three lists are sorted by X.
% (no internat sorting to save time!)
%
% Martin Knapmeyer, 13.12.2006

%%% how many samples?
sampleanz=length(x);

%%% init result
keepthese=1:sampleanz;


%%% concatenate the three lists into a matrix: each line of the matrix
%%% corresponds to one of the input lists
%%% first make vectors line vectors, if necessary
if size(x,1)~=1
    x=x(:)';
    y=y(:)';
    z=z(:)';
    transposed=1;
else
    transposed=0;
end; % if size(x,1)~=1
%%% the concatenate them into a matrix
xyz=[x; y; z];

%%% reduce matrix to unique elements
reduced=unique(xyz','rows')';

%%% decompose reduced matrix into output lists
newx=reduced(1,:);
newy=reduced(2,:);
newz=reduced(3,:);

%%% make result row vector again
if transposed==1
    newx=newx';
    newy=newy';
    newz=newz';
end; % if transposed==1





