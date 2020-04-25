function [ map ] = DcMapping( size, width )
%ConduitDiameter Mapping provides mapping vector for the expectation that the
%longest vessels are also the widest, in general
%   Longest vessels should be the
%   widest and widest vessels should possess the wider pores as a result of
% a stochastic process.
%   Here, a vector, that assigns every vessel a number taken from a sorted
%   length or diameter vector, is returned
%   size is the number of vessels or inter-vessel largest pores
%   width is the size of each rank: rank of size 5 means that the 5 longest
%   vessels obtain, randomly, the 5 largest pit pores

ranks = floor(size/width);
rest = mod(size,ranks);
shift = cumsum([0 ones(1,ranks-1).*width]);

% The trick is to add an extra unit to the width when the rest is non-zero.
% The size of rest is definitely less than ranks, so first, create a zero
% array of the size of ranks.
permSizeInc = zeros(1,ranks);
%Fill out 1 in as many locations as the number of 'rest' starting from the
%end.
permSizeInc(end-rest+1:end) = 1;
% Add 1 to each number in shift, starting from the end, until all rests
% have been accomodated.
shift(end-rest+2:end) = shift(end-rest+2:end) + ...
    cumsum(permSizeInc(end-rest+2:end));
map = zeros(1,size);
for k = 1:ranks
    map(shift(k)+(1:width+permSizeInc(k))) = ...
        randperm(width+permSizeInc(k))+shift(k);
end


end

