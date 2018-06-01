function [ map ] = DcMapping( size, width )
%RAREPITMAPPING provides mapping vector for rare-pit hypothesis 
%   To simulate the rare-pit hypothesis, longest vessels should be the
%   widest and widest vessels should possess the wider pores.
%   Here, a vector, that assigns every vessel a number taken from a sorted
%   length or diameter vector, is returned
%   size is the number of vessels or inter-vessel largest pores
%   width is the size of each rank: rank of size 5 means that the 5 longest
%   vessels obtain, randomly, the 5 largest pit pores

ranks=floor(size/width);
rest=mod(size,ranks);
shift=cumsum([0 ones(1,ranks-1).*width]);
% shift(2:end)=shift(2:end)-1;
permSizeInc=zeros(1,ranks);
permSizeInc(end-rest+1:end)=1;
shift(end-rest+2:end)=shift(end-rest+2:end)+cumsum(permSizeInc(end-rest+2:end));
map=zeros(1,size);
for k=1:ranks
    map(shift(k)+(1:width+permSizeInc(k))) = randperm(width+permSizeInc(k))+shift(k);
end


end

