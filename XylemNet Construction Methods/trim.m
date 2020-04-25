function [MatOut] = trim(Mat,num1,num2,reset)
%This function receives a conduit construction matrix Mat that has three
%values in it: num1 for conduit starts, num2 for conduit ends, and a value
%determined by reset that does not change the constuctor behavior at that
%node. 
%Mat provides each node with a toggle (num1,num2) or passive (reset) number
%that determines whether the constructor should start building a conduit at
%that node, carry on building an already existing conduit, or terminate a
%conduit.
%The goal of this function is to ensure that each column starts with a
%conduit initiation and ends with a conduit terminator
MatOut = zeros(size(Mat));
for j = 1 : size(Mat,3)
    
    [num1R,num1C] = find(Mat(:, :, j)==num1);
    [num2R,num2C] = find(Mat(:, :, j)==num2);
    temp = Mat(:, :, j);
    
    for i = 1:size(Mat(:, :, j),2)
        %If there exists no conduit starts or conduit ends in column i then the
        %construction matrix Mat should have all zeros in that column. This
        %overrides, for example, the conduit end inserted by the
        %generateConduits function at the end row and replaces it with the
        %value reset.
        if isempty(num1R(num1C == i)) || isempty(num2R(num2C == i))
            temp(:,i) = reset;
            continue
        elseif ~isempty(num1R(num1C == i))
            ind1 = num1R(num1C == i);
            %Find the first conduit start in that column
            ind1 = ind1(1);
            %Remove all conduit ends preceding the initial conduit start
            temp(1:ind1-1,i) = reset;
            ind2 = num2R(num2C == i);
            %Find the last conduit end in that column
            ind2 = ind2(end);
            %Remove all conduit starts after the final conduit end
            temp(ind2+1:end,i) = reset;
        end
    end
    MatOut(:, :, j) = temp;
end
%Note: Even if there are consecutive conduit starts or conduit ends, the
%generateConduits functions neglects the intermediate ones so there is no
%need to find and delete them.

end

