function [ counting_list ] = counting_sorted_list_xy( counting_list,sorted_col, col_to_write_in )

%UNTITLED2 Summary of this function goes here
% This functions counts how often one element appears in a sorted list.

z=1;
proof = counting_list(1,sorted_col); % proof is the first element of the previously sorted column
while z<=size(counting_list,1)
    count =0; n=0;
    while z>size(counting_list,1)==1 || counting_list(z,sorted_col)==proof; % If proof is in row=z in countinglist 
        if z>size(counting_list,1)==1
            break;
        end
        count = count+1*counting_list(z,2);
        n=n+1;
        z=z+1;
    end
    counting_list(z-n:z-1,col_to_write_in)=count;  % all rows containing the proof element obtain the number of counts
    if z>size(counting_list,1)==1 
        return;
    end
    proof=counting_list(z,sorted_col); % the new proof element is the next element in the list 
end

end

