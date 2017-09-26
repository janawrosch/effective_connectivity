function [ counting_list ] = generating_countinglist(counting_list,proof_element,xy,x1x,x )
%UNTITLED Summary of this function goes here
%produces a countinglist of combinations of k array cells
%proof_element [elemt which should be looked for in counting_list]

not_in=0;
t=1;
% % for t=1:size(counting_list(:,2))
% % if isempty(counting_list)==1
% if size(counting_list,1)==1 && sum(counting_list)==0
%     counting_list(1,1)= proof_element;
%     counting_list(t,2)= counting_list(t,2)+1;
%     counting_list(t,4)= xy;
%     counting_list(t,7)= x1x;
%     counting_list(t,10)=x ;
%     return;
% end

% not_in = sum(counting_list(:,2)==proof_element) > 0;
% count_proof_element = sum(data()==proof_element) > 0;

while  t<=length(counting_list(:,2))
    if counting_list(t,1)~=proof_element%rxyx1
        t=t+1;
        not_in=1;
    else
        not_in=0;
        t=t+1;
        break;
    end
end
if not_in==0
    counting_list(t-1,2)=counting_list(t-1,2)+1;
else
    counting_list(t,1)= proof_element;
    counting_list(t,2)= counting_list(t,2)+1;
    counting_list(t,4)= xy;
    counting_list(t,7)= x1x;
    counting_list(t,10)=x ;%wordsr_k_x;
    
end

end

