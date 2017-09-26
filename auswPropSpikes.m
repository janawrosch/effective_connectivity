function [ correlation_matrix ] = auswPropSpikes( traces )
%Counts how many spikes are propagated from one time frame to the next for
%all possible pairs of cells

n_n=size(traces,1);
correlation_matrix=zeros(n_n,n_n);
common_spikes=0;

for cellA=1:n_n
    for cellB=1:n_n
        if cellA~=cellB
            common_spikes=0;
            for frame=1:size(traces,2)-2 % go through all the frames
                if traces(cellA,frame)==1 % if the source cell has a spike
                    if sum(traces(cellB, frame+1))>0 % check if the target cell has a spike in the next time frame; To account for same bin interaction change to frame:frame+1
                        common_spikes(1,frame)=1;
                    end
                end
            end
            correlation_matrix(cellA,cellB)=sum(common_spikes);
        end
    end
end
end

