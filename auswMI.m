function [correlation_matrix]=auswMI(data, bin_size, time_shift)

%% Mutual Information
dispstat('','init')
%% 1. Input

% bin_size=5;
% time_shift=3;



number_of_cells=size(data,1);
number_of_frames=size(data,2);
correlation_matrix=zeros(number_of_cells,number_of_cells);
n_bins=ceil(number_of_frames/bin_size);



%% 2. Breakup time series into bins and count spike occurence

spike_count=zeros(number_of_cells,n_bins);

bin=1;
for frame=1:number_of_frames-bin_size
    spike_count(:,bin)=sum(data(:,frame:frame+bin_size),2);
    bin=bin+1;
end

maxbincount=max(spike_count(:));
dispstat('Spike count per bin generated')




%% 3. Calculate mutual information
mutual_information=zeros(number_of_cells, number_of_cells);

for source=1:number_of_cells
    for target=source+1:number_of_cells
        for source_count=1:maxbincount
            occurences_source_count=find(spike_count(source,:)==source_count);
            p_source_count=size(occurences_source_count,2)/(n_bins);
            for target_count=1:maxbincount
                occurences_target_count=find(spike_count(target,:)==target_count);
                p_target_count=size(occurences_target_count,2)/(n_bins);
                occurences_target_count_shifted_pos=occurences_target_count-(ones(size(occurences_target_count)).*time_shift);
                occurences_target_count_shifted_neg=occurences_target_count+(ones(size(occurences_target_count)).*time_shift);
                joint_occurences_pos=intersect(occurences_source_count, occurences_target_count_shifted_pos);
                joint_occurences_neg=intersect(occurences_source_count, occurences_target_count_shifted_neg);
                p_joint_pos=size(joint_occurences_pos,2)/(n_bins-time_shift);
                p_joint_neg=size(joint_occurences_neg,2)/(n_bins-time_shift);
                MI_to_add_pos=(p_joint_pos*log2(p_joint_pos/(p_source_count*p_target_count)));
                MI_to_add_neg=(p_joint_neg*log2(p_joint_neg/(p_source_count*p_target_count)));
                if isfinite(MI_to_add_pos)
                    mutual_information(source, target)=mutual_information(source, target)+MI_to_add_pos;
                end
                if isfinite(MI_to_add_neg)
                    mutual_information(target, source)=mutual_information(target, source)+MI_to_add_neg;
                end
            end
        end
        
    end
    
    dispstat(sprintf('Mutual information calculated for cell %i', source))
    
end


correlation_matrix=mutual_information;

dispstat('Mutual information reconstruction finished.', 'keepthis', 'timestamp')



