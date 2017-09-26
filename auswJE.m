function [correlation_matrix]=auswJE(spike_trains)


%% Bin size: which cISIs count as "the same"
bin_size=1; % 0.1ms von Garafalo et al.
k=3; % 3 bins


%% Joint Entropy
dispstat('','init')
%% 1. Input
number_of_cells=size(spike_trains,1);

correlation_matrix=zeros(number_of_cells,number_of_cells);


%% 2. Interspike Intervall berechnen

JE=zeros(number_of_cells, number_of_cells);

cISI2=zeros(number_of_cells, number_of_cells);
for source=1:number_of_cells
    for target=1:number_of_cells
        if source~=target
            cISI_count=0;
            for ref_spike=1:size(spike_trains{source,1},1)-1
                rf_st=spike_trains{source,1}(ref_spike,1);
                target_st=spike_trains{target,1}(find(spike_trains{target,1}>=rf_st, 1),1);
                if isfinite(target_st)  %<=spike_trains{source,1}(ref_spike+1,1);
                    cISI_count=cISI_count+1;
                    cISI2(source, target,cISI_count)=target_st-rf_st;
                end
            end
            cISI2(source, target,:)=sort(cISI2(source, target,:));
            start=find(cISI2(source, target,:)>0,1);
            cISI_step=start;
            while cISI_step<size(cISI2(source, target,:),3)
                for k_step=1:k
                    anzahl=sum(cISI2(source, target,cISI2(source, target,:)>=cISI2(source, target, cISI_step))<=cISI2(source, target, cISI_step)+k*bin_size);
                    p_cISI_step=anzahl/size(find(cISI2(source, target,:)>0),1);
                    JE(source, target)=JE(source, target)-(p_cISI_step*log2(p_cISI_step));
                end
                cISI_step=cISI_step+anzahl+1;
            end
            
        end
    end
    dispstat(sprintf('JE calculated for cell %i', source), 'timestamp')
end









%    for connection=1:number_of_connections
%        JE(connection)=0;
%        cISI_row=cISI(connection,:);
%        cISI_row_clean=cISI_row(cISI_row>0);
%        cISI_min=min(cISI_row_clean);
%        cISI_max=max(cISI_row_clean);
%        for k=cISI_min:cISI_max
%            p_cISI_k=(size(find(cISI_row_clean==k),2))/(size(cISI_row_clean,2));
%            if p_cISI_k>0
%            JE(connection)=JE(connection)+(p_cISI_k*log2(p_cISI_k));
%            end
%        end
%    end
%
%    for connection=1:number_of_connections
%        correlation_matrix(which_cells(connection,1), which_cells(connection,2))=-1*JE(connection);
%    end

%% Turn around the order of values to match the other algorithms

maxcorr=max(JE(:));
correlation_matrix_reordered=maxcorr-JE;
for cell=1:size(correlation_matrix_reordered,1)
    correlation_matrix_reordered(cell,cell)=0;
end



correlation_matrix=correlation_matrix_reordered;

dispstat('Joint Entropy reconstruction finished.', 'keepthis', 'timestamp')