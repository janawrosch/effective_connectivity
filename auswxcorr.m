function [correlation_matrix]=auswxcorr(data, windowsize)

%% Cross-correlation
dispstat('','init')

%% 1. Input

number_of_cells=size(data,1);
number_of_frames=size(data,2);
% windowsize=15;   %how many frames before and after the input spike are considered
correlation_matrix_nonorm=zeros(number_of_cells,number_of_cells);


%% 2. Xcorr der einzelnen Spur-Pärchen


for source=1:number_of_cells
    for target=1:number_of_cells
        if target~=source
            [r,lags]=xcorr(data(target,:), data(source,:),windowsize);
            [corr_value_temp]=cumprod(r(1,windowsize+1:2*windowsize+1));  
%             figure
%             plot(lags,r)
%             title(sprintf('Source %i, Target %i', source, target))
            correlation_matrix_nonorm(source,target)=corr_value_temp(1,end);

        end
    end
    dispstat(sprintf('Finished calculations for cell %i', source))
end


% 
% spike_frequency=sum(data,2)/number_of_frames;
% spikes_per_cell=sum(data,2);
% 
% for i=1:number_of_cells
%     for  j=1:number_of_cells
%         spike_frequency_matrix(i,j)=(spike_frequency(i)*spike_frequency(j))^2*1000;
%     end
% end

%% 3. Normierung der xcorr

for cell=1:number_of_cells
    power(cell)=max(xcorr(data(cell,:), data(cell,:)));
end

for source=1:number_of_cells
    for target=1:number_of_cells
        correlation_matrix_norm(source, target)=correlation_matrix_nonorm(source, target)/sqrt(power(source)*power(target));
    end
end


correlation_matrix=correlation_matrix_norm;
dispstat('Cross correlation reconstruction finished.', 'keepthis', 'timestamp')