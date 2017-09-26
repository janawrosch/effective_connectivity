function [ correlation_matrix_norm ] = norm_correlation( correlation_matrix )
%Normalize correlation matrix
idx = find(~eye(size(correlation_matrix)));

corrvalues=correlation_matrix(idx);
corrvalues_sorted=sort(corrvalues, 'descend');

min_val=min(corrvalues);

correlation_matrix_norm=(correlation_matrix-min_val);


ten_max=corrvalues_sorted(1:10)-min_val;
max_val=mean(ten_max);


correlation_matrix_norm=correlation_matrix_norm./max_val;

end

