function [path_length_parameter]=path_length_calculations(distance_matrix_bin, distance_matrix, correlation_matrix, show_figures)

%% 1. Path length distribution and characteristic path length on binary matrix
path_length_values_bin=nonzeros(distance_matrix_bin);

cpl_bin=mean(path_length_values_bin(path_length_values_bin>0));
cpl_bin_std=std(path_length_values_bin(path_length_values_bin>0));
[pl_bin_n, pl_bin_x]=hist(path_length_values_bin,0:1:10);
pl_bin_h=kstest(pl_bin_n);

%% 2. Path length distribution and characteristic path length on weighted matrix
pl_temp=nonzeros(distance_matrix);
path_length_values=pl_temp(isfinite(pl_temp));

cpl=mean(path_length_values(path_length_values>0));
cpl_std=std(path_length_values(path_length_values>0));
[pl_wei_n, pl_wei_x]=hist(path_length_values,50);
pl_wei_h=kstest(pl_wei_n);

%% 3. Global efficiency
global_efficiency=efficiency_wei(correlation_matrix,0);
if isempty(global_efficiency)==1
    global_efficiency=0;
end


%% 4. Vulnerability
n_n=size(correlation_matrix,1);
for i=1:n_n
    correlation_matrix_removed=correlation_matrix;
    correlation_matrix_removed(i,:)=[];
    correlation_matrix_removed(:,i)=[];
    node_efficiency=efficiency_wei(correlation_matrix_removed,0);
    vulnerability(i,1)=(global_efficiency-node_efficiency)/global_efficiency;
end

vulnerability_max=max(vulnerability);
vulnerability_mean=mean(vulnerability);



%% 5. Show figures
if show_figures==1
figure
bar(pl_bin_x, pl_bin_n)
xlabel('Shortest path length binary directed')
ylabel('count')
title(sprintf('Charachteristic path length binary directed: %f; SD: %f', cpl_bin, cpl_bin_std))

figure
bar(pl_wei_x, pl_wei_n)
xlabel('Shortest path length weighted directed')
ylabel('count')
title(sprintf('Charachteristic path length weighted directed: %f; SD: %f', cpl, cpl_std))
end

%% 6. Prepare output

path_length_parameter.bin.collection=path_length_values_bin;
path_length_parameter.bin.cpl=cpl_bin;
path_length_parameter.bin.std=cpl_bin_std;
path_length_parameter.bin.distribution.x=pl_bin_x;
path_length_parameter.bin.distribution.n=pl_bin_n;
path_length_parameter.bin.distribution.h=pl_bin_h;

path_length_parameter.wei.collection=path_length_values;
path_length_parameter.wei.cpl=cpl;
path_length_parameter.wei.std=cpl_std;
path_length_parameter.wei.distribution.x=pl_wei_x;
path_length_parameter.wei.distribution.n=pl_wei_n;
path_length_parameter.wei.distribution.h=pl_wei_h;

path_length_parameter.global_efficiency=global_efficiency;

path_length_parameter.vulnerability.collection=vulnerability;
path_length_parameter.vulnerability.max=vulnerability_max;
path_length_parameter.vulnerability.mean=vulnerability_mean;




end