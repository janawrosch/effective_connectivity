function [pp_parameter]=pp_analysis(correlation_matrix, correlation_matrix_bin, number_of_connections, number_of_cells, connectivity_degree, micrometer_per_pixel, region_properties, show_figures)

idx = find(correlation_matrix(:)>0);
%% 1. Propagation probability (=correlation-value) mean and distribution
[n_pp, x_pp]=hist(correlation_matrix(idx), ceil(number_of_connections/3));
pp_h=kstest(n_pp);

pp_values=correlation_matrix(idx);
mean_pp=mean(pp_values);
std_pp=std(pp_values);


%% 2. physical distance 

for source=1:number_of_cells
    for target=1:number_of_cells
        pd_matrix(source, target)=abs(micrometer_per_pixel*(sqrt( ( region_properties(target,1).Centroid(1,1)  -  region_properties(source,1).Centroid(1,1)  )^2+( region_properties(target,1).Centroid(1,2)  -  region_properties(source,1).Centroid(1,2)  )^2)));
    end
end


pd_values=pd_matrix(idx);
pd_values_all=pd_matrix(:);
[n_pd, x_pd]=hist(pd_values, ceil(number_of_connections/10));
[n_pd_all, x_pd_all]=hist(pd_values_all, ceil(number_of_cells*number_of_cells/10));
pd_h=kstest(n_pd);
pd_all_h=kstest(n_pd_all);

mean_pd=mean(pd_values);
mean_pd_all=mean(pd_values_all);
std_pd=std(pd_values);
std_pd_all=std(pd_values_all);

%% 3. Physical distance vs. propgations propability

[r,p]=corrcoef(pd_values, nonzeros(correlation_matrix(:)));

%% 4. Figures

if show_figures==1
figure
bar(x_pp, n_pp,'blue')
xlabel('correlation value (propagation propability)')
ylabel('count')
title(sprintf('Mean corr value: %f; SD: %f', mean_pp, std_pp))

figure
bar(x_pd, n_pd,  'blue')
xlabel('physical distance of connected cells [µm]')
ylabel('count')
title(sprintf('Mean distance: %fµm; SD: %fµm', mean_pd, std_pd))

figure
bar(x_pd_all, n_pd_all,  'blue')
xlabel('physical distance of all cells [µm]')
ylabel('count')
title(sprintf('Mean distance: %fµm; SD: %fµm', mean_pd_all, std_pd_all))

% figure
% scatter(nonzeros(pd_matrix_connected(:)), nonzeros(correlation_matrix(:)))
% xlabel('physical distance [µm]')
% ylabel('propagation propability')
% title(sprintf('linear correlation p: %f', p(1,2)))

end

%% 5. Prepare Output

pp_parameter.pp.collection=correlation_matrix(idx);
pp_parameter.pp.mean=mean_pp;
pp_parameter.pp.std=std_pp;
pp_parameter.pp.distribution.n=n_pp;
pp_parameter.pp.distribution.x=x_pp;
pp_parameter.pp.distribution.h=pp_h;

pp_parameter.pd.collection=pd_values;
pp_parameter.pd.mean=mean_pd;
pp_parameter.pd.std=std_pd;
pp_parameter.pd.distribution.n=n_pd;
pp_parameter.pd.distribution.x=x_pd;
pp_parameter.pd.distribution.h=pd_h;

pp_parameter.pd_all.collection=pd_values_all;
pp_parameter.pd_all.mean=mean_pd_all;
pp_parameter.pd_all.std=std_pd_all;
pp_parameter.pd_all.distribution.n=n_pd_all;
pp_parameter.pd_all.distribution.x=x_pd_all;
pp_parameter.pd_all.distribution.h=pd_all_h;

pp_parameter.connectivity_degree=connectivity_degree;


end