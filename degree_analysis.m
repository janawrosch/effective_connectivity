function [degree_parameter]=degree_analysis(number_of_cells, correlation_matrix_bin, region_properties, show_figures)
%% 1. In and out degree and their distributions based on correlation_matrix_bin
degree_number=zeros(number_of_cells,3);
for cell=1:number_of_cells;
    degree_number(cell,2)=sum(correlation_matrix_bin(cell,:)); %out-degree
    degree_number(cell,1)=sum(correlation_matrix_bin(:,cell)); %in-degree
    if degree_number(cell,1)>0 && degree_number(cell,2)>0
    degree_number(cell,3)=degree_number(cell,2)/degree_number(cell,1); % out per in
    end
end
in_degree_mean=mean(degree_number(degree_number(:,1)>0,1));
in_degree_std=std(degree_number(degree_number(:,1)>0,1));
in_degree_norm=degree_number(:,1)./number_of_cells;
out_degree_mean=mean(degree_number(degree_number(:,2)>0,2));
out_degree_std=std(degree_number(degree_number(:,2)>0,2));
out_degree_norm=degree_number(:,2)./number_of_cells;
opi_degree_mean=mean(degree_number(degree_number(:,3)>0,3));
opi_degree_std=std(degree_number(degree_number(:,3)>0,3));


[in_degree_n, in_degree_x]=hist(degree_number(degree_number(:,1)>0,1),0:max(degree_number(:,1:2)));
[out_degree_n, out_degree_x]=hist(degree_number(degree_number(:,2)>0,2),0:max(degree_number(:,1:2)));
[opi_degree_n, opi_degree_x]=hist(degree_number(degree_number(:,3)>0,3),10);

in_degree_h=kstest(in_degree_n);
out_degree_h=kstest(out_degree_n);
opi_degree_h=kstest(opi_degree_n);

if in_degree_h==1
    norm_in='no';
elseif in_degree_h==0
    norm_in='yes';
end
if out_degree_h==1
    norm_out='no';
elseif out_degree_h==0
    norm_out='yes';
end
if opi_degree_h==1
    norm_opi='no';
elseif opi_degree_h==0
    norm_opi='yes';
end
%% 2. Check for correlation with image properties
area=zeros(number_of_cells,1);
helligkeit=zeros(number_of_cells, 1);
for cell=1:number_of_cells
    area(cell,1)=size(region_properties(cell,1).PixelIdxList,1);
    helligkeit(cell,1)=region_properties(cell,1).mittlereHelligkeit(1,1);
end

[r_in_area, p_in_area]=corrcoef(degree_number(:,1), area(:,1));
[r_in_hell, p_in_hell]=corrcoef(degree_number(:,1), helligkeit(:,1));
[r_out_area, p_out_area]=corrcoef(degree_number(:,2), area(:,1));
[r_out_hell, p_out_hell]=corrcoef(degree_number(:,2), helligkeit(:,1));
% [r_in_sum, p_in_sum]=corrcoef(degree_number(:,1), spike_count(:,1));



%% 3. Show figures
if show_figures==1
figure
subplot(1,2,1)
bar(in_degree_x, in_degree_n, 'red')
xlabel('in-degree')
ylabel('count')
title(sprintf('Median in degree: %f; SD: %f\nDistribution normal: %s', in_degree_mean, in_degree_std, norm_in))

subplot(1,2,2)
bar(out_degree_x, out_degree_n,  'blue')
xlabel('out-degree')
ylabel('count')
title(sprintf('Median out degree: %f; SD: %f\nDistribution normal: %s',  out_degree_mean, out_degree_std, norm_out))

figure
bar(opi_degree_x, opi_degree_n,  'blue')
title(sprintf('Median out per in degree: %f; SD: %f; Distribution normal: %s', opi_degree_mean, opi_degree_std, norm_opi))
xlabel('out per in degree')
ylabel('count')

figure
subplot(1,2,1)
scatter(degree_number(:,1), area(:,1))
xlabel('in degree')
ylabel('region area')
title(sprintf('linear correlation p: %f', p_in_area(1,2)))
subplot(1,2,2)
scatter(degree_number(:,2), area(:,1))
xlabel('out degree')
ylabel('region area')
title(sprintf('linear correlation p: %f', p_out_area(1,2)))


figure
subplot(1,2,1)
scatter(degree_number(:,1), helligkeit(:,1))
xlabel('in degree')
ylabel('region mean intensity')
title(sprintf('linear correlation p: %f', p_in_hell(1,2)))
subplot(1,2,2)
scatter(degree_number(:,2), helligkeit(:,1))
xlabel('out degree')
ylabel('region mean intensity')
title(sprintf('linear correlation p: %f', p_out_hell(1,2)))

% figure
% subplot(1,1,1)
% scatter(degree_number(:,1), spike_count(:,1))
% xlabel('out degree')
% ylabel('region area')
% title(sprintf('linear correlation p: %f', p_in_sum(1,2)))



end

%% 4. Prepare output

degree_parameter.in_degree.mean=in_degree_mean;
degree_parameter.in_degree.std=in_degree_std;
degree_parameter.in_degree.collection=degree_number(:,1);
degree_parameter.in_degree.norm=in_degree_norm;
degree_parameter.in_degree.distribution.n=in_degree_n;
degree_parameter.in_degree.distribution.x=in_degree_x;
degree_parameter.in_degree.distribution.h=in_degree_h;
degree_parameter.out_degree.mean=out_degree_mean;
degree_parameter.out_degree.std=out_degree_std;
degree_parameter.out_degree.collection=degree_number(:,2);
degree_parameter.out_degree.norm=out_degree_norm;
degree_parameter.out_degree.distribution.n=out_degree_n;
degree_parameter.out_degree.distribution.x=out_degree_x;
degree_parameter.out_degree.distribution.h=out_degree_h;
degree_parameter.opi_degree.mean=opi_degree_mean;
degree_parameter.opi_degree.std=opi_degree_std;
degree_parameter.opi_degree.collection=degree_number(:,3);
degree_parameter.opi_degree.distribution.n=opi_degree_n;
degree_parameter.opi_degree.distribution.x=opi_degree_x;
degree_parameter.opi_degree.distribution.h=opi_degree_h;

end

