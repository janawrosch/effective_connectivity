function [noise] = smooth_noise( traces_raw, moving_average_span, number_of_cells,noise_std_factor, show_figures)


% Work with relative fluorescence

for cell=1:number_of_cells
    [counts, centers]=hist(traces_raw(cell,:),1000);
    fit_data=fit(centers', counts', 'gauss1');
    fit_mean=fit_data.b1;
    fit_sigma=fit_data.c1/sqrt(2);
    
    %% Threshold = mean+ n*sigma
    noise(cell,1)=fit_mean+noise_std_factor*fit_sigma; 



end


