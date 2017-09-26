function [ estimated_spikes ] = spike_estimation_th2( fluorescence_data, show_figures )
%spike estimation with our own method.
% 2 criteria for a spike:
%     1. derivate of the fluorescence curve must exceed a certain threshold (mean+n*sigma)
%     2. spike is located at the local maximum of the peak
% Parameter input:
moving_average_span=40;
noise_std_factor=2;
sprungbreite=100;
maxbreite=60;


number_of_cells=size(fluorescence_data,1);
number_of_images=size(fluorescence_data,2);

% Data smoothing and noise level calculation
[noise] = smooth_noise( fluorescence_data, moving_average_span, number_of_cells,noise_std_factor, show_figures); 
% Peak detection
[ estimated_spikes ] = peakdetection_threshold2( fluorescence_data, noise, number_of_images, number_of_cells, sprungbreite, maxbreite, show_figures);


end

