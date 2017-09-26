% close all


rec_nr_start=13; % First recording to be analyzed when running batches
rec_nr_stop=22; % Last recording to be analyzed when running batches
data_path='V:\AG_Neurophotonik\Projekte\CONNECT-Sulfasalazine\Data\Networks\'; % Data path to find the results of part B
show_figures=0; % Show figures?
clear_data=1; %clear simulation data after each run (to save memory); will always be on if more than 20 simulations are batch processed
data_filename='connect_sulfa'; % filename stem of restuls of part B

data=[];
if rec_nr_stop-rec_nr_start>20
    clear_data=1;
    show_figures=0;
end




dispstat('','init')
%%
for rec_nr=rec_nr_start:rec_nr_stop
    dispstat(sprintf('Loading data for recording %i',rec_nr), 'keepthis', 'timestamp')
    filename_spikes=sprintf('%s%s_%i_estimated_spikes.mat',data_filename, data_path, rec_nr);
    
    spikes_temp=load(filename_spikes);
    data(rec_nr,1).estimated_spikes=spikes_temp;
    
    dispstat('Data passed to network reconstruction', 'keepthis', 'timestamp')
    %% Network reconstruction with the different algorithms and analysis of accuracy
    [data(rec_nr,1).xcorr.correlation_matrix]=auswxcorr(data(rec_nr,1).estimated_spikes.raster_plot, 1);
    [data(rec_nr,1).xcorr.correlation_matrix_norm]=norm_correlation(data(rec_nr,1).xcorr.correlation_matrix);
    
    [data(rec_nr,1).MI.correlation_matrix]=auswMI(data(rec_nr,1).estimated_spikes.raster_plot, 1,1);
    [data(rec_nr,1).MI.correlation_matrix_norm]=norm_correlation(data(rec_nr,1).MI.correlation_matrix);
    
    [data(rec_nr,1).JE.correlation_matrix]=auswJE(data(rec_nr,1).estimated_spikes.spike_trains); 
    [data(rec_nr,1).JE.correlation_matrix_norm]=norm_correlation(data(rec_nr,1).JE.correlation_matrix);
    
    [data(rec_nr,1).TE.correlation_matrix]=auswTE16(data(rec_nr,1).estimated_spikes.raster_plot);
    [data(rec_nr,1).TE.correlation_matrix_norm]=norm_correlation(data(rec_nr,1).TE.correlation_matrix);
    
    [data(rec_nr,1).GTE.correlation_matrix]=auswGTE16(data(rec_nr,1).estimated_spikes.raster_plot);
    [data(rec_nr,1).GTE.correlation_matrix_norm]=norm_correlation(data(rec_nr,1).GTE.correlation_matrix);
    
    [data(rec_nr,1).SC.correlation_matrix] = auswPropSpikes(data(rec_nr,1).estimated_spikes.raster_plot);
    [data(rec_nr,1).SC.correlation_matrix_norm]=norm_correlation(data(rec_nr,1).SC.correlation_matrix);
   
    [data(rec_nr,1).PP.correlation_matrix] = auswPropSpikes_ratio(data(rec_nr,1).estimated_spikes.raster_plot);    
    [data(rec_nr,1).PP.correlation_matrix_norm]=norm_correlation(data(rec_nr,1).PP.correlation_matrix);
    
    [data(rec_nr,1).PP_samebin.correlation_matrix] = auswPropSpikes_ratio_samebin(data(rec_nr,1).estimated_spikes.raster_plot); 
    [data(rec_nr,1).PP_samebin.correlation_matrix_norm]=norm_correlation(data(rec_nr,1).PP_samebin.correlation_matrix);
    
    
    
    network_xcorr=data(rec_nr,1).xcorr.correlation_matrix_norm;
    network_MI=data(rec_nr,1).MI.correlation_matrix_norm;
    network_JE=data(rec_nr,1).JE.correlation_matrix_norm;
    network_TE=data(rec_nr,1).TE.correlation_matrix_norm;
    network_GTE=data(rec_nr,1).GTE.correlation_matrix_norm;
    network_SC=data(rec_nr,1).SC.correlation_matrix_norm;
    network_PP=data(rec_nr,1).PP.correlation_matrix_norm;
    network_PP_samebin=data(rec_nr,1).PP_samebin.correlation_matrix_norm;
    
    save(sprintf('%s%s_%i_reconstructed_network_norm.mat',data_filename ,data_path, rec_nr),  'network_xcorr', 'network_MI', 'network_JE','network_TE','network_GTE','network_SC', 'network_PP', 'network_PP_samebin');
    
    network_xcorr=data(rec_nr,1).xcorr.correlation_matrix;
    network_MI=data(rec_nr,1).MI.correlation_matrix;
    network_JE=data(rec_nr,1).JE.correlation_matrix;
    network_TE=data(rec_nr,1).TE.correlation_matrix;
    network_GTE=data(rec_nr,1).GTE.correlation_matrix;
    network_SC=data(rec_nr,1).SC.correlation_matrix;
    network_PP=data(rec_nr,1).PP.correlation_matrix;
    network_PP_samebin=data(rec_nr,1).PP_samebin.correlation_matrix;
    
    save(sprintf('%s%s_%i_reconstructed_network.mat', data_filename,data_path, rec_nr),  'network_xcorr', 'network_MI', 'network_JE','network_TE','network_GTE','network_SC', 'network_PP', 'network_PP_samebin');
    
    
    
    
    
    if clear_data==1
        data(rec_nr,:)=[];
    end
    
    
    dispstat(sprintf('Network reconstrucion completed for recording %i (Currently running recording %i through %i)\n', rec_nr, rec_nr_start, rec_nr_stop), 'keepthis', 'timestamp')
    
    
end





