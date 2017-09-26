close all
clearvars

show_control_figures=1; % show images for cell detection?                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          ;
saving=1; % save results?



filename_stem = '20170815_SAS_25_X';   % insert the stem of tif-file stack filenames
numb_images=31156;  % number of tif files
stimstapel=15; % In which tif-file stack is the electrical stimulation?
image_path='W:\psychiatrie\neurophotonik\2017\CONNECT_sulfasalazine\Networks\20170815_SAS_25\'; % data path to the image files. Make sure to include \ at the end
savepath='V:\AG_Neurophotonik\Projekte\CONNECT-Sulfasalazine\Data\Networks\'; % data path to where results should be saved
name='connect_sulfa_22';                           % file name for saved results


%%%%% Note additional input parameters in lines 47, 48 and 54!!!!!!



%________________________________________________________________________________________________________________________________________________%
%% 1. Cell detection

savename=sprintf('%s_rawdata.mat', name);
loadname=sprintf('%s_meanstack.mat', name);


dispstat('', 'init')
stimstapel_name=sprintf('%s%s%i.tif',image_path, filename_stem, stimstapel);

warning('off', 'all')
InfoImage=imfinfo(stimstapel_name);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
stapel_size=length(InfoImage);
stapel=zeros(nImage,mImage,stapel_size,'uint16');
TifLink=Tiff(stimstapel_name, 'r');


for bild=1:stapel_size
    TifLink.setDirectory(bild);
    stapel(:,:,bild) = TifLink.read();
end
TifLink.close();

% 
image_vorstim=mean(stapel(:,:,895:900),3); %% Adjust image stack layers according to the time point of electrical stimulation. Choose around 5 frames to be averaged before onset of the stimulation
image_stim=mean(stapel(:,:,1200:1205),3);  %% Choose about 5 frames to be averaged during/ after the electrical stimulation


diff_image=image_stim-image_vorstim;


[region_properties] = peakdetect_sbalzarini(diff_image,7,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cell detection: first number = cell radius, second number = intensity parameter
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         Trouble shooting:
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % If weak intensity cell bodies are not found try increasing the intensity parameter.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % If small cell bodies are not found, reduce the approximate cell radius.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % If many small “cells” are found within bigger bright spots in the recording, increase the cell radius.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % If many “cells” are found in the background noise, reduce the intensity parameter.
numb_regions=size(region_properties,1);
drawnow
dispstat('Regions detected', 'timestamp', 'keepthis')

if saving==1
    region_name=sprintf('%s_region_properties.mat', name);
    save([savepath region_name], 'region_properties', 'diff_image');
    dispstat('Region properties of detected regions saved','timestamp', 'keepthis')
end



%% 2. Read images
meanstack=zeros(numb_regions,numb_images);

for frame=1:stapel_size
    temp=stapel(:,:,frame);
    for region=1:numb_regions
        pixelstack(:,frame)=temp(region_properties(region,1).PixelIdxList);
        meanstack(region,31156-526-2042+frame)=mean(pixelstack(:,frame));
        clear pixelstack
    end
    switch frame
        case{2*round(stapel_size/10),  4*round(stapel_size/10),  6*round(stapel_size/10),  8*round(stapel_size/10), stapel_size}
            dispstat(sprintf('Finished extracting fluorescence traces from image %i of %i in the cell detection stack\n', frame, stapel_size), 'timestamp')
    end
end

warning('off', 'all')
for stapelnr=1:16
    if stapelnr==stimstapel
    else
        stapel_name=sprintf('%s%s%i.tif',image_path, filename_stem, stapelnr);
        
        InfoImage=imfinfo(stapel_name);
        mImage=InfoImage(1).Width;
        nImage=InfoImage(1).Height;
        stapel_size=length(InfoImage);
        stapel=zeros(nImage,mImage,stapel_size,'uint16');
        TifLink=Tiff(stapel_name, 'r');
        
        
        for bild=1:stapel_size
            TifLink.setDirectory(bild);
            stapel(:,:,bild) = TifLink.read();
            switch bild
                case{2*round(stapel_size/10),  4*round(stapel_size/10),  6*round(stapel_size/10),  8*round(stapel_size/10), stapel_size}
                    dispstat(sprintf('Finished loading image %i of %i in stack %i\n', frame, stapel_size, stapelnr), 'timestamp')
            end
        end
        TifLink.close();
        
        for frame=1:stapel_size
            temp=stapel(:,:,frame);
            for region=1:numb_regions
                pixelstack(:,1)=temp(region_properties(region,1).PixelIdxList);
                meanstack(region,(frame-2042+stapelnr*2042))=mean(pixelstack(:,1));
                clear pixelstack
            end
            switch frame
                case{2*round(stapel_size/10),  4*round(stapel_size/10),  6*round(stapel_size/10),  8*round(stapel_size/10), stapel_size}
                    dispstat(sprintf('Fluorescence trace extracted for image %i of %i in stack %i\n', frame, stapel_size, stapelnr), 'timestamp')
            end
        end
    end
end

warning('on', 'all')
dispstat('Fluorescence traces extracted', 'timestamp', 'keepthis')



if saving==1
    firstname=sprintf('%s_meanstack.mat', name);
    save([savepath firstname], 'meanstack');
    dispstat(sprintf('Meanstack of %s saved as %s', filename_stem, name), 'timestamp', 'keepthis')
end


if show_control_figures==1
    mean_trace=mean(meanstack,1);
    
    figure
    hold on
    for i=1:numb_regions
        plot(meanstack(i,:), 'Color', rand(1,3))
    end
    plot(mean_trace, 'Color', 'black', 'LineWidth', 2)
    title('Traces of all cells and mean trace')
    
    
    
    figure
    hold on
    for i=1:numb_regions
        plot(smooth(meanstack(i,:),55), 'Color', rand(1,3))
    end
    plot(mean_trace, 'Color', 'black', 'LineWidth', 2)
    title('Smoothened traces of all cells and mean trace')
end

% load chirp, sound(y,Fs)