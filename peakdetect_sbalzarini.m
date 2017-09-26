function [region_properties] = peakdetect_sbalzarini(img, w, pth)


figure
imagesc(img)
axis xy
title('Search image')


% Correlation length of camera noise (usu. set to unity)
lambdan = 1;

% Normalize image:
img = double(img);
img = (img-min(img(:)))/(max(img(:))-min(img(:)));

% some often used quantities
idx = [-w:1:w];     % index vector
dm = 2*w+1;         % diameter
im = repmat(idx',1,dm);
jm = repmat(idx,dm,1);
imjm2 = im.^2+jm.^2;
siz = size(img);   % image size

%======================================================================
% STEP 1: Image restoration
%======================================================================

% build kernel K for background extraction and noise removal
% (eq. [4])
B = sum(exp(-(idx.^2/(4*lambdan^2))));
B = B^2;
K0 = 1/B*sum(exp(-(idx.^2/(2*lambdan^2))))^2-(B/(dm^2));
K = (exp(-(imjm2/(4*lambdan^2)))/B-(1/(dm^2)))/K0;

% apply convolution filter
filtered = conv2(img,K,'same');



%======================================================================
% STEP 2: Locating particles
%======================================================================

% determining upper pth-th percentile of intensity values
pth = 0.01*pth;
% das gefilterte bild ist genauso groﬂe wie das eigentlich bild, d.h. der
% filter wird auch am bildrand angewendet. da kommen aber zwangsweise
% unbrauchbare werte raus.
[m,n]=size(filtered);


[cnts,bins] = imhist(filtered(1+w:m-w ,1+w:n-w));
l = length(cnts);
k = 1;
while sum(cnts(l-k:l))/sum(cnts) < pth,
    k = k + 1;
end;
thresh = bins(l-k+1);



% generate circular mask of radius w
mask = zeros(dm,dm);
mask(find(imjm2 <= w*w)) = 1;

% identify individual particles as local maxima in a
% w-neighborhood that are larger than thresh
dil = imdilate(filtered,mask);
[Rp,Cp] = find((dil-filtered)==0);
temp_r=Rp;
temp_c=Cp;
count=1;
clearvars Rp Cp
for cand=1:size(temp_r,1)
if temp_r(cand,1)>w && temp_c(cand,1)>w && temp_r(cand,1)<size(img,2)-w && temp_c(cand,1)<size(img,1)-w % Cut off the edges
Rp(count,1)=temp_r(cand,1);
Cp(count,1)=temp_c(cand,1);
count=count+1;
end
end



particles = zeros(siz);
V = find(filtered(sub2ind(siz,Rp,Cp))>thresh);
R = Rp(V);
C = Cp(V);
particles(sub2ind(siz,R,C)) = 1;

[columnsInImage rowsInImage] = meshgrid(1:size(img,2), 1:size(img,1));
radius=7;

for region=1:size(R,1)
region_properties(region,1).Centroid(1,1)=C(region,1);
region_properties(region,1).Centroid(1,2)=R(region,1);
centerX=C(1);
centerY=R(1);
region_properties(region,1).Circle=(rowsInImage - region_properties(region,1).Centroid(1,2)).^2 + (columnsInImage - region_properties(region,1).Centroid(1,1)).^2 <= radius.^2;
region_properties(region,1).PixelIdxList=find(region_properties(region,1).Circle>0);
region_properties(region,1).mittlereHelligkeit=mean(img(region_properties(region,1).PixelIdxList(:)));
end



filter=fspecial('average',5);
img_f=imfilter(img,filter);




figure('units','normalized','outerposition',[0 0 1 1])
imagesc(img_f)
colormap('jet')
axis xy
hold on
scatter(C,R, 'filled', 'MarkerFaceColor', 'black')
title('moving average filtered image with found regions')

control_img=zeros(size(img));
for region=1:size(region_properties,1)
control_img(region_properties(region,1).PixelIdxList(:))=region;
end
figure
subplot(1,2,1)
imagesc(img)
% axis xy
title('Search image')
subplot(1,2,2)
imagesc(control_img)
% axis xy
title(sprintf('Found regions: %i', size(region_properties,1)))

figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
imagesc(img_f)
% axis xy
title('Smoothed search image')
subplot(1,2,2)
imagesc(control_img)
% axis xy
title(sprintf('Found regions: %i', size(region_properties,1)))

end