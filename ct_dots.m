close all
clear all

%% load folder and images

folder = 'thesis_data/CT_0_6MM_RETAKE';
[image_ct,info_ct] = dicomfolder(folder);
%%
dimx = round(info_ct{1}.PixelSpacing(1)*size(image_ct,1));
dimy = round(info_ct{1}.PixelSpacing(2)*size(image_ct,2));
dimz = round(info_ct{1}.SliceThickness*size(image_ct,3));

%resize voxels to be 1mm x 1mm x 1mm
slices_ct = imresize3(image_ct,[dimx,dimy,dimz]);

slices_ct = slices_ct(35:294,:,72:336); 

%make 3d template of grid cross
dim = 23; 
middle = round(dim/2);
mask_3 = zeros(dim,dim,dim);
mask_3(:,middle-1:middle+1,middle) = 1;
mask_3(:,middle,middle-1:middle+1) = 1;
mask_3(middle,:,middle-1:middle+1) = 1;
mask_3(middle-1:middle+1,:,middle) = 1;
mask_3(middle,middle-1:middle+1,:) = 1;
mask_3(middle-1:middle+1,middle,:) = 1;

%% crosscorr in 3 dimensions to find crosses (time intensive)
cross_corr_3 = normxcorr3(mask_3,slices_ct,'same'); 
Mcc = max(max(max(cross_corr_3)));
cross_corr_3 = cross_corr_3 / Mcc;

% filter
filter_spot=zeros(size(cross_corr_3));
filter_spot(cross_corr_3>0.4)=1;
   
%% find centroids
CC = bwconncomp(filter_spot);
S = regionprops(CC,'Centroid');
centroids = cat(1, S.Centroid); %x,y,z of centroids, (cat struct->array) 

central = round(size(slices_ct,3)/2);

%display images around center to find central slice
figure
dim=size(slices_ct(:,:,central-3:central+3));
M=max(max(max(slices_ct)));
im_mount=reshape(slices_ct(:,:,central-3:central+3),[dim(1) dim(2) 1 dim(3)]);
montage(im_mount/M)

central_z = input('enter the index of the middle slice (most prominent grid): ')+central - 3;

% choose control points to determine middle in xy plane
circ = [ones(3,8);[ones(1,3) 0 0 ones(1,3)];[1 1 0 0 0 0 1 1];[ones(1,3) 0 0 ones(1,3)];ones(3,8)];
%choose centre, save fixed point to fixedPoints:
cpselect(circ,slices_ct(:,:,central_z)); 

%% mod filter coordinates (run this section again after selecting the centre)
%(only taking points that are close to expected coordinates into filtered)
iso = [(fixedPoints) central_z]; %iso centre defined by usier input
delta = (centroids-iso);

delta_switch(:,1) = delta(:,3);
delta_switch(:,2) = delta(:,2);
delta_switch(:,3) = -delta(:,1);
delta = delta_switch;



mod_z = mod(delta,20);
filtered = [];
for n = 1:size(centroids,1)
            if abs(mod_z(n,3)-20)<3 || mod_z(n,3)<3
                filtered = cat(1,filtered,delta(n,:));
            end        
end


% normal grid

load grid_norm
grid_norm = grid_norm(146:end,:);
%
[Idx] = knnsearch(filtered,grid_norm);
sorted_ct = filtered(Idx,:);
sorted_ct = sorted_ct - sorted_ct(778+145,:); %adjust middle to 0,0,0

% registering using pointCloud objects to use CV toolbox
pc_ct = pointCloud(sorted_ct);
pc_norm = pointCloud(grid_norm);

[tform,ct_reg] = pcregistericp(pc_ct,pc_norm);

[~,D] = knnsearch(ct_reg.Location,grid_norm); %calculate distance of points

D_red=D;
D_red(D>10) = -1; 
D = D_red;
 
map = jet;

m = 63/max(D);
b = 1;
%
figure; hold on; title('Distortion in mm')
pcshow(pointCloud([0 0 0]));
colormap jet; colorbar;caxis([0 max(D)]);
colorbar('Ticks',[0,1,max(D)],...
         'TickLabels',{0,1,1.8385},'Color','w')
for point = 1:size(D)
    if D(point)>=0
    color = round(m*D(point)+b);
    plot3(grid_norm(point,1),grid_norm(point,3),-grid_norm(point,2),'+','Color',map(color,:),'MarkerSize',8,'LineWidth',2)
    end
end

% dist. from iso center vs distortion

[~,D_centre] = knnsearch([0 0 0],grid_norm);
D_10=zeros(size(D_centre));
D_10(D_centre<=100)=D(D_centre<=100);
'maximal distortion in r = 10cm: ', max(D_10)

D_5=zeros(size(D_centre));
D_5(D_centre<=50)=D(D_centre<=50);
'maximal distortion in r = 5cm: ', max(D_5)

D_20=zeros(size(D_centre));
D_20(D_centre<=200)=D(D_centre<=200);
'maximal distortion in r = 20cm: ', max(D_20)



% Scatter Plot
figure;
scatter(D_centre(D~=-1),D(D~=-1),'x')
title('Distortion in CT');
xlabel('Distance from centre in mm');ylabel('Distortion in mm');
grid on;






    


