close all
clear all

%% load folder and images

folder = 'CT_STEREO_SCHICHT_0_6MM_0003';
[image_ct,info_ct] = dicomfolder(folder);
%%
dimx = round(info_ct{1}.PixelSpacing(1)*size(image_ct,1));
dimy = round(info_ct{1}.PixelSpacing(2)*size(image_ct,2));
dimz = round(info_ct{1}.SliceThickness*size(image_ct,3));

%resize voxels to be 1mm x 1mm x 1mm
slices_ct = imresize3(image_ct,[dimx,dimy,dimz]);

slices_ct = slices_ct(35:294,:,72:336); %to generalize,insert option to set these parameters


%make 3d template of grid cross
dim = 23; %has to be uneven, >=23 so there are no artefacts in between
middle = round(dim/2);
mask_3 = zeros(dim,dim,dim); %make 3d template of grid cross
mask_3(:,middle-1:middle+1,middle) = 1;
mask_3(:,middle,middle-1:middle+1) = 1;
mask_3(middle,:,middle-1:middle+1) = 1;
mask_3(middle-1:middle+1,:,middle) = 1;
mask_3(middle,middle-1:middle+1,:) = 1;
mask_3(middle-1:middle+1,middle,:) = 1;

%% crosscorr in 3 dimensions to find crosses (time intensive)
cross_corr_3 = normxcorr3(mask_3,slices_ct,'same'); 

% normalize
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

mod_xy = mod(delta,20); %not needed for ct so far
mod_z = mod(delta,20);
filtered = [];
for n = 1:size(centroids,1)
%     if abs(mod_xy(n,1)-20)<5 || mod_xy(n,1)<5 %ckecking all the dim seperate
%         if abs(mod_xy(n,2)-20)<5 || mod_xy(n,2)<5
            if abs(mod_z(n,3)-20)<3 || mod_z(n,3)<3
                filtered = cat(1,filtered,delta(n,:));
            end 
%         end
%     end
        
end
%3dplot of grid crosses
% figure;plot3(filtered(:,1),filtered(:,2),filtered(:,3),'o');
% axis equal

%make sure that all possible points have been found (all pillars to 13):
%histogram2(filtered(:,1),filtered(:,2),[15,11])
%edge points will not be detected because of the shape of the mask (might
%correct later)

% normal grid
slice_coord = readmatrix('slice-coord');
grid_norm = [];

for n = -6:6
    grid_norm = cat(1,grid_norm,[slice_coord n*20*ones(145,1)]);  
end
%
[Idx] = knnsearch(filtered,grid_norm);
sorted_ct = filtered(Idx,:);
sorted_ct = sorted_ct - sorted_ct(923,:); %adjust middle to 0,0,0

% registering using pointCloud objects to use CV toolbox
pc_ct = pointCloud(sorted_ct);
pc_norm = pointCloud(grid_norm);

[tform,ct_reg] = pcregistericp(pc_ct,pc_norm);

[~,D] = knnsearch(ct_reg.Location,grid_norm); %calculate distance of points
% signx = sign(sorted(:,1)-grid_norm(:,1));
% signy = sign(sorted(:,2)-grid_norm(:,2));
% signz = sign(sorted(:,3)-grid_norm(:,3));
% D_sign = D.*(signz);

% create 3d map
% D = transpose(D);

% D_slides = [];
% D_zeros = [];
% 
% for n = 1:13
%     D_slides = D((n-1)*145+1:n*145);
%     D_temp = [D_slides(1:15*7) 0 D_slides(15*7+1:15*7+13) 0 0 0 ...
%         D_slides(15*7+14:15*8+9) 0 0 0 0 0 D_slides(15*8+10:15*9+3) ...
%         zeros(1,7) D_slides(15*9+4:15*9+10) 0 0 0 0];
%     D_slice = [];
%     for row=1:11
%         D_slice=cat(1,D_slice,D_temp((row-1)*15+1:row*15));
%     end
%     D_zeros = cat(3,D_zeros,D_slice);
% end
% Vq = interp3(D_zeros,'cubic'); %use this to make plot smoother
% % Vq = D_zeros;
% figure;
% slice(Vq,1:size(Vq,2),1:size(Vq,1),1:size(Vq,3));
% %choose color and add transparency to slices
% colormap jet;
% alpha('color');
% alphamap('rampup');
% alphamap('decrease',0);
% colorbar

 
map = jet;

m = 63/(max(D));
b = 1;
%
figure; hold on; title('3D plot')
pcshow(pointCloud([0 0 0]));
colormap jet; colorbar;caxis([min(D) max(D)]);
colorbar('Ticks',[0,1,max(D)],...
         'TickLabels',{0,1,max(D)},'Color','w')
for point = 1:size(D)
    color = round(m*D(point)+b);
    plot3(grid_norm(point,1),grid_norm(point,3),-grid_norm(point,2),'+','Color',map(color,:),'MarkerSize',8,'LineWidth',2)
    
end


%%


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



%% Scatter Plot
figure;
scatter(D_centre(D~=-1),D(D~=-1),'x')
title('Scatter Plot of Distortion in CT');
xlabel('Distance from centre in mm');ylabel('Distortion in mm');
grid on;






    


