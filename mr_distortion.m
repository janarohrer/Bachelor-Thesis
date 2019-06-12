close all
clear all

%% load folder and images
folder = 'thesis_data/SHIFT_X_40MM_RETAKE';
[image_mr,info_mr] = dicomfolder(folder);
%%
dimx = round(info_mr{1}.PixelSpacing(1)*size(image_mr,1));
dimy = round(info_mr{1}.PixelSpacing(2)*size(image_mr,2));
dimz = round(info_mr{1}.SliceThickness*size(image_mr,3));
iso = transpose(-round(info_mr{1}.ImagePositionPatient+[0;0;47]));
%resize voxels to be 1mm x 1mm x 1mm
slices_mr = imresize3(image_mr,[dimx,dimy,dimz]);
slices_mr = slices_mr(1:284,:,48:312); 

%% denoise and gradient mr (VERY time intensive)
mr_den = zeros(size(slices_mr));
for z_coord = 1:size(slices_mr,3)
    M=max(max(slices_mr(:,:,z_coord)));  
    mr_den(:,:,z_coord) = wdenoise2(slices_mr(:,:,z_coord)/M,8,'CycleSpinning',3,'DenoisingMethod','UniversalThreshold');
end
mr_grad = imgradient3(mr_den);
%% using imdiffusefilt
mr_den = imdiffusefilt(slices_mr,'NumberOfIterations',3);
mr_den = mr_den(:,:,:)/max(max(max(mr_den)));
mr_imd = imgradient3(mr_den);
mr_grad = mr_imd/max(max(max(mr_imd)));
%% make 3d template of grid cross
dim = 23; 
middle = round(dim/2);
mask_3 = zeros(dim,dim,dim); %make 3d template of grid cross
mask_3(:,middle-1:middle+1,middle) = 1;
mask_3(:,middle,middle-1:middle+1) = 1;
mask_3(middle,:,middle-1:middle+1) = 1;
mask_3(middle-1:middle+1,:,middle) = 1;
mask_3(middle,middle-1:middle+1,:) = 1;
mask_3(middle-1:middle+1,middle,:) = 1;

mask_3(middle,middle,:) = 0;
mask_3(middle,:,middle) = 0;
mask_3(:,middle,middle) = 0;

mask_left = mask_3(:,11:end,:);
mask_right = mask_3(:,1:13,:);


%% crosscorr in 3 dimensions to find crosses (time intensive)
cross_corr_3 = real(normxcorr3(mask_3,mr_grad,'same')); %use mr_wd for bad contrast or big shifts, mr_imd else
cross_corr_left = real(normxcorr3(mask_left,mr_grad(:,1:100,:),'same'));
cross_corr_right = real(normxcorr3(mask_right,mr_grad(:,end-99:end,:),'same'));


cross_corr_3 = cross_corr_3 / max(max(max(cross_corr_3)));
cross_corr_left = cross_corr_left(:,8:52,:)/max(max(max(cross_corr_left(:,8:52,:))));
cross_corr_right = cross_corr_right(:,end-51:end-7,:)/max(max(max(cross_corr_right(:,end-51:end-7,:))));

cross_corr_3(:,1:45,:) = (cross_corr_3(:,1:45,:)+cross_corr_left)/max(max(max(cross_corr_3(:,1:45,:)+cross_corr_left)));
cross_corr_3(:,end-44:end,:) = (cross_corr_3(:,end-44:end,:)+cross_corr_right)/max(max(max(cross_corr_3(:,end-44:end,:)+cross_corr_right)));


%% smooth and filter
cross_corr_4 = medfilt3(cross_corr_3);
cross_corr_4 = cross_corr_4/max(max(max(cross_corr_4)));
filter_spot=zeros(size(cross_corr_4));
filter_spot(cross_corr_4>0.5)=1;
   
% find centroids
CC = bwconncomp(filter_spot);
S = regionprops(CC,'Centroid','Area');
s_area = [S.Area];
centroids = cat(1, S.Centroid); 

%% manually find center of phantom
central = round(size(slices_mr,3)/2);
x_shift = input('enter shift in x in mm: ');
z_shift = input('enter shift in z in mm: ');
central_shift = central + z_shift; %z coord of center

%display images around center to find central slice
pad = 4 + min([mod(z_shift,20) mod(-z_shift,20)]);
figure;
dim=size(slices_mr(:,:,central_shift-pad:central_shift+pad));
M=max(max(max(slices_mr)));
im_mount=reshape(slices_mr(:,:,central_shift-pad:central_shift+pad),[dim(1) dim(2) 1 dim(3)]);
montage(im_mount/M)
central_z = input('enter the index of the middle slice (most prominent grid): ')- pad+central_shift;

circ = zeros(8);
% choose control points to determine middle in xy plane
%choose centre and save fixed point to fixedPoints:
[movingPoints,fixedPoints] = cpselect(circ,slices_mr(:,:,central_z)/max(max(slices_mr(:,:,central_z))),[1 1],iso(1:2),'Wait',true);

%% mod filter coordinates 
%only taking points that are close to expected coordinates
center_phantom = [fixedPoints(1)+x_shift fixedPoints(2) central_z-z_shift]; %centre of phantom defined by user input
delta = (centroids-center_phantom);

mod_xyz = mod(delta,20); 
filtered_mr = [];
area = [];
for n = 1:size(centroids,1)
    if (abs(mod_xyz(n,1)-20)<5 || mod_xyz(n,1)<5) && abs(delta(n,1))< 7*20+10 %ckecking all the dim seperate
        if (abs(mod_xyz(n,2)-20)<5 || mod_xyz(n,2)<5) && delta(n,2)> -(7*20+10) && delta(n,2)<3*20+10
            if abs(mod_xyz(n,3)-20)<3 || mod_xyz(n,3)<3 
                filtered_mr = cat(1,filtered_mr,delta(n,:));
                area = cat(1,area,s_area(n));    
            end 
        end
    end     
end

load grid_norm

[Idx] = knnsearch(filtered_mr,grid_norm);
sorted_mr = filtered_mr(Idx,:);
sorted_mr = sorted_mr - sorted_mr(923,:); %adjust middle to 0,0,0
sorted_area = area(Idx);

% registering using pointCloud objects to use CV toolbox
pc_mr = pcdenoise(pointCloud(sorted_mr));
pc_norm = pointCloud(grid_norm);

[tform,mr_reg] = pcregistericp(pc_mr,pc_norm);
[~,D] = knnsearch(mr_reg.Location,grid_norm); %calculate distance of points

%points that are further away than 10 should not be taken into account -->
%exclude points that have not been found in the mr
D_red=D;
D_red(D>10) = -1; 
D_red(sorted_area<60 & D>2.5) = -1; %exclude small points with big distortion
D = D_red;
%% plots
map = jet;

m = 63/(max(D));
b = 1;
maxd = max(D);
%
figure; hold on; title('Distortion in mm')
pcshow(pointCloud([0 0 0]));
colormap jet; colorbar;caxis([0 maxd]);
cb = colorbar('Ticks',[0,1,2,maxd],...
         'TickLabels',{0,1,2,maxd},'Color','w');
for point = 1:size(D)
    if D(point) >=0 %excluding points that haven't been found
    color = round(m*D(point)+1);
    plot3(grid_norm(point,1),grid_norm(point,3),-grid_norm(point,2),'+','Color',map(color,:),'MarkerSize',8,'LineWidth',2)
    end
    
end


% dist. from iso center vs distortion
%iso-center_phantom: coordinates of isocenter in rel. to center of phantom:
[~,D_iso] = knnsearch(iso-center_phantom,grid_norm); 
D_5=zeros(size(D_iso));
D_5(D_iso<=50)=D(D_iso<=50);
D_10=zeros(size(D_iso));
D_10(D_iso<=100)=D(D_iso<=100);
D_20=zeros(size(D_iso));
D_20(D_iso<=200)=D(D_iso<=200);
'maximal distortion in r = 5cm: ', max(D_5)
'maximal distortion in r = 10cm: ', max(D_10)
'maximal distortion in r = 20cm: ', max(D_20)
'points that were not found: ', size(find(D==-1),1)
'points not found within r=10cm: ', size(find(D_10==-1),1)

% Scatter Plot
figure;
scatter(D_iso(D~=-1),D(D~=-1),'x')
title('Control Rotation 0°');
xlabel('Distance from iso centre in mm');ylabel('Distortion in mm');
grid on;

%% Traceback
%click on scatter plot point and save it to cursor_info --> reverse
%pipeline, it also offers the option to delete the point

D = traceback(cursor_info.Position,D,D_iso,grid_norm,sorted_mr,map,slices_mr,center_phantom);

%% final plots and statistics

% m = 63/(max(D)); %uncomment to adjust to new maximum
figure; hold on; title('Distortion in mm')
pcshow(pointCloud([0 0 0]));
colormap jet; colorbar;caxis([0 maxd]);
cb = colorbar('Ticks',[0,1,2,maxd],...
         'TickLabels',{0,1,2,maxd},'Color','w');
     
for point = 1:size(D)
    if D(point) >=0 %excluding points that haven't been found
    color = round(m*D(point)+1);
    plot3(grid_norm(point,1),grid_norm(point,3),-grid_norm(point,2),'+','Color',map(color,:),'MarkerSize',8,'LineWidth',2)
    end
    
end

D_5=zeros(size(D_iso));
D_5(D_iso<=50)=D(D_iso<=50);
D_10=zeros(size(D_iso));
D_10(D_iso<=100)=D(D_iso<=100);
D_20=zeros(size(D_iso));
D_20(D_iso<=200)=D(D_iso<=200);
'maximal distortion in r = 5cm: ', max(D_5)
'maximal distortion in r = 10cm: ', max(D_10)
'maximal distortion in r = 20cm: ', max(D_20)
'points that were not found: ', size(find(D==-1),1)
'points not found within r=10cm: ', size(find(D_10==-1),1)
% Scatter Plot
figure;
scatter(D_iso(D~=-1),D(D~=-1),'x')
title('Rotation 90°');
xlabel('Distance from iso centre in mm');ylabel('Distortion in mm');
grid on;