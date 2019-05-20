function [D] = traceback(info,D,D_iso,grid_norm,sorted_mr,map,slices_mr,center_phantom)
%TRACEBACK INPUT:  [dist_iso distortion] 
% OUTPUT: 3Dplot with marked point at MR coord
%   Detailed explanation goes here
m = 63/(max(D));
figure; hold on; title('3D plot with chosen point in pink')
xlabel('x');ylabel('y'),zlabel('z')
pcshow(pointCloud([0 0 0]));
colormap jet; colorbar;caxis([0 max(D)]);
colorbar('Ticks',[0,1,2,max(D)],...
         'TickLabels',{0,1,2,max(D)})
for point = 1:size(D)
    if D(point) >=0 %excluding points that haven't been found
    color = round(m*D(point)+1);
    plot3(grid_norm(point,1),grid_norm(point,3),-grid_norm(point,2),'+','Color',map(color,:),'MarkerSize',8,'LineWidth',2)
    end
    
end
distortion = info(2);
distance_iso = info(1);
number = find(D==distortion & D_iso==distance_iso);
coord = sorted_mr(number,:);
plot3(coord(1),coord(3),-coord(2),'Marker','*','MarkerSize',12,'Color',[1 0 1],'LineWidth',2)

%% printing images where point is on
coord = round(coord+center_phantom);

figure;title('Montage of Images around Point');
dim=size(slices_mr(:,:,coord(3)-4:coord(3)+4));
M=max(max(max(slices_mr)));
im_mount=reshape(slices_mr(:,:,coord(3)-4:coord(3)+4),[dim(1) dim(2) 1 dim(3)]);
montage(im_mount/M)

figure;hold on;title('Slice with Selected Point');
imshow(slices_mr(:,:,coord(3))/M);hold on
plot(coord(1),coord(2),'Marker','*','Color',[1 0 1],'MarkerSize',12)

if input('Do you want to delete this point? ','s') == 'yes'
    D(number) = -1;
    figure;
    scatter(D_iso(D~=-1),D(D~=-1),'x')
    title('Distortion ');
    xlabel('Distance from iso centre in mm');ylabel('Distortion in mm');
    grid on;
else
    D = D(:);
end


end

