clear all
close all
clc
%% Reading Images and Canny edge detection
im_left = imread('D:\Acads\8th sem\CV\cones-png-2\cones\im2.png');
figure(),   imshow(im_left);
title('Original Left Image')

left_qua_edge_map = im2double(edge(imresize(rgb2gray(im_left),0.25) ,'canny'));
figure(),imshow(left_qua_edge_map)
title('Edge Map of Left Image')

im_right = imread('D:\Acads\8th sem\CV\cones-png-2\cones\im6.png');
figure(),imshow(im_right);
title('Original Right Image')

right_qua_edge_map = im2double(edge(imresize(rgb2gray(im_right),0.25),'canny'));
figure(),imshow(right_qua_edge_map)
title('Edge Map of Right Image')

%% Finding Correaltion and hence disparity for all edges in edge map of left quarter size image
Corr = zeros(1,200);
index = zeros(size(left_qua_edge_map));

pad_leftqua = padarray(left_qua_edge_map,[2,2]);
pad_leftqua = padarray(pad_leftqua,[0,100],'post');
pad_rightqua = padarray(right_qua_edge_map,[2,2]);
pad_rightqua = padarray(pad_rightqua,[0,100],'post');

for x= 3: size(left_qua_edge_map,1)+2
    for y = 3:size(left_qua_edge_map,2)+2
        for d = 1:100;
            for i= -2:1:2
                for j = -2:1:2
                   if(pad_leftqua(x,y)==1)
                       Corr(d) = Corr(d) + pad_leftqua(x+i,y+j)*pad_rightqua(x+i,y+d+j);
                   end
                end
            end
        end
        
[value, index(x-2,y-2)] = max(Corr); % index storing disparities where correlation is max ;"sparse disparity map" 
Corr = zeros(1,200);
   end
end

figure(),imshow(uint8(index),[0 255])
title('Sparse Disparity Map');
index = index- ones(size(index)); % for making disparities of non edge points as 0 for convenience

%% Interpolation from sparse map to dense map
disp_pts = [];  % storing disparities at edged points
edge_pts = [];  % storing coloumn of edges , row wise

for x= 1: size(left_qua_edge_map,1)
    k = x ;l = 1 ;
    
    edge_pts{k,1} = [];
    disp_pts{k,1} = [];
    
    for y = 1:size(left_qua_edge_map,2)
       if(left_qua_edge_map(x,y)==1)
           disp_pts{k,1}(l) = index(x,y);
           edge_pts{k,1}(l) = y;
           l = l + 1;
       end
   end
end

l = 1; count = 1;

dense_map = zeros(size(index));
for k= 2:(size(index,1)-1)
    if (abs(disp_pts{k,1}(1) - disp_pts{k,1}(2)) > 10)
           dense_map(k,1:edge_pts{k,1}(1)) = disp_pts{k,1}(1); 
            l = 2;
    elseif (abs(disp_pts{k,1}(1) - disp_pts{k,1}(2)) <= 10)
           x  = [edge_pts{k,1}(1),edge_pts{k,1}(2)];
           v  = [index(k,edge_pts{k,1}(1)),index(k,edge_pts{k,1}(2))];
           xq = edge_pts{k,1}(1):edge_pts{k,1}(2);
           
           vq = interp1(x,v,xq);     %interpolating between edge pairs
           dense_map(k,edge_pts{k,1}(1):edge_pts{k,1}(2)) = vq ; 
           l = 3 ;
    end
        
     count = 0;
    while (l < size(disp_pts{k,1},2))      
       if(abs(disp_pts{k,1}(l) - disp_pts{k,1}(l+1)) <= 10)
           x = [edge_pts{k,1}(l),edge_pts{k,1}(l+1)];
           v = [index(k,edge_pts{k,1}(l)),index(k,edge_pts{k,1}(l+1))];
           xq= edge_pts{k,1}(l):edge_pts{k,1}(l+1);
           
           vq= interp1(x,v,xq); %interpolating between edge pairs
           dense_map(k,edge_pts{k,1}(l):edge_pts{k,1}(l+1)) = vq ; 
           l = l+2 ;
        
        elseif (abs(disp_pts{k,1}(l) - disp_pts{k,1}(l+1)) > 10)
           dense_map( k , edge_pts{k,1}(l) : edge_pts{k,1}(l+1)) = disp_pts{k,1}(l) ;
           l = l+1 ;
       end
    count   = count +1;
    end
end

dense_map1 = 255*(dense_map./max(max(dense_map)));
figure(),imshow(uint8(dense_map1),[0 255])
title('Disparity Dense Map')

depth_dense_map = 1./dense_map;
figure(),imshow(uint8(depth_dense_map),[0 255])
title('Depth Dense Map')