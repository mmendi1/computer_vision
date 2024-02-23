%% Comp 558: Assignment 1
%% Question 2
%% Author: Malena Mendilaharzu
 
%% Question 2 Part a
 
sigma = 1; 
y_sigma = 2;
 
syms x
formula = normpdf(x,0,sigma);
second_diff = diff(formula,x,2);
 
% Computing second derivatives in the x direction
derivatives = {};
for coord=-2*sigma:2*sigma
    derivatives{end+1} = subs(second_diff,{x},coord);
end
 
% Populating the filter with the second derivatives
filter = zeros((4*sigma)+1,(4*sigma)+1);
for i=1:(4*sigma+1)
    filter(:,i) = derivatives{i};
end

% Images
%figure;					
%imagesc(filter)
%colormap jet;
%axis image
%axis off;
%result1 = filter;
%figure;
%surfc(result1);
%colorbar

%% Question 2 Part b
 
% Multiplying entries by a normalized Gaussian in the y direction
current = 2*sigma;
for i=1:(4*sigma+1)
    a = filter(1,:);
    filter(i,:) = filter(i,:)*normpdf(current,0,y_sigma);
    current = current - 1;
end
 
% Images
%figure;					
%imagesc(filter)
%colormap jet;
%axis image
%axis off;
%result2 = filter;
%figure;
%surfc(result2)
%colorbar
 
 
%% Question 2 Part c
 
% Creating a bigger array so that we don't lose entries when performing the
% rotation
z0 = zeros((4*sigma)+1,(4*sigma)+1);
all = [z0 z0 z0
    z0 filter z0 
    z0 z0 z0];
 
% Following matrix simulates coordinates with (0,0) located in the center of the
% array
for x=1:(12*sigma)+3
    for y=1:(12*sigma)+3
        coordinates(y,x)= mat2cell(zeros(2,1),2,1);
        coordinates{y,x} = [-6*sigma+(x-2),6*sigma-(y-2)]; 
    end
end

% Following matrix rotates the coordinates of the previous matrix by theta
% Notice that the procedure works for all theta
% Only pi/4 is provided for demonstration
theta = pi/4;
for x=1:(12*sigma)+3
    for y=1:(12*sigma)+3
        rot_x = ((coordinates{y,x}(1)*cos(theta) - coordinates{y,x}(2)*sin(theta)));
        rot_y = ((coordinates{y,x}(1)*sin(theta) + coordinates{y,x}(2)*cos(theta)));
        rotated_coordinates{y,x} = [round(rot_x),round(rot_y)];
    end
end
 
% Rotating the filter
% Please refer to function in the end 

rotation0 = all;
rotation45 = rotate_matrix(all,coordinates, rotated_coordinates);
rotation90 = all';
rotation135 = rotate_matrix(rotation90,coordinates, rotated_coordinates);

% Rotation 0 
%figure;					
%imagesc(rotation0)
%colormap jet;
%axis image
%axis off;
%figure;
%surfc(rotation0)
%colorbar

% Rotation45
%rotation45 = smooth_matrix(rotation45);
%figure;					
%imagesc(rotation45)
%colormap jet;
%axis image
%axis off;
%figure;
%surfc(rotation45)
%colorbar

% Rotation 90
%figure;					
%imagesc(rotation90)
%colormap jet;
%axis image
%axis off;
%figure;
%surfc(rotation90)
%colorbar

% Rotation 135
%figure;					
%imagesc(rotation135)
%colormap jet;
%axis image
%axis off;
%rotation135 = smooth_matrix(rotation135);
%figure;
%surfc(rotation135)
%colorbar
 
 
%% Question 2 Part d
 
original = imread('skyscrapers.jpg');
grayscale = rgb2gray(original);

% Convolving an image with 
 
I0 = conv2(rotation0,grayscale);
I45 = conv2(rotation45,grayscale);
I90 = conv2(rotation90,grayscale);
I135 = conv2(rotation135,grayscale);

%figure;
%imshow(I0);
%figure;
%imshow(I45);
%figure;
%imshow(I90);
%figure;
%imshow(I135);
 
% Finding zero crossings
% For each: % x-direction, y-direction, diagonal
% Then fusing then ploting

figure;
%imshow(I0);
zc0 = zeros(size(I0));
zc0(:,:,1) = (sign(circshift(I0,-1,1))~=sign(I0));
zc0(:,:,2) = (sign(circshift(I0,-1,2))~=sign(I0));
zc0(:,:,3) = sign(circshift(I0,[-1,-1]))~=sign(I0);
PT1 = imfuse(zc0(:,:,1),zc0(:,:,2));
PT2 = imfuse(PT1,zc0(:,:,3));
%imshow(PT2);

figure;
%imshow(I45);
zc45 = zeros(size(I45));
zc45(:,:,1) = (sign(circshift(I45,-1,1))~=sign(I45));
zc45(:,:,2) = (sign(circshift(I45,-1,2))~=sign(I45));
zc45(:,:,3) = sign(circshift(I45,[-1,-1]))~=sign(I45);
PT1 = imfuse(zc45(:,:,1),zc45(:,:,2));
PT2 = imfuse(PT1,zc45(:,:,3));
%imshow(PT2);

figure;
%imshow(I90);
zc90 = zeros(size(I90));
zc90(:,:,1) = (sign(circshift(I90,-1,1))~=sign(I90));
zc90(:,:,2) = (sign(circshift(I90,-1,2))~=sign(I90));
zc90(:,:,3) = sign(circshift(I90,[-1,-1]))~=sign(I90);
PT1 = imfuse(zc90(:,:,1),zc90(:,:,2));
PT2 = imfuse(PT1,zc90(:,:,3));
%imshow(PT2);

figure;
%imshow(I135);
zc135 = zeros(size(I135));
zc135(:,:,1) = (sign(circshift(I135,-1,1))~=sign(I135));
zc135(:,:,2) = (sign(circshift(I135,-1,2))~=sign(I135));
zc135(:,:,3) = sign(circshift(I135,[-1,-1]))~=sign(I135);
PT1 = imfuse(zc135(:,:,1),zc135(:,:,2));
PT2 = imfuse(PT1,zc135(:,:,3));
%imshow(PT2);


%% Question 2 Part e


LoG = fspecial('log',40,6);
           
ILoG = conv2(grayscale,LoG, 'same');
temp = zeros(size(grayscale));
%circshift in x-direction   
temp(:,:,1) = (sign(circshift(ILoG,-1,1))~=sign(ILoG));
%circshit in y-direction
temp(:,:,2) = (sign(circshift(ILoG,-1,2))~=sign(ILoG));
%circshift in diagonal-direction
temp(:,:,3) = (sign(circshift(ILoG,[-1,-1]))~=sign(ILoG));
   

figure('Name','Laplacian of Gaussian - Zero crossings');
%subplot(2,2,1),imshow(temp(:,:,1));
%subplot(2,2,2),imshow(temp(:,:,2));
%subplot(2,2,3),imshow(temp(:,:,3));

PT1 = imfuse(temp(:,:,1),temp(:,:,2));
PT2 = imfuse(PT1,temp(:,:,3));
%imshow(PT2);

%LoG=fspecial('log',40,6);
%LoG90=LoG';
%figure;
%surfc(LoG);
%colorbar
%ILoG = conv2(grayscale,LoG,'same');
%temp(:,:,1) = sign(circshift(ILoG,[-1,-1]))~=sign(ILoG);

 
%% Functions
 
% This function rotates a matrix by using the coordinates and the
% rotated_coordinates matrices.
% It takes the original coordinates, the finds the corresponding rotated
% coordinates and finally updates the values of the matrix.
% It "reorganizes" the matrices.

function rot = rotate_matrix(matrix, coords, rotcoords)
    [a,b] = size(matrix);
    rot = zeros(a,b);
    for y=1:b
        for x=1:a
            current = matrix(x,y);
            element = rotcoords{x,y};
            indexx=1;
            indexy=1;
            for c1=1:a
               for c2=1:a
                   if(element == cell2mat(coords(c2,c1)))
                       indexx = c2;
                       indexy = c1;
                   end
               end
            end
            rot(indexx,indexy) = current;
        end
    end
end



% This function was used to smooth the filters
% It performs a variation of a 2D interpolation to get rid of undesired
% peaks.
% Performance improved in the appearance of the filter but not that much in
% the results. 

% function smoothed = smooth_matrix(matrix)
% [~,c] = size(matrix);
% smoothed = {};
% %smoothed{1} = (matrix(1,:)+matrix(2,:))/2;
% smoothed{1} = (matrix(1,:)+matrix(2,:))/2;
% %smoothed{1} = (matrix(:,1)+matrix(:,2))/2;
% for i=2:c-1
%     %smoothed{i} = (matrix(i-1,:)+matrix(i,:)+matrix(i+1,:))/3;
%     smoothed{i} = (matrix(i,:)+matrix(i+1,:))/2;
%     %smoothed{i} = (matrix(:,i)+matrix(:,i+1))/2;
% end
% %smoothed{c} = (matrix(c-1,:)+matrix(c,:))/2;
% smoothed{c} = (matrix(c,:)+matrix(c-1,:))/2;
% %smoothed{c} = (matrix(:,c)+matrix(:,c-1))/2;
% smoothed = cell2mat(smoothed);
% smoothed = reshape(smoothed,c,c);
% end
