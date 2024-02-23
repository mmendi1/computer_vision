% Smoothing frames
% A line of code was added in demo_optical_flow.m so that all frames are
% smoothed before using them. The smoothing factor used is 0.5

function [Vx,Vy] = compute_LK_optical_flow(frame_1,frame_2)

im1 = single(rgb2gray(frame_1));
im2 = single(rgb2gray(frame_2));

% Initializing Parameters
max_level = 6; % total number of levels in pyramids
window_size = 13; 

% Defining kernels
kernel_x = [-1/2 1/2; -1/2 1/2];
kernel_y = [-1/2 -1/2; 1/2 1/2];
kernel_t = [1/2 1/2; 1/2 1/2];

% Creating pyramids for each image
initial_size1 = size(im1,1);
initial_size2 = size(im1,2);

% Pyramid for im1
pyramid_im1(1:initial_size1,1:initial_size2,1) = im1;
pyramid_im1(1:initial_size1/2,1:initial_size2/2,2) = impyramid(pyramid_im1(1:initial_size1,1:initial_size2,1),'reduce');
pyramid_im1(1:initial_size1/4,1:initial_size2/4,3) = impyramid(pyramid_im1(1:initial_size1/2,1:initial_size2/2,2), 'reduce');
pyramid_im1(1:initial_size1/8,1:initial_size2/8,4) = impyramid(pyramid_im1(1:initial_size1/4,1:initial_size2/4,3), 'reduce');
pyramid_im1(1:initial_size1/16,1:initial_size2/16,5) = impyramid(pyramid_im1(1:initial_size1/8,1:initial_size2/8,4), 'reduce');
pyramid_im1(1:initial_size1/32,1:initial_size2/32,6) = impyramid(pyramid_im1(1:initial_size1/16,1:initial_size2/16,5), 'reduce');

% Pyramid for im2
pyramid_im2(1:initial_size1,1:initial_size2,1) = im2;
pyramid_im2(1:initial_size1/2,1:initial_size2/2,2) = impyramid(pyramid_im2(1:initial_size1,1:initial_size2,1), 'reduce');
pyramid_im2(1:initial_size1/4,1:initial_size2/4,3) = impyramid(pyramid_im2(1:initial_size1/2,1:initial_size2/2,2), 'reduce');
pyramid_im2(1:initial_size1/8,1:initial_size2/8,4) = impyramid(pyramid_im2(1:initial_size1/4,1:initial_size2/4,3), 'reduce');
pyramid_im2(1:initial_size1/16,1:initial_size2/16,5) = impyramid(pyramid_im2(1:initial_size1/8,1:initial_size2/8,4), 'reduce');
pyramid_im2(1:initial_size1/32,1:initial_size2/32,6) = impyramid(pyramid_im2(1:initial_size1/16,1:initial_size2/16,5), 'reduce');


% Iterating from max_level to first and saving corresponding dimensions per level
level_info = {};
for level = max_level:-1:1
    dimension1 = round(initial_size1/2^(level-1));
    dimension2 = round(initial_size2/2^(level-1));
    level_info{end+1} = [dimension1 dimension2 level];
end

% Saving corresponding image at each level
images_p1 = {};
images_p2 = {};
for ctr = 1:max_level
    current_p1 = pyramid_im1(1:level_info{ctr}(1), 1:level_info{ctr}(2), level_info{ctr}(3));
    current_p2 = pyramid_im2(1:level_info{ctr}(1), 1:level_info{ctr}(2), level_info{ctr}(3));
    images_p1{end+1} = current_p1;
    images_p2{end+1} = current_p2;
end

% At largest scale, we initialize the vectors to be zero
final_size1 = initial_size1/32;
final_size2 = initial_size2/32;
hx=zeros(final_size1,final_size2);
hy=zeros(final_size1,final_size2);

for ctr = 1:max_level-1 % looping through the levels of the pyramids
    current_im1 = images_p1{ctr};
    current_im2 = images_p2{ctr};
    current_size = size(current_im1);
     
    % Interpolating estimates from next largest scale
    hx = imresize(hx,[level_info{ctr}(1) level_info{ctr}(2)]);
    hy = imresize(hy,[level_info{ctr}(1) level_info{ctr}(2)]);
 
    % Sliding through the whole image making sure window will fit (corners
    % are cropped so that coordinates are not out of bounds)
     for i = 1+window_size:current_size(1)-window_size
         for j = 1+window_size:current_size(2)-window_size
             
         % get window boundaries
            x = j-floor(window_size/2):j+floor(window_size/2);
            y = i-floor(window_size/2):i+floor(window_size/2);
  
         % access image and shifted image at window
            image = current_im1(y,x);
            shifted_image = current_im2(y + round(hy(i,j)), x + round(hx(i,j)));
            
         % GRADIENTS
%             [Gx_image,Gy_image] = imgradientxy(image);
%             [Gx_shifted, Gy_shifted] = imgradientxy(shifted_image);
%             Gx = Gx_image + Gx_shifted;
%             Gy = Gy_image + Gy_shifted
            
            Ix_image = conv2(image,kernel_x);
            Iy_image = conv2(image,kernel_y);
            Ix_shifted = conv2(shifted_image, kernel_x);
            Iy_shifted = conv2(shifted_image, kernel_y);
            
            ins = 2:window_size-2;
            
            Ix = reshape((Ix_image(ins,ins) + Ix_shifted(ins,ins))',[],1);
            Iy = reshape((Iy_image(ins,ins) + Iy_shifted(ins,ins))',[],1);
            
            A = [Ix Iy];
            
            lhs = A'*A; %AtA left hand side of equation
            
            It_image = conv2(image,kernel_t);
            It_shifted = conv2(shifted_image, kernel_t);
            It = reshape((It_shifted(ins,ins) - It_image(ins,ins))',[],1);
            
            rhs = A'*It; %Atb right hand side of equation
            
            new_shift = lhs\rhs; %lhs inverse times rhs
            hx(i,j)=hx(i,j)+new_shift(1); % update x shift
            hy(i,j)=hy(i,j)+new_shift(2); % update y shift
        end
    end
end
    % return final estimates
    Vx = hx;
    Vy = hy; 
end




