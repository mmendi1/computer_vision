
% Pair 1 
I1_double = esto;
I2_double = pair1rectifiedright;

% Pair 2
%I1_double = Pair2Left;
%I2_double = Pair2Right;

f1 = 31;
f2 = 32;
Tx = 100;
f = f2;

window = 11;
half_window = ceil(window/2);
shift = 40;

left_image = I1_double;
right_image = I2_double;

% going through the whole image
for y=1+window:500-window                                              
    for x=1+window:719-window                                           
        
        % creating boundaries for the window
        y_bound = y-half_window:y+half_window;
        x_bound = x-half_window:x+half_window;
        
        % storing left pixel values at window
        left_pixels=left_image(y_bound,x_bound);
        
        if (x+shift<719-half_window)
            ctr = 1;
            % shifting the window in the x direction in the right image
            for x1=x:x+shift
                % updating x boundary values for the window
                shifted_bound = x1-half_window:x1+half_window;
                % storing right pixel values at window
                right_pixels=right_image(y_bound,shifted_bound,1); 
                % computinng the sum of absolute difference between pixels
                % at window in the left and right image
                sums(ctr) = sum(abs(left_pixels-right_pixels),'all'); 
                ctr = ctr+1;
            end
        end
        % retrieving location with smallest absolute difference
        [~,i] = min(sums);  
        % creating a map of the disparities
        disparity_map(y,x) = i/shift;  
        % calculating depth from disparity
        depth = ceil((f*Tx)/disparity_map(y,x));
        % creating a map for the depth
        depth_map(y,x) = depth;
        
    end
end
                                    
depth_uint = zeros(493,712,1,'uint16');
for y=1:493
    for x=1:712
        depth_uint(y,x) = depth_map(y,x);
    end
end

% smoothing maps
disparity_map = imgaussfilt(disparity_map,2);
depth_uint = imgaussfilt(depth_uint,2);

figure;
imshow(disparity_map);  

figure;
imshow(depth_uint);