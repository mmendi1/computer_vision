
leftFile = "p11.jpg";
rightFile = "p12.jpg";
IL = imread(leftFile);
IR = imread(rightFile);

rectified_images = zeros(1200, 2000, 3, 'uint8');

F = -FundamentalRefinement(leftFile, rightFile); %ransac

disp("Matrix F for stereo pair 2:");
disp(F);
disp("Matrix rank:");
disp(rank(F));

e1 = null(F); % left
e2 = null(F'); % right

%[U,S,V] = svd(F);
%e1 = V(:,end); % left
%e2 = U(:,end); % right

e1(3) = 1;
e2(3) = 1;

H1 = [1 0 0; -e1(2)/e1(1) 1 0; -1/e1(1) 0 1]; % left homography
H2 = [1 0 0; -e2(2)/e2(1) 1 0; -1/e2(1) 0 1]; % right homography


for i=1:500 % y pixel
    for j=1:720 % x pixel
        homogeneous = [j i 1]';
        
        left_transformed =  H1 * homogeneous ; % H1
        right_transformed = H2 * homogeneous; % H2
        
        xl = round(left_transformed(1)); %x
        yl = round(left_transformed(2)); %y
        
        rectified_images(yl+400,xl+150,1) = IL(i,j,1);
        rectified_images(yl+400,xl+150,2) = IL(i,j,2);
        rectified_images(yl+400,xl+150,3) = IL(i,j,3);
        
        xr = round(right_transformed(1)); %x
        yr = round(right_transformed(2)); %y
        
        rectified_images(yr+350,xr+900,1) = IR(i,j,1); %250
        rectified_images(yr+350,xr+900,2) = IR(i,j,2); %250
        rectified_images(yr+350,xr+900,3) = IR(i,j,3); %250
        
    end
end

figure;
imshow(rectified_images);










