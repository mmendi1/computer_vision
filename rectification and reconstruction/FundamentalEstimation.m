function F = FundamentalEstimation()

% Pixel sampling
%imtool("p21.jpg");
%imtool("p22.jpg");

% Stereo pair 1:
%left_points = [[120 397]; [279 185]; [360 322]; [506 111]; [313 118]; [684 315]; [479 483]; [195 148]];
%right_points = [[21 378]; [230 183]; [310 315]; [548 97]; [468 107]; [720 314]; [202 478]; [352 142]];

% Stereo pair 2:
left_points = [[560 338];[546 456];[302 240];[478 68];[212 330];[660 112];[78 382];[140 86]];
right_points = [[512 338];[408 458];[288 240];[466 68];[170 330];[654 110];[50 378];[136 96]];


x1s = left_points(:,1);
y1s = left_points(:,2);

x2s = right_points(:,1);
y2s = right_points(:,2);

% Data normalization
x1_mean = mean(x1s);
y1_mean = mean(y1s);
x2_mean = mean(x2s);
y2_mean = mean(y2s);

variance1 = 1/16 * sum( (x1s-x1_mean).^2 + (y1s-y1_mean).^2);
variance2 = 1/16 * sum( (x2s-x2_mean).^2 + (y2s-y2_mean).^2);
std1 = sqrt(variance1);
std2 = sqrt(variance2);

M1 = [1/std1 0 -x1_mean/std1; 0 1/std1 -y1_mean/std1; 0 0 1];
M2 = [1/std2 0 -x2_mean/std2; 0 1/std2 -y2_mean/std2; 0 0 1];

x1_norm = (x1s - x1_mean)/std1;
y1_norm = (y1s - y1_mean)/std1;
x2_norm = (x2s - x2_mean)/std2;
y2_norm = (y2s - y2_mean)/std2;

% Setting up least squares problem
A = [x2_norm.*x1_norm x2_norm.*y1_norm x2_norm y2_norm.*x1_norm y2_norm.*y1_norm y2_norm x1_norm y1_norm ones(8,1)]; % least squares

% Least squares problem solution
[U1,S1,V1] = svd(A);
F_normalized = reshape(V1(:,end),3,3)';

% Denormalization
F_denormalized = (M2')*F_normalized*(M1); 

% Rank resolution
[U2,S2,V2] = svd(F_denormalized); 
S2(3,3) = 0; 
F = U2*S2*(V2');

% To make it match with the built in
F = -F/norm(F);

F_rank = rank(F);

fNorm8Point = estimateFundamentalMatrix(left_points,right_points,'Method','Norm8Point');

disp("Matrix F for stereo pair 2:");
disp(F);
disp("Matrix rank:");
disp(F_rank);
disp("It is similar to the one calculated by the built in function:");
disp(fNorm8Point);

end

