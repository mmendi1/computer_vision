
ILeft1 = "p11.jpg";
IRight1 = "p12.jpg";

ILeft2 = "p21.jpg";
IRight2 = "p22.jpg";

K1 = [31/(31.2e-6) 0 720/2; 0 31/(31.2e-6) 500/2; 0 0 1]; % stereo pair 1
K2 = [32/(31.2e-6) 0 720/2; 0 32/(31.2e-6) 500/2; 0 0 1]; % stereo pair 2


F1 = FundamentalRefinement(ILeft1,IRight1);
F2 = FundamentalRefinement(ILeft2,IRight2);

% getting E from F
E1 = K1' * F1 * K1;
E2 = K2' * F2 * K2;

% finding the extrinsics from E
[U1,S1,V1] = svd(E1);
[U2,S2,V2] = svd(E1);

W = [0 -1 0; 1 0 0; 0 0 1];

R11 = U1*W*V1';
R12 = U1*W'*V1';
R21 = U2*W*V2';
R22 = U2*W'*V2';

t11 = U1(:,end);
t12 = -U1(:,end);

t21 = U2(:,end);
t22 = -U2(:,end);

R1 = R12;
T1 = t11;

R2 = R21;
T2 = t21;


% 5 depth levels - computer - table front - wall(switch) - risk - bonaparte
left_points1 = [[698 108] [482 454] [316 120] [472 232] [122 398]];  
right_points1 = [[738 82] [202 448] [470 112] [420 226] [22 380]]; 
% expected ordering from smallest to largest 
% depth2
% depth5
% detph4
% depth1 
% depth3

% Calculating depths for image pair 1
pl1 = [698 108 31]';
pr1 = [738 82 31]';
depth1 = get_depth(pl1,pr1,R1,T1);

pl2 = [482 454 31]';
pr2 = [202 448 31]';
depth2 = get_depth(pl2,pr2,R1,T1);

pl3 = [316 120 31]';
pr3 = [470 112 31]';
depth3 = get_depth(pl3,pr3,R1,T1);

pl4 = [472 232 31]';
pr4 = [420 226 31]';
depth4 = get_depth(pl4,pr4,R1,T1);

pl5 = [122 398 31]';
pr5 = [22 380 31]';
depth5 = get_depth(pl5,pr5,R1,T1);

% 5 depth levels - ball - telephone - head - poster on wall - mug
left_points2 = [[560 338] [546 456] [302 240] [478 68] [212 330]];
right_points2 = [[512 338] [408 458] [288 240] [466 68] [170 330]];
% expected ordering from smallest to largest 
% depth2
% depth5
% detph1
% depth3 
% depth4

% Calculating depths for image pair 2
pl12 = [560 338 32]';
pr12 = [512 338 32]';
depths21 = get_depth(pl12,pr12,R2,T2);

pl22 = [546 456 32]';
pr22 = [408 458 32]';
depths22 = get_depth(pl22,pr22,R2,T2);

pl32 = [302 240 32]';
pr32 = [288 240 32]';
depths23 = get_depth(pl32,pr32,R2,T2);

pl42 = [478 68 32]';
pr42 = [466 68 32]';
depths24 = get_depth(pl42,pr42,R2,T2);

pl52 = [212 330 32]';
pr52 = [170 330 32]';
depths25 = get_depth(pl52,pr52,R2,T2);


function depth = get_depth(pl, pr, R2, T)
    pleft = pl;
    pright = R2' * pr;
    w = cross(pleft,pright); % cross producted
    matrix = [pleft pright w]; % solving system of linear equations
    unknown = matrix\T; %inv(matrix)
    a = unknown(1);
    b = -unknown(2);
    endpoint1 = a*pleft; % origin of left camera is (0,0,0)
    endpoint2 = T + b*pright;
    z_midpoint = (endpoint1(3)+endpoint2(3))/2;
    depth = z_midpoint;
end