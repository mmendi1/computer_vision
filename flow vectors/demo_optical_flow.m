% Only changes made were to make the code compile
% A line of code was also added to smooth all frames prior to using them

function [Vx,Vy] = demo_optical_flow(folder_name,frame_number_1,frame_number_2)

if(nargin == 0)
    folder_name = 'Backyard';
    frame_number_1 = 7;
    frame_number_2 = frame_number_1 + 1;
elseif(nargin == 1)
    frame_number_1 = 7;
    frame_number_2 = frame_number_1 + 1;
elseif(nargin ==2)
    frame_number_2 = frame_number_1 + 1;
end

addpath(folder_name);

% This line was added to smooth all frames before calling compute_LK_optical_flow
smooth_frames(folder_name,0.5,7,14);

frame_1 = read_image(folder_name,frame_number_1);
frame_2 = read_image(folder_name,frame_number_2);

[Vx,Vy] = compute_LK_optical_flow(frame_1,frame_2);

plotflow(Vx,Vy);
title('Quiver plot');

end

function plotflow(Vx,Vy)

s = size(Vx);
step = max(s)/80; %max(s)/100
[X, Y] = meshgrid(1:step:s(2), s(1):-step:1);
u = interp2(Vx, X, Y);
v = interp2(Vy, X, Y);

%figure;
%hold on
%quiver(X,Y,u,v,1,'cyan','LineWidth',1);
quiver(X, -Y, u, -v, 1, 'cyan', 'LineWidth', 1);

axis image;

end

function I = read_image(folder_name,index)

if(index < 10)
    I = imread(fullfile(folder_name,strcat('frame0',num2str(index),'.png')));
else
    I = imread(fullfile(folder_name,strcat('frame',num2str(index),'.png')));
end

end
