%% Comp 558 : Assignment 1
%% Question 3
%% Author: Malena Mendilaharzu


%% Script

im = imread('skyscrapers.jpg');
im = imresize(rgb2gray(im), 0.5);
Iedges = edge(im,'canny');
[~,grad_dir]=imgradient(im);
grad_dir = - grad_dir;
%imshow(Iedges)
[row, col] = find(Iedges);
edges = [col, row, zeros(length(row),1), zeros(length(row),1) ];
for k = 1:length(row)
     edges(k,3) = cos(grad_dir(row(k),col(k))/180.0*pi);
     edges(k,4) = sin(grad_dir(row(k),col(k))/180.0*pi);
end


%% Question 3: Part a

M = 6000;
threshold_distance = 5;
threshold_angle= 0.1;
c_min= 70;
line_models = {};


for i=1:M
    disp(i)
    rs = datasample(edges,1);
    m = rs(4)/rs(3); 
    rs_line = m*(0-rs(1))+rs(2);
    c = rs_line - m*rs(1);
    consensus_set = [];
    
    for j=1:length(edges)
        if (edges(j,:)==rs)
            continue
        end
        
        x1 = 0;
        y1 = rs_line;
        x2 = rs(1);
        y2 = rs(2);
        x3 = edges(j,1);
        y3 = edges(j,2);
        
        dist_var = calc_distance(x1,y1,x2,y2,x3,y3); % this is the point line distance
        theta_var = abs(acos(rs(3))-acos(edges(j,3)));
        %theta_var = abs(atan2(rs(4),rs(3)) - atan2(edges(j,4),edges(j,3)));
        
        if((dist_var<=threshold_distance) && (theta_var<=threshold_angle))
            consensus_set = [consensus_set; edges(j,1) edges(j,2) acos(edges(j,3)) j];
        end
    end
    
    if (length(consensus_set) > c_min)
        disp("adentro");
        average_x = mean(consensus_set(:,1));
        average_y = mean(consensus_set(:,2));  
        average_angle = mean(consensus_set(:,3));
        line_model = [average_x average_y cos(average_angle) j];
        line_models{end+1} = line_model;
        edges = setdiff(edges,consensus_set,'rows');
        edges = [edges; line_model];
    end
    
end

upperbound = length(line_models);
lowerbound = upperbound - 10 + 1;

figure;
imshow(im)
hold on
for k=lowerbound:upperbound
    x0 = line_models{k}(1);
    y0 = line_models{k}(2);
    %slope = sin(line_models{k}(3))/cos(line_models{k}(3));
    angle = acos(line_models{k}(3));
    x = 0 : 250;
    y = tan(angle)*(x-x0) + y0; % it was x-x0
    plot(x,y,'r')
end

function distance = calc_distance(x1,y1,x2,y2,x3,y3)
n = abs((x2 - x1) * (y1 - y3) - (x1 - x3) * (y2 - y1));
d = sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2);
distance = n ./ d;
return; 
end
