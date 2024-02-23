function F = FundamentalRefinement(leftFile,rightFile)

I1 = imread(leftFile);
I2 = imread(rightFile);
I1 = rgb2gray(I1);
I2 = rgb2gray(I2);

ptsLeft  = detectSURFFeatures(I1); 
ptsRight = detectSURFFeatures(I2);
[featuresLeft,  validPtsLeft]  = extractFeatures(I1,  ptsLeft);
[featuresRight, validPtsRight] = extractFeatures(I2, ptsRight);
indexPairs = matchFeatures(featuresLeft, featuresRight);
matchedLeft  = validPtsLeft(indexPairs(:,1));
matchedRight = validPtsRight(indexPairs(:,2));

% Select 40 random SURF keypoints
cL = {};
cR = {};
indices = randi([1 length(matchedLeft)], 40,1);
for i=1:40
    cL{end+1} = matchedLeft(indices(i));
    cR{end+1} = matchedRight(indices(i));
end
 
iterations = 1;
max_iterations = 100;
Fs = {};
counts = [];

while(iterations < max_iterations)

    % select 8 random points
    u = randi([1 40], 8,1);
    left_points = [[cL{u(1)}.Location];[cL{u(2)}.Location];[cL{u(3)}.Location];[cL{u(4)}.Location];[cL{u(5)}.Location];[cL{u(6)}.Location];[cL{u(7)}.Location];[cL{u(8)}.Location]];
    right_points = [[cR{u(1)}.Location];[cR{u(2)}.Location];[cR{u(3)}.Location];[cR{u(4)}.Location];[cR{u(5)}.Location];[cR{u(6)}.Location];[cR{u(7)}.Location];[cR{u(8)}.Location]];
    
    % dont need normalization
    
    x1s = left_points(:,1);
    y1s = left_points(:,2);
    x2s = right_points(:,1);
    y2s = right_points(:,2);
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
    
    A = [x2_norm.*x1_norm x2_norm.*y1_norm x2_norm y2_norm.*x1_norm y2_norm.*y1_norm y2_norm x1_norm y1_norm ones(8,1); zeros(1,9)]; % least squares
    [U1,S1,V1] = svd(A);
    F_normalized = reshape(V1(:,end),3,3)';
    F_denormalized = (M2')*F_normalized*(M1); % denormalization

    [U2,S2,V2] = svd(F_denormalized); 
    S2(3,3) = 0; % rank resolution
    F = U2*S2*(V2');
    
    count = 0;
    for i=1:40 
        if any(u(:) == i)
            continue
        elseif in_consensus(cL{i}, cR{i}, F)
            count = count + 1;
        else
            continue
        end
    end
    
    Fs{end+1} = F/norm(F);
    counts(end+1) = count;
    
    iterations = iterations + 1;
end


[M,I] = max(counts);
F = Fs{1,I};
end


function boolean = in_consensus(left_keypoint, right_keypoint, F)
    boolean = false;
    x1 = [left_keypoint.Location 1]';
    x2 = [right_keypoint.Location 1]';
    total = x2' * (F * x1);
    if (abs(total) < 0.005)
        boolean = true;
    end
end