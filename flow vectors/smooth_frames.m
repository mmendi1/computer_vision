% Only changes made were to make the code compile

function smooth_frames(folder_name,smoothing_factor,start_index,end_index)

frame_1 = single(rgb2gray(read_image(folder_name,start_index)));

stack(:,:,:)=zeros(size(frame_1,1),size(frame_1,2),(end_index - start_index)+1);

for i=start_index:1:end_index
    frame = single(rgb2gray(read_image(folder_name,i)));
    stack(:,:,i)=frame;
end

smoothstack = smoothdata(stack,3,'gaussian','SmoothingFactor',smoothing_factor);

for i=start_index:1:end_index
    imwrite(mat2gray(stack(:,:,i)),fullfile(folder_name,strcat('image_smoothed_',num2str(i),'.png')));
end

end


function I = read_image(folder_name,index)

if(index < 10)
    I = imread(fullfile(folder_name,strcat('frame0',num2str(index),'.png')));
elseif(index < 100)
    I = imread(fullfile(folder_name,strcat('frame',num2str(index),'.png')));
else    
    I = imread(fullfile(folder_name,strcat('frame',num2str(index),'.png')));
end

end

