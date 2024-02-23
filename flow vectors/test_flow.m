outfile = 'flow_vectors';
video_obj = VideoWriter(outfile);
open(video_obj);
for i=7:13
    figure;
    demo_optical_flow('Backyard',i,i+1);
    frame = getframe(gcf);
    writeVideo(video_obj,frame);
    fprintf('Frame: %d\n',i)
end
close(video_obj);



