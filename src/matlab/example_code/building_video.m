% Source : https://se.mathworks.com/matlabcentral/answers/384848-get-figures-and-use-them-to-build-a-video-avi

for i = 1:N
    figure(1)  
    imshow(processo(:,:,1,i))
      hold on
      plot(X,Y,'o')
      plot(X0,Y0,'o')
      plot(X1,Y1,'o')
      plot(X2,Y2,'o')
      plot(X3,Y3,'o')
      hold off
      F(i) = getframe(gcf) ;
      drawnow
    end
  % create the video writer with 1 fps
  writerObj = VideoWriter('myVideo.avi');
  writerObj.FrameRate = 10;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);