
dims = [2 1 3]; % Video dimensions (One video per cell)
dims = 2;

for dim = dims

% Set up the figure and video properties
vidObj = VideoWriter(fullfile(output_dir, [subj_id '_tussim_skull_3D_' transducer '_dim' int2str(dim) '.mp4']), 'MPEG-4');
vidObj.FrameRate = 10; % Adjust the frame rate as desired
open(vidObj);

numFrames = size(t1_img, dim);
t1_img = flip(t1_img, 1);

if dim == 1
    order = [1 2 3];
elseif dim == 2
    order = [2 1 3];
elseif dim == 3
    order = [3 1 2];
end
% order = circshift(1:3, -dim+1);
t1_new = permute(t1_img, order);
p_new = permute(p, order);

% Get colorbar limits
color_limits = [0, max(p(:))*1e-6];

% Loop over each frame and create the plot
for frameIndex = 1:numFrames % Replace numFrames with the total number of frames in your sequence

%     fig = figure;
    fig = figure('visible','off');
    ax1 = axes;
    imagesc(ax1, imrotate(squeeze(t1_new(frameIndex, :,:)),90), [50,500]);
%     imagesc(ax1, imrotate(squeeze(t1_new(mx, :,:)),90), [50,500]);
    hold all;
    ax2 = axes;
    im2 = imagesc(ax2, imrotate(squeeze(p_new(frameIndex, :,:))*1e-6,90));
%     im2 = imagesc(ax2, imrotate(squeeze(p_new(mx, :,:))*1e-6,90));
    im2.AlphaData = 0.5;
    linkaxes([ax1,ax2]); ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = [];
    colormap(ax1,'gray')
    colormap(ax2,'turbo')
    set([ax1,ax2],'Position',[.17 .11 .685 .815]);
    cb2 = colorbar(ax2,'Position',[.85 .11 .0275 .815]);
    clim(color_limits)
    xlabel(cb2, '[MPa]');
    title(ax1,'Acoustic Pressure Amplitude')
    
    % Save the current frame as an image file
    frameFile = sprintf('frame_%04d.png', frameIndex);
    saveas(fig, frameFile);
    
    % Add the frame to the video object and delete the image file
    frameData = imread(frameFile);
    writeVideo(vidObj, frameData);
    delete(frameFile);

    disp(['Frame ' int2str(frameIndex) '/' int2str(numFrames) ' processed'])
end

disp('Complete!')

% Close the video object
close(vidObj);

end
