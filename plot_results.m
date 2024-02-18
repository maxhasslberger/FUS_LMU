% Load results file in advance!

disp("Max. ISPPA: ", Isppa)
disp("ISPPA at Focus: ", Isppa_focus)
disp("Max. ISPTA: ", Ispta)

figure;
ax1 = axes;
imagesc(ax1, imrotate(squeeze(t1_img(mx,:,:)),90), [50,500]);
hold all;
ax2 = axes;
im2 = imagesc(ax2, imrotate(squeeze(p(mx,:,:))*1e-6,90));
im2.AlphaData = 0.5;
linkaxes([ax1,ax2]); ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = [];
colormap(ax1,'gray')
colormap(ax2,'turbo')
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb2 = colorbar(ax2,'Position',[.85 .11 .0275 .815]);
xlabel(cb2, '[MPa]');
title(ax1,'Acoustic Pressure Amplitude')
% saveas(gcf, fullfile(output_dir, [subj_id '_tussim_skull_3D_' transducer '_sag.jpg']));

figure;
ax1 = axes;
imagesc(ax1, imrotate(squeeze(t1_img(:,my,:)),90), [50,500]);
hold all;
ax2 = axes;
im2 = imagesc(ax2, imrotate(squeeze(p(:,my,:))*1e-6,90));
im2.AlphaData = 0.5;
linkaxes([ax1,ax2]); ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = [];
colormap(ax1,'gray')
colormap(ax2,'turbo')
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb2 = colorbar(ax2,'Position',[.85 .11 .0275 .815]);
xlabel(cb2, '[MPa]');
title(ax1,'Acoustic Pressure Amplitude')
% saveas(gcf, fullfile(output_dir, [subj_id '_tussim_skull_3D_' transducer '_cor.jpg']));

figure;
ax1 = axes;
imagesc(ax1, imrotate(squeeze(t1_img(:,:,mz)),90), [50,500]);
hold all;
ax2 = axes;
im2 = imagesc(ax2, imrotate(squeeze(p(:,:,mz))*1e-6,90));
im2.AlphaData = 0.5;
linkaxes([ax1,ax2]); ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = [];
colormap(ax1,'gray')
colormap(ax2,'turbo')
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb2 = colorbar(ax2,'Position',[.85 .11 .0275 .815]);
xlabel(cb2, '[MPa]');
title(ax1,'Acoustic Pressure Amplitude')
% saveas(gcf, fullfile(output_dir, [subj_id '_tussim_skull_3D_' transducer '_ax.jpg']));