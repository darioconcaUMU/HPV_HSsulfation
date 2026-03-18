function [particle_frame, frame_props, surf_coverage, norm_slope] = surface_coverage_stuck(track_props, metadata, label)

% Extract info from metadata and preallocates memory for the following loop
dt = metadata.deltaT;
number_frames = metadata.totalFrames;
total_time = number_frames*dt;
frames = 1:number_frames;
particle_frame = zeros(number_frames,1);

if nargin<3
    label = 'Full video';
end

for i=1:length(track_props)
    start_track = track_props(i).beginning+1;
    end_track = track_props(i).ending+1;
    if end_track > number_frames
        end_track = number_frames;
    end
    particle_frame(start_track:end_track) = particle_frame(start_track:end_track)+1;
end

% Calculates the particle coverage
surf_coverage.mean = mean(particle_frame);
surf_coverage.initial = mean(particle_frame(round(0.15*number_frames):round(0.15*number_frames+10)));
surf_coverage.final = mean(particle_frame(end-13:end-3));
surf_coverage.maximum = max(particle_frame);
surf_coverage.minimum = min(particle_frame);

% Loop to extract the number of particles in each frame and infomation
% about the particle intensity
for i=1:number_frames
    frame_props(i).num_particles = particle_frame;
    temp = particle_frame(2:end)-particle_frame(1:end-1);
    frame_props(i).new = sum(temp(temp>=0));
%     frame_props(i).mean_int = mean(particle_frame{i}(:,5));
%     frame_props(i).std_int = std(double(particle_frame{i}(:,5)));
% This last two will need to extract data from tracks, not track_props, ADD
% IT WHEN YOU HAVE TIME!!
end

% Fits the number of particles per frame with a streight line to determine
% if there is a significant variatio across the video
fit_coverage = fit(frames'*dt,particle_frame,'poly1');
norm_slope = fit_coverage.p1*total_time/surf_coverage.maximum;

figure;
plot(frames*dt,particle_frame, '.');
hold on
plot([1 total_time], fit_coverage([1 total_time]),'r');
hold off
ylabel('Number of particles');
xlabel('Time (s)');
xlim([0 total_time]);
if surf_coverage.maximum > 0
    ylim([0 surf_coverage.maximum*1.1]);
end
title(['Surface coverage - ' label]);

drawnow

% Adds a text to the image, indicating the average surface coverage and the
% variation (as percentage of the maximum coverage) in the particle number
% during the video. A warning is issued if the variation excedes 10%.
txt = ['Average surface coverage: ' num2str(surf_coverage.mean) ' particles/frame'];
if norm_slope < -0.1
    txt_warn = ['Significant reduction in the particle number: ' num2str(norm_slope*100) '%'];
    warning(['Possible significant bleaching detected. '...
        'The number of particles is reduced by ' num2str(norm_slope*100) ...
        '% during the video']);
elseif norm_slope > 0.1
    txt_warn = ['Significant increase in the particle number: ' num2str(norm_slope*100) '%'];
    warning(['The measurement might have been performed before reaching equilibrium. '...
        'The number of particles is increased by ' num2str(norm_slope*100) ...
        '% during the video']);
else
    txt_warn = ['Particle variation: ' num2str(norm_slope*100) '%'];
end
text(number_frames*0.1, surf_coverage.maximum*0.1, {txt, txt_warn}, 'VerticalAlignment','bottom');
 
