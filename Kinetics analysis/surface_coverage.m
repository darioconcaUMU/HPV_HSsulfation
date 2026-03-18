function [particle_frame, frame_props, surf_coverage] = surface_coverage(spot_record, metadata)

% Extract info from metadata and preallocates memory for the following loop
dt = metadata.deltaT;
number_frames = metadata.totalFrames;
total_time = number_frames*dt;
frames = 1:number_frames;
particle_frame = cell(number_frames, 1);

%spot_record = spot_record(spot_record(:,6)==0,:);

% Loop to extract the number of particles in each frame and infomation
% about the particle intensity
for i=1:number_frames
    particle_frame{i} = spot_record(spot_record(:,2) == frames(i),:);
    frame_props(i).num_particles = size(particle_frame{i},1);
    temp = particle_frame{i}(:,1)-particle_frame{i}(:,2);
    frame_props(i).new = sum(temp>=0);
    frame_props(i).mean_int = mean(particle_frame{i}(:,5));
    frame_props(i).std_int = std(double(particle_frame{i}(:,5)));
end

% Calculates the particle coverage
surf_coverage.mean = mean([frame_props(:).num_particles]);
surf_coverage.initial = frame_props(3).num_particles;
surf_coverage.final = frame_props(end-3).num_particles;
surf_coverage.maximum = max([frame_props(:).num_particles]);
surf_coverage.minimum = min([frame_props(:).num_particles]);

% Fits the number of particles per frame with a streight line to determine
% if there is a significant variatio across the video
fit_coverage = fit(frames'*dt,[frame_props(:).num_particles]','poly1');
norm_slope = fit_coverage.p1*total_time/surf_coverage.maximum;

% Plots the number of particles against the elapsted time.
figure;
plot(frames*dt,[frame_props(:).num_particles], '.');
hold on
plot([1 total_time], fit_coverage([1 total_time]),'r');
hold off
ylabel('Number of particles');
xlabel('Time (s)');
xlim([0 total_time]);
ylim([0 surf_coverage.maximum]);
title('Surface coverage');

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
 
