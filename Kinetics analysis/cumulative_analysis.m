function cumulative_analysis

% Funciton to use to merge data from several videos and perform the
% cumulative analysis on them

defaultDirectory = 'K:\220803\tiff\SARS\kinetics\registered';
file_cum = {};
cum_track = [];
cum_spot = [];
index = 1;
warning on;

% Double expoential analysis
double_exponential = 1;

while true
    [file_name,pathName]=uigetfile('*.mat','Select the Matlab files. Press cancel to stop'...
        ,defaultDirectory);
    
    if file_name == 0
        break
    elseif any(strcmp(file_cum,file_name))
        warning('The file was already loaded, please choose a differnt file or cancel to end the selection');
    else
        load([pathName file_name]);
        cum_track = [cum_track track_props];
        cum_spot = [cum_spot; spot_record];
        file_cum{index} = file_name;
        defaultDirectory = pathName;
        index = index+1;
        dir_name = pathName;
    end
end

close all

[particle_frame, frame_props, surf_coverage] = surface_coverage_stuck(cum_track, metadata);
for i = 1:length([cum_track(:).bleached])
    if isempty(cum_track(i).bleached)
        cum_track(i).bleached = 1;
    end
end
[duration,cum_duration,arrival,cum_arrival] = track_analysis(cum_track,metadata);

% Analysis dissociaiton with double exponential as well
if double_exponential == 1
    [fit_dissociation_2exp, gof_2edx] = doubleExpFit(cum_duration, duration, metadata, fit_dissociation);
end
    

old_dir = cd;
cd(dir_name)
mkdir(strcat('Cumulative_', file_cum{1}(1:end-4)));
cd(strcat('Cumulative_', file_cum{1}(1:end-4)));

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    %FigName = get(FigHandle, 'Name');
    savefig(FigHandle, strcat('cumulative_', file_cum{1}(1:end-4), '_', num2str(iFig), '.fig'));
    saveas(FigHandle, strcat('cumulative_', file_cum{1}(1:end-4), '_', num2str(iFig), '.png'));
end

close all

save([file_cum{1}(1:end-4) '_cumulative.mat'])
cd(old_dir)

