function cumulative_analysis_batch(nFiles)

if nargin<1
    nFiles = 3;
end

defaultDirectory = 'K:\220803\tiff\SARS\kinetics\registered';

% Double expoential analysis
double_exponential = 1;

warning on;

selPath = uigetdir(cd, 'Select folder');
if selPath == 0
    warning("No folder was selected. Analisis is cancelled");
    return
end
selPath = strcat(selPath, '\');
old_dir_cum = cd;

cd(selPath)
files = dir('*SpotSize*.mat');
fileNameCum = {files(:).name};
cd(old_dir_cum);
dir_name = selPath;

if isempty(fileNameCum)
    warning("No file was selected. Analisis is cancelled");
    return
end

cycles = floor(length(fileNameCum)/nFiles);

for ii = 1:cycles
    
    file_cum = {};
    cum_track = [];
    cum_spot = [];

    for jj = 1:nFiles
        file_name = fileNameCum{(ii-1)*nFiles+jj};
        load([selPath file_name]);
        cum_track = [cum_track track_props];
        cum_spot = [cum_spot; spot_record];
        file_cum{jj} = file_name;
        defaultDirectory = selPath;
    end
    
    close all
    cd(old_dir_cum)
    
    [particle_frame, frame_props, surf_coverage] = surface_coverage_stuck(cum_track, metadata);
    for z = 1:length([cum_track(:).bleached])
        if isempty(cum_track(z).bleached)
            cum_track(z).bleached = 1;
        end
    end
    [duration,cum_duration,arrival,cum_arrival,gof_analysis,...
        fit_dissociation,fit_association] = track_analysis(cum_track,metadata);
    
    % Analysis dissociaiton with double exponential as well
    if double_exponential == 1
        [fit_dissociation_2exp, gof_2edx] = doubleExpFit(cum_duration, duration, metadata, fit_dissociation);
    end
    
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
    
    save([file_cum{1}(1:end-4) '_cumulative.mat']);
    
    clear file_cum cum_track cum_spot 
    
end

cd(old_dir_cum)