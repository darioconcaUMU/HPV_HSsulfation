function cumulative_analysis_batch_2v0(nFiles)

if nargin<1
    nFiles = 3;
end

defaultDirectory = 'G:\Data\Datadrive D\Fouzia\2024\20241218\BouncingonPEGPOPC\TIFF';

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
        if exist('fraction','var')
            num_fractions = length(fraction);
        else
            num_fraction = 0;
        end
        cum_track = [cum_track track_props];
        cum_spot = [cum_spot; spot_record];
        file_cum{jj} = file_name;
        defaultDirectory = selPath;

    end
    
    close all
    cd(old_dir_cum)

    if isfield(parameters, 'intensityFraction')
        [fraction, parameters] = analysis_fract_intensity(parameters, intensity_fraction, cum_track, cum_spot, metadata);
    end
    
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
        parameters.double_exp = 1;
        if length(duration)>5
            [fit_dissociation_2exp, gof_2edx] = doubleExpFit(cum_duration, duration, metadata, fit_dissociation);
            
        else
            warning('Not enough points to perform 2 exponent fit');
        end
    else
        parameters.double_exp = 0;
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

    if exist('fraction','var')
        for i=1:length(fraction)


        end
    end
    
end

cd(old_dir_cum)