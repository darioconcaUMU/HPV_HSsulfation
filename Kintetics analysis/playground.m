selPath = uigetdir(cd, 'Select folder');
if selPath == 0
    warning("No folder was selected. Analisis is cancelled");
    return
end
selPath = strcat(selPath, '\');

old_dirMAT = cd;
cd(selPath);

files = dir('*SpotSize*.mat');
fileNameMAT = {files(:).name};

if isempty(fileNameMAT)
    warning("No file was selected. Analisis is cancelled");
    return
end


for i=1:length(fileNameMAT)
    
    cd(selPath);
    
    close all
    clear spot_record metadata tracks track_props duration cum_duration...
        arrival cum_arrival particle_frame frame_props surf_coverage
    
    currentfileMAT = fileNameMAT{i};
    
    load(currentfileMAT);
    
    cd(old_dirMAT);
    
    selpath = selPath;
        
    % Analyze remaining tracks
    try
    [duration,cum_duration,arrival,cum_arrival,gof_analysis,...
        fit_dissociation,fit_association] = track_analysis(track_props,metadata);
    catch
        a = 1;
    end

    cd(selpath);
    
    clear selPath
    
    %save(strcat(matFileName, '.mat'));
    save(currentfileMAT);
    
    selPath = selpath;

end

close all
cd(old_dirMAT);