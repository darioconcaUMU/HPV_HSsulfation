function intensity_histogram

batch_processing = false;

if batch_processing
    
     % Display a dialog box from where to the image-file to be analyzed is chosen.
    % Default directory is User-dependent.
    
    selpath = uigetdir(defaultDirectory, 'Select folder');
    if selpath == 0
        warning("No folder was selected. Analisis is cancelled");
        return
    end
    if ~strcmp(selpath(end),'/')
        selpath = strcat(selpath, '/');
    end
    cd(selpath);
    
    files = dir('*.mat');
    fileName = {files(:).name};
    
    if isempty(fileName)
        warning("No file was selected. Analisis is cancelled");
        return
    end
    
else
    [file,selpath] = ...
        uigetfile('*.mat','Select the MAT file', "G:\Data\Datadrive D\Fouzia\2025\20250102\Kinetics\Test");
    if file == 0
       warning("No file was selected. Analisis is cancelled");
       return
    end
    fileName{1} = file;
end

spot_intensity = [];
particle_intensity = [];
old = cd;

for i=1:length(fileName)
    
    clear spot_record metadata tracks track_props duration cum_duration...
       arrival cum_arrival particle_frame frame_props surf_coverage
   
   cd(selpath)
   
   currentfile = fileName{i};
   load(currentfile);
   
   spot_intensity = [spot_intensity spot_record(:,8)];
   start_int = [track_props(:).initial_prominence];
   particle_intensity = [particle_intensity start_int];
   
%    % This code will calculate the intensity. The prominence is much more
%    interesting but it requires the analysis with the last version of the
%    software 3.1
% 
%    spot_intensity = [spot_intensity spot_record(:,5)];
%    
%    for j=1:length(tracks)
%        start_int(j) = mean(tracks{j}(1:3,5));
%    end
%    
%    particle_intensity = [particle_intensity start_int];
    
end 

figure
histogram(spot_intensity);
title('Histrogram spot intensities');
xlabel('Intensity spot');
ylabel('Counts');
set(gca,'YScale','log')

figure
histogram(particle_intensity);
title('Histrogram initial intensity detected particles');
xlabel('Intensity particle');
ylabel('Counts');
set(gca,'YScale','log')

cd(old);