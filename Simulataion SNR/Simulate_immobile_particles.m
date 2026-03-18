% simulate_particles.m
% -------------------------------------------------------------------------
% Simulates imaging immobile fluorescent particles on a surface and reports
% three Signal‑to‑Noise Ratios (SNRs), each defined as:
%            average particle intensity  /  σ_totalNoise
% where σ_totalNoise is the standard deviation of *all* noise added on top of
% the noise‑free reference image (background fluctuation + acquisition noise).
%   • Theoretical    – uses input parameters only.
%   • EffectiveLow   – measured in frame 1.
%   • EffectiveHigh  – measured from the full movie (all frames).
%
% OUTPUTS (names include the theoretical SNR value, e.g. “…_SNR23p50…”)
%   1) simulation_data_SNRXXpXX[_N].mat     – particle data + metadata + SNRs
%   2) simulation_video_SNRXXpXX[_N].tif    – synthetic movie
% -------------------------------------------------------------------------
% USER‑CONFIGURABLE PARAMETERS
params.imgWidth            = 704;      % image width  (pixels)
params.imgHeight           = 704;      % image height (pixels)
params.numParticles        = 150;      % number of particles
params.meanIntensity       = 170;     % mean particle peak intensity  (a.u.)
params.intensityStd        = 0;     % standard deviation of particle intensities
params.psfSigma            = 1.2;      % Gaussian PSF sigma (pixels)
params.backgroundMean      = 145;      % mean uniform background  (a.u.)
params.backgroundStd       = 0;       % std of background fluctuations (a.u.)
params.noiseStd            = 8.3;       % acquisition (read) noise std (a.u.)
params.numFrames           = 100;      % total number of frames in the movie
params.bitDepth            = 16;       % 8 or 16; determines output data type
params.randomSeed          = 42;       % for reproducibility (set [] to random)

% -------------------------------------------------------------------------
% INITIALISATION
if ~isempty(params.randomSeed); rng(params.randomSeed); end
switch params.bitDepth
    case 8,  imgClass = 'uint8';
    case 16, imgClass = 'uint16';
    otherwise, error('bitDepth must be 8 or 16');
end
maxVal = double(intmax(imgClass));

% Generate sub‑pixel particle positions (double)
margin     = ceil(3*params.psfSigma);
particleX  = margin + (params.imgWidth  - 2*margin)  * rand(params.numParticles,1);
particleY  = margin + (params.imgHeight - 2*margin)  * rand(params.numParticles,1);

% Generate particle peak intensities
particleI  = max(params.meanIntensity + params.intensityStd*randn(params.numParticles,1),0);

% -------------------------------------------------------------------------
% BUILD NOISE‑FREE REFERENCE IMAGE (static signal)
referenceFrame = zeros(params.imgHeight, params.imgWidth, 'double');
radius = ceil(3*params.psfSigma);
for p = 1:params.numParticles
    xMin = max(floor(particleX(p)) - radius, 1);
    xMax = min(floor(particleX(p)) + radius, params.imgWidth);
    yMin = max(floor(particleY(p)) - radius, 1);
    yMax = min(floor(particleY(p)) + radius, params.imgHeight);
    [xLocal, yLocal] = meshgrid(xMin:xMax, yMin:yMax);
    dx = xLocal - particleX(p);
    dy = yLocal - particleY(p);
    psf = exp(-(dx.^2 + dy.^2) / (2*params.psfSigma^2));
    referenceFrame(yMin:yMax, xMin:xMax) = referenceFrame(yMin:yMax, xMin:xMax) + ...
        particleI(p) * psf;
end
referenceFrame = max(min(referenceFrame, maxVal), 0);

% -------------------------------------------------------------------------
% SIMULATE FRAMES (add background + noise)
videoStack = zeros(params.imgHeight, params.imgWidth, params.numFrames, imgClass);
for f = 1:params.numFrames
    frame = referenceFrame + params.backgroundMean + params.backgroundStd*randn();
    frame = frame + params.noiseStd * randn(params.imgHeight, params.imgWidth);
    frame = max(min(frame, maxVal), 0);
    videoStack(:,:,f) = cast(frame, imgClass);
end

% -------------------------------------------------------------------------
% SNR CALCULATIONS
sigmaTotal_theoretical = sqrt(params.backgroundStd^2 + params.noiseStd^2);
SNR_theoretical        = params.meanIntensity / sigmaTotal_theoretical;

% Frame‑1 noise image (double)
noiseFrame1 = double(videoStack(:,:,1)) - referenceFrame - params.backgroundMean;
SNR_effectiveLow = mean(interp2(double(videoStack(:,:,1)), particleX, particleY, 'linear')) / std(noiseFrame1(:));

% All‑frames noise array
noiseAll = double(videoStack) - referenceFrame - params.backgroundMean;
SNR_effectiveHigh = mean(mean(interp2(double(videoStack(:,:,1)), particleX, particleY, 'linear'))) ; % placeholder, will calc below
measI_all = zeros(params.numParticles, params.numFrames);
for f = 1:params.numFrames
    measI_all(:,f) = interp2(double(videoStack(:,:,f)), particleX, particleY, 'linear');
end
SNR_effectiveHigh = mean(mean(measI_all,2)) / std(noiseAll(:));

% -------------------------------------------------------------------------
% AUTOGENERATE FILENAMES WITH THEORETICAL SNR
strSNR  = strrep(sprintf('%.2f', SNR_theoretical), '.', 'p');
baseVideo = ['simulation_video_SNR' strSNR '.tif'];
baseData  = ['simulation_data_SNR' strSNR '.mat'];
videoName = baseVideo; dataName = baseData; suffix = 1;
while exist(videoName, 'file') || exist(dataName, 'file')
    suffix = suffix + 1;
    videoName = [baseVideo(1:end-4) '_' num2str(suffix) '.tif'];
    dataName  = [baseData(1:end-4)  '_' num2str(suffix) '.mat'];
end

% -------------------------------------------------------------------------
% WRITE VIDEO
for f = 1:params.numFrames
    if f == 1
        imwrite(videoStack(:,:,f), videoName, 'Compression', 'none');
    else
        imwrite(videoStack(:,:,f), videoName, 'WriteMode', 'append', 'Compression', 'none');
    end
end

% -------------------------------------------------------------------------
% SAVE DATA & METADATA (INCLUDE SNRs)
metadata              = params;
metadata.creationDate = datestr(now);
metadata.SNR_theoretical   = SNR_theoretical;
metadata.SNR_effectiveLow  = SNR_effectiveLow;
metadata.SNR_effectiveHigh = SNR_effectiveHigh;

save(dataName, 'particleX', 'particleY', 'particleI', 'metadata');

% -------------------------------------------------------------------------
% REPORT
fprintf('Simulation complete.\n');
fprintf('  > Data saved to   %s\n', dataName);
fprintf('  > Video saved to  %s (%d frames)\n', videoName, params.numFrames);
fprintf('Theoretical SNR: %.2f\n', SNR_theoretical);
fprintf('Effective  SNR (frame1): %.2f\n', SNR_effectiveLow);
fprintf('Effective  SNR (all)   : %.2f\n', SNR_effectiveHigh);
