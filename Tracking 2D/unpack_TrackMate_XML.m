function tracks = unpack_TrackMate_XML(xmlfile)
% Unpacks a TrackMate XML, output: tracks{n} = [frame x y]
% Usage: tracks = unpack_TrackMate_XML('yourfile.xml')

if ~endsWith(xmlfile,'xml')
    xmlfile = strcat(xmlfile,'.xml');
end
% Read XML document
doc = xmlread(xmlfile);

% Find all 'particle' elements (each is a track)
particles = doc.getElementsByTagName('particle');
nParticles = particles.getLength;
tracks = cell(nParticles, 1);

for i = 0:nParticles-1
    thisParticle = particles.item(i);
    detections = thisParticle.getElementsByTagName('detection');
    nDetections = detections.getLength;
    dat = nan(nDetections, 3); % [frame x y]
    for j = 0:nDetections-1
        det = detections.item(j);
        t = str2double(det.getAttribute('t'));
        x = str2double(det.getAttribute('x'));
        y = str2double(det.getAttribute('y'));
        dat(j+1, :) = [t, x, y];
    end
    tracks{i+1} = dat;
end

tracks = tracks';

end
