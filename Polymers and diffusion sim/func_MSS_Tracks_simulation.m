function [transDiffResults_all,Track_MSS_all] = func_MSS_Tracks_simulation(track)

[dim1, dim2, dim3] = size(track);

probDim=2; %dimension of tracks. 
plotRes=0;%1 to plot tracks colored by mobility type
peakAlpha=[];%Confidence level for choosing peaks when initially segmenting track. Default: 95 (

warning('off','all')
for i=1:dim1
    for j = 1:dim2
        for k = 1:25
            Track{i,j,k}(:,1) = 1:size(track{i,j,k}.positions,1);
            Track{i,j,k}(:,2:3) = track{i,j,k}.positions;
        end
        Track_MSS_all{i,j}=subfunc_MSS_convertData(squeeze(squeeze(Track(i,j,:))));
        transDiffResults_all{i,j} = basicTransientDiffusionAnalysisv1(Track_MSS_all{i,j},probDim,plotRes,peakAlpha);
        if mod(j,5) == 0
            display(['Analysing row ' num2str(i) ' and column ' num2str(j)]);
        end
    end
end
warning('on','all')


function [Track_MSS]=subfunc_MSS_convertData(trackcell)

size_all=zeros(size(trackcell,1),1); %to get the size of each track
for kk=1:size(trackcell,1)
    size_all(kk,1) = size(trackcell{kk,1},1) ;
end
                                 
Track_MSS=NaN(size(trackcell,1),max(size_all)*8); %this is the way the tracks need to be for DC-MSS
rest=3:8; %the other columns that need to be zeros
for k=1:size(trackcell,1)
    trackmate_track=trackcell{k,1};
    for pp=1:size_all(k)
        xi=1+8*(pp-1);
        yi=2+8*(pp-1);
        resti=rest+8*(pp-1);
        Track_MSS(k,xi)=trackmate_track(pp,2);
        Track_MSS(k,yi)=trackmate_track(pp,3);
        Track_MSS(k,resti(1):1:resti(6))= 0;
    end
end