function mitoIntensityRaful(pathGFP,pathRFP)
no_points=100;
GFPfilesList=getAllFiles(pathGFP, '*.csv');
RFPfilesList=getAllFiles(pathRFP, '*.csv');
no_files=length(GFPfilesList);
distProfileGFP=zeros(no_points, no_files);
distProfileRFP=zeros(no_points, no_files);
for i=1:length(GFPfilesList)
    cm=hsv(2);
    GFPFile=importdata(GFPfilesList{i});
    intDataGFP=GFPFile.data;
    xGFP=intDataGFP(:,1)/max(intDataGFP(:,1));
    intensityGFP = (intensityGFP - min(intensityGFP))/(max(intensityGFP) - min(intensityGFP));
    %hold on, plot(x,intensity,'Color',cm(j,:));
    [n, bin]=histc(xGFP, linspace(0,1,no_points));%diving each cell's length into 100 equally spaced points-1
    vals=sparse(1:length(xGFP), bin, intensityGFP);%diving each cell's length into 100 equally spaced points-2
    mu=full(sum(vals)./sum(vals~=0));
    distProfileGFP(:,i)=mu;
    RFPFile=importdata(RFPfilesList{i});
    intDataRFP=RFPFile.data;
    xRFP=intDataRFP(:,1)/max(intDataRFP(:,1));
    intensityRFP=intDataRFP(:,2)/max(intDataRFP(:,2));
    %hold on, plot(x,intensity,'Color',cm(j,:));
    [n, bin]=histc(xRFP, linspace(0,1,no_points));%diving each cell's length into 100 equally spaced points-1
    vals=sparse(1:length(xRFP), bin, intensityRFP);%diving each cell's length into 100 equally spaced points-2
    mu=full(sum(vals)./sum(vals~=0));
    distProfileRFP(:,i)=mu;
    
end
meanProfileGFP=zeros(no_points,1);
stdProfileGFP=zeros(no_points,1);
meanProfileRFP=zeros(no_points,1);
stdProfileRFP=zeros(no_points,1);
for j=1:no_points
    meanProfileGFP(j)=mean(distProfileGFP(j,:));
    stdProfileGFP(j)=std(distProfileGFP(j,:));
    meanProfileRFP(j)=mean(distProfileRFP(j,:));
    stdProfileRFP(j)=std(distProfileRFP(j,:));
end
seProfileGFP=stdProfileGFP./sqrt(no_files);
seProfileRFP=stdProfileRFP./sqrt(no_files);
shadedErrorBar(1:no_points,meanProfileGFP,seProfileGFP,'g');
hold on; plot(1:no_points, meanProfileGFP, 'Color', 'g', 'LineWidth',2);
shadedErrorBar(1:no_points,meanProfileRFP,seProfileRFP,'r');
plot(1:no_points, meanProfileRFP, 'Color', 'r', 'LineWidth',2);
axis square, axis([0,100,0,1]);
adjustedCrosscor(meanProfileGFP,meanProfileRFP);
end


