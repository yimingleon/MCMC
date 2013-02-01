% this function is to calcaluate and draw 2 dimensional histogram.
function samples = twoDbin
load PTsampleing.mat;

X = chain(:,1);
Y = chain(:,2);

xbins = -1.5:1:1.5;
% xbins(find(xbins==0))=[];
ybins = -1.5:1:1.5;
% ybins(find(ybins==0))=[];
xNumBins = numel(xbins);
yNumBins = numel(ybins);

Xi = round(interp1(xbins,1:xNumBins,X,'linear','extrap'));
Yi = round(interp1(ybins,1:yNumBins,Y,'linear','extrap'));

Xi = max( min(Xi,xNumBins), 1);
Yi = max( min(Yi,yNumBins), 1);

H = accumarray([Yi(:) Xi(:)],1,[yNumBins xNumBins]);

%figure
%imagesc(xbins, ybins, H), axis on %# axis image
%colormap hot; colorbar
%hold on, plot(X, Y,'b.','MarkerSize',1), hold off

H(2,:) = H(2,:)+H(3,:);
H(3,:) = [];
H(:,2) = H(:,2)+H(:,3);
H(:,3) =[]

samples = reshape(H,1,[]);
samples = samples/sum(samples);
return;
clear
