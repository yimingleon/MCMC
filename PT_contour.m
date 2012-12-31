% this routine is to realize the contour for the MCMC simulation result
% It uses the former simulation result sampling.mat
% also, it use a code from http://www.mathworks.com/matlabcentral/fileexchange/1487
%
% Yiming Hu, Oct, 2012

load PTsampleing;
chains = permute(chains,[2,3,1]);
chains = chains(:,:,1);
mYX = [chains(1,:)',chains(2,:)'];
vXEdge = linspace(0,4,400); 
vYEdge = linspace(0,2,50); 
% be care that in some case it can go to (-1,-1)and (-1,-3)region.
mHist2d = hist2d(mYX,vYEdge,vXEdge); 
  
nXBins = length(vXEdge); 
nYBins = length(vYEdge); 
vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins)); 
vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins)); 
%pcolor(vXLabel, vYLabel,mHist2d); colorbar
%shading flat
%shading interp

%contour(vXLabel,vYLabel,mHist2d,100);
%mesh(vXLabel,vYLabel,mHist2d);

surf(vXLabel,vYLabel,mHist2d); colorbar
shading flat
shading interp
colormap(flipud(bone))
zlim([0,8000]);
clear
