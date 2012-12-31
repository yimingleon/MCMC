%===================================================
% this code is to draw contours from the data
% points are sampled with mixed MCMC method, and saved in sampling.mat
% Yiming Hu, Sep, 2012
%==================================================

load PTsampleing.mat;

chains = permute(chains,[2,3,1]);
chains = chains(:,:,1);
% this is to get only the lowest temperature chain, and make a 3D matrix to 2D.
chain = [chains',chi2(:,1)];
sizeofdata = length(chain);
NoPara = length(chain(1,:));

sorted = sortrows(chain,NoPara);
sigma1 = round(sizeofdata*0.683);
sigma2 = round(sizeofdata*0.954);
sigma3 = round(sizeofdata*0.9973);

figure
hold on
plot(sorted(1:sigma1,1),sorted(1:sigma1,2),'.','Color','b');
plot(sorted(sigma1+1:sigma2,1),sorted(sigma1+1:sigma2,2),'.','Color','g');
plot(sorted(sigma2+1:sigma3,1),sorted(sigma2+1:sigma3,2),'.','Color','r');
xlabel('amplitude');
ylabel('\omega');
hold off

base = sorted(1,NoPara);
fprintf('Delta chi2 value for sigma1, 2 and 3 are %g,%g and %g\n',sorted(sigma1,NoPara)-base,sorted(sigma2,NoPara)-base,sorted(sigma3,NoPara)-base);

clear
return
