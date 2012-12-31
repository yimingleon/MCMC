%===================================================
% this code is to draw contours from the data
% points are sampled with mixed MCMC method, and saved in sampling.mat
% Yiming Hu, Sep, 2012
%==================================================

load sampleing.mat;

sizeofdata = length(chain);
NoPara = length(chain(1,:));

sorted = sortrows(chain,NoPara-1);
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

base = sorted(1,NoPara-1);
fprintf('Delta chi2 value for sigma1, 2 and 3 are %g,%g and %g\n',sorted(sigma1,NoPara-1)-base,sorted(sigma2,NoPara-1)-base,sorted(sigma3,NoPara-1)-base);

ChainNumber = sortrows(chain(:,NoPara-2));
[sum(ChainNumber==1) sum(ChainNumber==2) sum(ChainNumber==3) sum(ChainNumber==4)]

count
clear
return
