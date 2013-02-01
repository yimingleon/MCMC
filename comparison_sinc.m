% this routine will load the previously saved .mat files, and compare the output of mixed MCMC and parallel tempering MCMC for the problem of 2D sinc likelihood
% Yiming Hu, Jan, 2013
%============================================================
load average.mat
%to read the variable times, as the total number of .mat files.

%here I will use loop instead of i because that there will be variable i in the .mat files, and may induce conflict
for loop = 1:times
	fprintf('%d\t',loop);
	load(sprintf('mix%03d.mat',loop));
	sizeofdata = length(chain);
	NoPara = length(chain(1,:));

	sigma1 = round(sizeofdata*0.683);
	sigma2 = round(sizeofdata*0.954);
	sigma3 = round(sizeofdata*0.9973);

	sorted = sortrows(chain,NoPara-1);
	base = sorted(1,NoPara-1);

	delta_chimix(loop,:) = [sorted(sigma1,NoPara-1)-base,sorted(sigma2,NoPara-1)-base,sorted(sigma3,NoPara-1)-base];
	%clearvars -except delta_chimix,delta_chiPT,loop,times,sizeofdata,NoPara,sigma1,sigma2,sigma3
	load(sprintf('PT%03d.mat',loop));
	delta_chiPT(loop,:) = [sorted(sigma1,NoPara)-base,sorted(sigma2,NoPara)-base,sorted(sigma3,NoPara)-base];
	%clearvars -except delta_chimix,delta_chiPT,loop,times,sizeofdata,NoPara,sigma1,sigma2,sigma3

end

clearvars -except delta_chimix delta_chiPT
