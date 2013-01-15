%===================================================
% This routine is to run mixed MCMC and PT mcmc several times on the same data
% the purpose is to check if the average results are similar in mean and variance
% Yiming Hu, Sep, 2012
%==================================================
times = 10;

for i = 1:times
	chivalue(i,:)=PTmcmc;
	%mixchi(i,:)= mixMCMC;
end
meanchiPT=mean(chivalue)
%meanchimix=mean(mixchi)

save average.mat

return
