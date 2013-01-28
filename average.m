%===================================================
% This routine is to run mixed MCMC and PT mcmc several times on the same data
% the purpose is to check if the average results are similar in mean and variance
% Yiming Hu, Sep, 2012
%==================================================
times = 100;

for i = 1:times
	chivalue(i,:)=PT_sinc;
	movefile('PTsampleing.mat',sprintf('PT%03d.mat',i));
	mixchi(i,:)= mixMCMC_sinc;
	movefile('sampleing.mat',sprintf('mix%03d.mat',i));
	%weight(i,:) = check_weight;
end
meanchiPT=mean(chivalue)
medianchiPT = median(chivalue);
maxchiPT = chivalue(find(chivalue(:,5)==max(chivalue(:,5))),:);
minchiPT = chivalue(find(chivalue(:,5)==min(chivalue(:,5))),:);

meanchimix=mean(mixchi)
medianchimix = median(mixchi);
maxchimix = mixchi(find(mixchi(:,5)==max(mixchi(:,5))),:);
minchimix = mixchi(find(mixchi(:,5)==min(mixchi(:,5))),:);

%meanchimix=mean(mixchi)
%medianchimix = median(mixchi);
%maxchimix = mixchi(find(mixchi(:,5)==max(mixchi(:,5))),:);
%minchimix = mixchi(find(mixchi(:,5)==min(mixchi(:,5))),:);
%
%meanweight=mean(weight)
%maxweight=max(weight)
%minweight=min(weight)

save average.mat

return
