% this function is written to evaluate the likelihood distribution with a grid search.
% aim: fit a a1.sin(w1.t)+a2.sin(w2.t) data with a.sin(t) likelihood.
% as a reference to the mixMCMC method.
% Yiming Hu, Sep, 2012

load data.mat
i = 1;
amplitude = 0:0.1:4;
frequency = 0:0.1:5;

for A = amplitude
	j = 1;
	for omega = frequency
		grids(i,j) = likelihood([A,omega],t,y,sigma);
		j = j+1;
	end
	i = i+1;
end

grids = grids - ones(size(grids))*min(grids(:));
%grids = exp(-grids/2);
surface(amplitude,frequency,grids')
xlabel('amlitude');
ylabel('frequency');
zlabel('likelihood');
return
