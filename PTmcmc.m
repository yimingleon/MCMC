%===================================================
% this code is to check mixed MCMC
% using parallel tempering to sample, to see if the result are the same as mixed MCMC.
% Yiming Hu, Sep, 2012
%==================================================

function [chivalue] = PTmcmc;
load data.mat;

% q is the ratio of temperature between 2 neighber chains. according to christian Rover's PhD thesis.
q=6.98;

maxlength = 1e5;
chi2_expect = length(t);
Tmax = 1e3;
No_chain = ceil(log(Tmax)/log(q));
cycle = 1e2;
target = 0.25;
% target acceptance rate
successive = 5;
% number of cycles to successively tweak the same parameter.

%in this case, parameter 1 is amplitude, parameter 2 is frequency

boundary = [0,2;0,4];

sig(1,:) = (boundary(:,2)-boundary(:,1))/50;
% this is the sigma of the proposed distribution of the normal noise

% thanks to the structure of the matlab 3dimentional matrix, when I access boundary(1,1,:), I got a size 1 1 2 matrix rather than an array. So I have to rearrange the size.

for i = 1:No_chain
	chains(i,:,1) = rand(1,2).*(boundary(:,2)-boundary(:,1))'+boundary(:,1)';
	T(i)=q^(i-1);
	sig(i,:) = sig(1,:)*sqrt(T(i));
end

for i = 1:No_chain
	chi2(1,i) = likelihood(chains(i,:,1),t,y,variance) - chi2_expect;
	likeli(1,i) = exp(-chi2(1,i)/2);
end

%chi2 = chi2 - ones(size(chi2))*max(chi2);
%here to reduce calculation and to maintain a fix value, using chi2_expected instead of choosing the max chi2 all the time.

n = 2;

tic;
switch_params = 0;
accept = 0;
boxsize = size(sig(1,:));

while(n < maxlength)
	if (~mod(n,cycle))
		switch_params = switch_params + 1;
		accep_rate = accept/cycle;
		%fprintf('accep_rate is %g\t',accep_rate);
		if switch_params<successive
			sig(1,1) = sig(1,1)*exp(-target+accep_rate);
			%fprintf('amplitude changed to %g\n',sig(1,1));
		else 
			sig(1,2) = sig(1,2)*exp(-target+accep_rate);
			%fprintf('frequency changed to %g\n',sig(1,1));
			if switch_params == 2*successive
				switch_params = 0;
			end
		end

		for i = 1:No_chain
			sig(i,:) = sig(1,:)*sqrt(T(i));
		end
		accept = 0;

		% this is to adapt the variance of the different chains
		% I plan to change variance sperately and successively, so that the process is effectively independant, and could converge to the real situation.
	end


	for i = 1:No_chain
		new =  chains(i,:,n-1)+normrnd(0,sig(i,:),boxsize);
		new_chi2 = likelihood(new,t,y,variance) - chi2_expect;
		new_likeli = exp(-new_chi2/2);

		r = new_likeli/likeli(n-1,i);
		if (r<1)	
			if (r<rand)
				chains(i,:,n) = chains(i,:,n-1);
				chi2(n,i) = chi2(n-1,i);
				likeli(n,i) = likeli(n-1,i);
				if (i==1)
					accept = accept - 1;
				end
			else
				chains(i,:,n) = new;
				chi2(n,i) = new_chi2;
				likeli(n,i) = new_likeli;
			end
		else
			chains(i,:,n) = new;
			chi2(n,i) = new_chi2;
			likeli(n,i) = new_likeli;
		end

		if (i==1)
			accept = accept + 1;
		end
	end
	
	for i=2:No_chain
		w = (likeli(n,i)/likeli(n,i-1))^(1/T(i-1)-1/T(i));
		% probability to accept a swap;
		if w>rand
			new = chains(i,:,n);
			chains(i,:,n) = chains(i-1,:,n);
			chains(i-1,:,n) = new;
			
			new_chi2 = chi2(n,i);
			chi2(n,i) = chi2(n,i-1);
			chi2(n,i-1) = new_chi2;
			
			new_likeli = likeli(n,i);
			likeli(n,i) = likeli(n,i-1);
			liekli(n,i-1) = new_likeli;
		end
	end
	n=n+1;
end

save PTsampleing.mat;

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

base = sorted(1,NoPara);

chivalue = [sorted(sigma1,NoPara)-base,sorted(sigma2,NoPara)-base,sorted(sigma3,NoPara)-base]

toc;
%clear
return;
