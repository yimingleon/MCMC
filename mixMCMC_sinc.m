%===================================================
% this code is to realize mixed MCMC
% the method is motivated to explore multiple mode in parameter space
% while using only single chain (parallel tempering not needed.)
% In future version, K-means would be required, and the boundary would be more sofisticated.
% for this toy model, a manul boundary is enough.
% this modified code is specifically fopcused on the validity of mixed MCMC on situations where Likelihood modes have different peaks.
% Yiming Hu, Sep, 2012
%==================================================


function [mixchi] = mixMCMC_sinc;
%load data.mat;

maxlength = 1e5;
%chi2_expect = length(t);
chi2_expect = 0;
% this is just a manually chosen bias. It's for the sake of computing accuracy
cycle = 1e3;

%structure of boundary(x,y,z,...), x represent number of subchain
%y represent parameter 1 :[lower limit, upper limit]
%z represent parameter 2 :[lower limit, upper limit]
%
%in this case, parameter 1 is amplitude, parameter 2 is frequency

%boundary(1,:,:) = [0,2;0,2];
%boundary(2,:,:) = [0,2;2,4];

boundary(1,:,:) = [-2,-1;-2,-1];
boundary(2,:,:) = [-1,1;-2,-1];
boundary(3,:,:) = [1,2;-2,-1];
boundary(4,:,:) = [-2,-1;-1,1];
boundary(5,:,:) = [-1,1;-1,1];
boundary(6,:,:) = [1,2;-1,1];
boundary(7,:,:) = [-2,-1;1,2];
boundary(8,:,:) = [-1,1;1,2];
boundary(9,:,:) = [1,2;1,2];

sig = (boundary(:,:,2)-boundary(:,:,1))/10;
boxsize = size(sig(1,:));
% this is the sigma of the proposed distribution of the normal noise

bound = permute(boundary,[2,3,1]);
% thanks to the structure of the matlab 3dimentional matrix, when I access boundary(1,1,:), I got a size 1 1 2 matrix rather than an array. So I have to rearrange the size.

%subchain(1,:) = [1,1];
%subchain(2,:) = [1,3];

subchain(1,:) = [-1.5,-1.5];
subchain(2,:) = [0,-1.5];
subchain(3,:) = [1.5,-1.5];
subchain(4,:) = [-1.5,0];
subchain(5,:) = [0,0];
subchain(6,:) = [1.5,0];
subchain(7,:) = [-1.5,1.5];
subchain(8,:) = [0,1.5];
subchain(9,:) = [1.5,1.5];



noise = normrnd(0,sig,size(sig));
subchain = subchain + noise;

% if the initial point is too far from the real peak, the algorithm will just ignore the subchain

t=1;
y=1;
variance=1;
%they are no longer needed, but changing the whole code will be to messy, and dont' worth it, so I just randomly set some value.
c = rand;
for i = 1:length(subchain(:,1))
	chi2(i) = likelihood(subchain(i,:),t,y,variance);
end


fun = @(x,y) abs(sinc(x)).*abs(sinc(y));
volume(1) = quad2d(fun,-2,-1,-2,-1);
volume(2) = quad2d(fun,-1,1,-2,-1);
volume(3) = quad2d(fun,1,2,-2,-1);
volume(4) = quad2d(fun,-2,-1,-1,1);
volume(5) = quad2d(fun,-1,1,-1,1);
volume(6) = quad2d(fun,1,2,-1,1);
volume(7) = quad2d(fun,-2,-1,1,2);
volume(8) = quad2d(fun,-1,1,1,2);
volume(9) = quad2d(fun,1,2,1,2);

volume =volume/sum(volume);



%chi2 = chi2 - ones(size(chi2))*max(chi2);
chi2 = chi2 - ones(size(chi2))*chi2_expect;
%here to reduce calculation and to maintain a fix value, using chi2_expected instead of choosing the max chi2 all the time.

likeli = exp(-chi2/2);
normlikeli = likeli/sum(likeli);
p = 0;
leng = length(likeli);
propose = zeros(1,leng);
% this is the number of diff subchains

for i=1:leng
	p = p + normlikeli(i);
	%p_ref is the reference of the normalised probability
	%in the asymmetric Metropolis ratio, it may play an important role
	%The meaning is the normalised probability for the subchain which has the newest point in main chain
	if (c<p), break, end;
end
p_ref = normlikeli(i);
i_ref = i;
propose(i) = propose(i)+1;

chain(1,:) = [subchain(i,:), i,chi2(i),likeli(i)];
 	%points(1,:) = chain(1,:);
% structure of chain:
% parameters; subchain number; chi2 value; likelihood
num_likeli = length(chain(1,:));
num_chi = num_likeli - 1;

count = 1;
n = 2;

tic;
while(n <= maxlength)
	if (~mod(n,cycle))
		fprintf('%d\t',n);
	end

		%if (max(chi2)-min(chi2)<3)
			%normlikeli = likeli/sum(likeli);
		%else
			%normlikeli = exp(-chi2/(max(chi2)-min(chi2)));
			% effectively 'tempering' it when generate candidate.
			% since the MCMC method doesn't care how the candidate is generated, it's OK to do so.
			%normlikeli = normlikeli/sum(normlikeli);
		%end

		%this is for the purpose of test
		%normlikeli = 1/leng*ones(1,leng);
		normlikeli = volume;

	p = 0;
	% p is the accumulating normalized probability
	% this algorithm will pick chain according to their probability
	
	c = rand;
	for i=1:leng
		p = p + normlikeli(i);
		if (c<p)
			%fprintf('%d',i);
			break
	       	end;
	end
	% pick the i^th subchain to generate a candidate
	propose(i) = propose(i)+1;

	while(true)
		q = normrnd(0,sig(i,:),boxsize);
		% here q is the proposal distribution.
		new =  subchain(i,:)+q;
		%new
		%bound(:,2,i)'
		%bound(:,1,i)'
		%i
		if([new<=bound(:,2,i)',new>=bound(:,1,i)']),break,end;
		% the if condition looks quite weird, but it's workable and simple
		% it will be judged true only if every point is within the boundary.
	end


	new_chi2 = likelihood(new,t,y,variance) - chi2_expect;
	new_likeli = exp(-new_chi2/2);

	%if i == i_ref
	%	coefficient = p_ref/normlikeli(i);
	%coefficient = exp(1/2*sum(q.^2./sig(i,:).^2));
	%note that the expression within bracket is positive instead of negative, since the original expression is 1/exp(-1/2...)
	%previous expression has been proved to be wrong	

	%coefficient = p_ref/normlikeli(i);
	%else
		coefficient = p_ref/normlikeli(i)*exp(1/2*sum(q.^2./sig(i,:).^2))*2*pi*sqrt(prod(sig(i,:)));
	%end

	r = new_likeli/chain(n-1,num_likeli)*coefficient;
	%fprintf('%.2g\t',coefficient);
	if (r< rand)	
		chain(n,:) = chain(n-1,:);
		count = count - 1;
	else
		subchain(i,:) = new;
		chi2(i) = new_chi2;
		likeli(i) = new_likeli;
		chain(n,:) = [subchain(i,:), i,chi2(i),likeli(i)];
		p_ref = normlikeli(i);
		i_ref = i;
	end
		%points(n,:) = [new,i,new_chi2,new_likeli];
	%compare with former point

	n=n+1;
	count = count + 1;
end

save sampleing.mat;

sizeofdata = length(chain);
NoPara = length(chain(1,:));


sorted = sortrows(chain,NoPara-1);
sigma1 = round(sizeofdata*0.683);
sigma2 = round(sizeofdata*0.954);
sigma3 = round(sizeofdata*0.9973);
%
%figure
%hold on
%%axis([-1.5,1.5,2.5,7])
%plot(sorted(1:sigma1,1),sorted(1:sigma1,2),'.','Color','b');
%plot(sorted(sigma1+1:sigma2,1),sorted(sigma1+1:sigma2,2),'.','Color','g');
%plot(sorted(sigma2+1:sigma3,1),sorted(sigma2+1:sigma3,2),'.','Color','r');
%xlabel('amplitude');
%ylabel('\omega');
%hold off

%figure
%plot(sorted(1:sigma3,NoPara-1))
%%ylim([-45,-34]);
%xlabel('iteration');
%ylabel('\chi^2');
%title('95% of chi-squared value for mixed MCMC');

ChainNumber = sortrows(chain(:,NoPara-2));
for i=1:leng
	subchain_distribution(i) = sum(ChainNumber==i);
end
subchain_distribution;

count
base = sorted(1,NoPara-1);

%mixchi = [sorted(sigma1,NoPara-1)-base,sorted(sigma2,NoPara-1)-base,sorted(sigma3,NoPara-1)-base];

mixchi(1,:) = subchain_distribution;
%mixchi(2,:) = volume;

mixchi = mixchi/maxlength;


fprintf('the first line is the actual sample probability in different sub-chains\n while the second line is the proposed sample probability in different sub-chains\n');% the third line is the theoretical probability for all sub-chains\n');

toc;
return
clear
