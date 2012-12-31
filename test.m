sig =1:1;
boxsize = size(sig);
mu = zeros(boxsize);
%randomly pick up some number
for i=1:10000
	q=normrnd(mu,sig,boxsize);
	random(i)=exp(-1/2*sum((q-mu).^2./sig.^2));
end

