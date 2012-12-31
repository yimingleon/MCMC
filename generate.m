%=========================================================
% This code is to generate a data set randomlized aruond a sin function.
% purpose is to check the ability of the mixed MCMC method.
% it's just a toy model.
% Yiming Hu, Sep, 2012
%==========================================================

a1 = 1;
a2 = 1;
omega1 = 1;
omega2 = 3;
phi1 = 0;
phi2 = 0;

sigma = 1;

t = 0: 0.01: 2*pi;
y = a1*sin(omega1*t+phi1)+a2*sin(omega2*t+phi2);
y = normrnd(0,sigma,1,length(y))+y;

%plot(t,y)

save('data.mat','y','t','sigma');

%likelihood([a1,omega1],t,y,sigma)-likelihood([a2,omega2],t,y,sigma)
