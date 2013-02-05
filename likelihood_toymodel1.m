function chi2 = likelihood_toymodel1 (params,t,f_of_t,sigma)

% this code will calculate the likelihood
% with input parameters and data
% normal noise is assumed with sigma known
% t is the variable and f_of_t is the function based on t (array)
% Yiming Hu, Sep, 2012

%y = params(1)*sin(params(2)*t);%+params(3))%+params(4)*sin(params(5)*t+params(6));

y = params(1)*sin(params(2)+t);

% the generation of data are actually a combination of two sin function
% but using only one sin funciton to fit
% purpose: make the likelihood function have different peaks with different likelihood.
chi2 = sum(((f_of_t - y)/sigma).^2);
return 
