%%%%%%%% pixel wise fitting temporal or spatial g2 %%%%%%%%%
% Input: 
%     tauaxis       - timelag axis  
%     g2            - spatial or temporal g2(tau) as 1D vector [number of timelags, 1]
%     functiontype  - 't': temperal; 's': spatial
% Output: 
%     fittingresults - fitting results saved as a structure with 3 fields
%         fittingresults.varFit: [beta, rho, tauC]
%         fittingresults.g2Fit: fitted g2(tau) as 1D vector [number of timelags, 1]  
%         fittingresults.R: R^2 score fitting accuracy 
% Example: fittingresults = pixelfitg2(tauaxis, g2, 't')
%
% Note:1) this function is based on g1 model exp(-tau/tauC)^n with n=1; 
%         For n=0.5 or 2, one can directly modify on 'fun', for example: case 't' with n=0.5, 
%         fun = fittype(('1+beta*(rho^2*exp[(-tau/tauC)^0.5]^2+2*rho*(1-rho)*exp[(-tau/tauC)^0.5])'),'indep','tau');
%      2) To fix any parameter in [beta rho tauC] for fitting, 
%         one can set the upper and lower boundary the same, 
%         for example: fix beta=0.5, lb = [0.5 0 10e-3], ub = [0.5 1 20].
%
% last modified by Bingxue Liu, 08/28/2023

function fittingresults = pixelfitg2(tauaxis, g2, functiontype)
    switch functiontype
        case 't'
            fun = fittype(('1+beta*(rho^2*exp(-tau/tauC)^2+2*rho*(1-rho)*exp(-tau/tauC))'),'indep','tau');
        case 's'
            fun = fittype(('1+beta*(rho^2*exp(-tau/tauC)^2+2*rho*(1-rho)*exp(-tau/tauC)+(1-rho)^2)'),'indep','tau');
    end
    lb = [1 0 10e-3]; % [beta rho tauC]
    ub = [1 1 20]; % ms
    Fopts = fitoptions(fun);
    Fopts.Lower = lb;
    Fopts.Upper = ub;
    Fopts.Display='off';
    Fopts.TolFun=1e-6;
    Fopts.Robust='LAR';
    Fopts.StartPoint=[0.5 0.5 1];% [beta rho tauC]
    est = fit(tauaxis',g2,fun,Fopts);
    fittingresults.varFit = [est.beta, est.rho, est.tauC];
    fittingresults.g2Fit = fun(est.beta, est.rho, est.tauC, tauaxis');
    fittingresults.R = (1-sum(abs(g2-fittingresults.g2Fit).^2)./sum(abs(g2-mean(g2)).^2));
end