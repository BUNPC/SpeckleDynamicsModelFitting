%%%%%%%% pixel wise fitting temporal or spatial g2 with mixed dynamics%%%%%%%%%
% Input: 
%     tauaxis       - timelag axis  
%     g2            - spatial or temporal g2(tau) as 1D vector [number of timelags, 1]
%     functiontype  - 't': temperal; 's': spatial
% Output: 
%     fittingresults - fitting results saved as a structure with 3 fields
%         fittingresults.varFit: [beta, rho1, rho2, tauC1, tauC2]
%         fittingresults.g2Fit: fitted g2(tau) as 1D vector [number of exposure times, 1]  
%         fittingresults.R: R^2 score fitting accuracy 
% Example: fittingresults = pixelfitg2_2tau(tauaxis, g2, 't')
%
% Note:1) this function is based on g1 model exp(-(tau/tauC)^n) with n=1; 
%         For n=0.5 or 2, one can directly modify on 'fun', for example:
%         case 's' with n=2 for tauC1, n=0.5 for tauC2
%         fun = fittype(('1+beta*(rho1*exp(-(tau/tauC1)^2)+rho2*exp(-(tau/tauC2)^0.5)+(1-rho1-rho2))^2'),'indep','tau');
%      2) To fix any parameter in [beta rho1 rho2 tauC1 and tauC2] for fitting, 
%         one can set the upper and lower boundary the same, 
%         for example: fix beta=0.5, tauC1 = 5, lb = [0.5, 0, 0, 5, 5]; ub = [0.5, 1, 1, 5, 10e3];
%
% last modified by Bingxue Liu, 08/28/2023

function fittingresults = pixelfitg2_2tau(tauaxis, g2, functiontype)
    switch functiontype
        case 't'
            fun = fittype(('1+beta*((rho1*exp(-tau/tauC1)+rho2*exp(-tau/tauC2)+(1-rho1-rho2))^2-(1-rho1-rho2)^2)'),'indep','tau');
        case 's'
            fun = fittype(('1+beta*(rho1*exp(-tau/tauC1)+rho2*exp(-tau/tauC2)+(1-rho1-rho2))^2'),'indep','tau');
    end
    lb = [1 0 0 10 100]; % [beta rho1, rho2, tauC1 tauC2]
    ub = [1 0.5 0.5 10 10000]; % ms
    Fopts = fitoptions(fun);
    Fopts.Lower = lb;
    Fopts.Upper = ub;
    Fopts.Display='off';
    Fopts.TolFun=1e-8;
    Fopts.Robust='LAR';
    Fopts.StartPoint=[0.5 0.3 0.3 10 1000];%  [beta rho1, rho2, tauC1, tauC2]
    est = fit(tauaxis',g2,fun,Fopts);
    fittingresults.varFit = [est.beta, est.rho1, est.tauC1, est.rho2, est.tauC2];
    fittingresults.g2Fit = fun(est.beta, est.rho1, est.rho2, est.tauC1, est.tauC2, tauaxis');
    fittingresults.R = (1-sum(abs(g2-fittingresults.g2Fit).^2)./sum(abs(g2-mean(g2)).^2));
end