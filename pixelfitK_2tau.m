%%%%%%%% pixel wise fitting temporal or spatial contrast K with mixed dynamics%%%%%%%%%
% Input: 
%     T             - multiple exposure times  
%     K             - spatial or temporal K(T) as 1D vector [number of exposure times, 1]
%     functiontype  - 't': temperal; 's': spatial
% Output: 
%     fittingresults - fitting results saved as a structure with 3 fields
%         fittingresults.varFit: [beta, rho1, rho2, tauC1, tauC2]
%         fittingresults.g2Fit: fitted g2(tau) as 1D vector [number of exposure times, 1]  
%         fittingresults.R: R^2 score fitting accuracy 
% Example: fittingresults = pixelfitK_2tau(T, K, 't')
%
% Note:1) this function is based on g1 model exp(-(tau/tauC)^n) with n=1; 
%         For n=0.5 or 2, one can directly modify on 'fun', for example:
%         case 's' with n=2 for tauC1, n=0.5 for tauC2
%         fun = @(beta, rho1, rho2, tauC1, tauC2, T) (integral(@(x) 2*beta./T.*(1-x./T).*(rho1*exp(-(x./tauC1)^2)+rho2*exp(-(x./tauC2)^0.5)+1-rho1-rho2).^2, 0, T, 'ArrayValued', true)).^0.5;
%      2) To fix any parameter in [beta rho1 rho2 tauC1 and tauC2] for fitting, 
%         one can set the upper and lower boundary the same, 
%         for example: fix beta=0.5, tauC1 = 5, lb = [0.5, 0, 0, 5, 5]; ub = [0.5, 1, 1, 5, 10e3];
%
% last modified by Bingxue Liu, 08/28/2023

function fittingresults = pixelfitK_2tau(T, K, functiontype)
    switch functiontype
        case 't'
            fun = @(beta, rho1, rho2, tauC1, tauC2, T) (integral(@(x) 2*beta./T.*(1-x./T).*((rho1*exp(-x./tauC1)+rho2*exp(-(x./tauC2))+1-rho1-rho2).^2-(1-rho1-rho2).^2), 0, T, 'ArrayValued', true)).^0.5;
        case 's'
            fun = @(beta, rho1, rho2, tauC1, tauC2, T) (integral(@(x) 2*beta./T.*(1-x./T).*(rho1*exp(-x./tauC1)+rho2*exp(-(x./tauC2))+1-rho1-rho2).^2, 0, T, 'ArrayValued', true)).^0.5;
    end
    funs = @(beta, rho1, rho2, tauC1, tauC2, T) arrayfun(fun, beta, rho1, rho2, tauC1, tauC2, T);

    % Define the objective function for lsqcurvefit
    objfun = @(varFit,T) funs(varFit(1)*ones(length(T),1),varFit(2)*ones(length(T),1),varFit(3)*ones(length(T),1),varFit(4)*ones(length(T),1),varFit(5)*ones(length(T),1),T);

    % Fit the function using lsqcurvefit
    x0 = [1, 0.3, 0.3, 10, 1000]; % Initial guess for [beta rho1 rho2 tauC1 and tauC2]
    lb = [1, 0, 0, 10, 5]; % lower boundary for [beta rho1 rho2 tauC1 and tauC2]
    ub = [1, 1, 1, 10, 10e3]; % upper boundary for [beta rho1 rho2 tauC1 and tauC2]
    options = optimoptions('lsqcurvefit','FunctionTolerance',1e-6,'Display','off');
    [varFit, resnorm, residual, exitflag, output] = lsqcurvefit(objfun, x0, T', K, lb, ub,options);
    fittingresults.varFit = varFit;
    fittingresults.KFit = objfun(varFit, T');
    fittingresults.R = (1-sum(abs(K-fittingresults.KFit).^2)./sum(abs(K-mean(K)).^2));
end