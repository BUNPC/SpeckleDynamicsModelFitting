%%%%%%%% pixel wise fitting temporal or spatial contrast K %%%%%%%%%
% Input: 
%     T             - multiple exposure times  
%     K             - spatial or temporal K(T) as 1D vector [number of exposure times, 1]
%     functiontype  - 't': temperal; 's': spatial
% Output: 
%     fittingresults - fitting results saved as a structure with 3 fields
%         fittingresults.varFit: [beta, rho, tauC]
%         fittingresults.g2Fit: fitted g2(tau) as 1D vector [number of exposure times, 1]  
%         fittingresults.R: R^2 score fitting accuracy 
% Example: fittingresults = pixelfitK(T, K, 't')
%
% Note:1) this function is based on g1 model exp(-(tau/tauC)^n) with n=1; 
%         For n=0.5 or 2, one can directly modify on 'fun', for example: case 's' with n=0.5, 
%         fun = @(beta, rho, tauC, T) (integral(@(x) 2*beta./T.*(1-x./T).*(rho*exp(-(x./tauC)^0.5)+1-rho).^2, 0, T, 'ArrayValued', true)).^0.5;
%      2) To fix any parameter in [beta rho tauC] for fitting, 
%         one can set the upper and lower boundary the same, 
%         for example: fix beta=0.5, lb = [0.5 0 10e-3], ub = [0.5 1 20].
%
% last modified by Bingxue Liu, 08/28/2023

function fittingresults = pixelfitK(T, K, functiontype)
    switch functiontype
        case 't'
            fun = @(beta, rho, tauC, T) (integral(@(x) 2*beta./T.*(1-x./T).*((rho*exp(-x./tauC)+1-rho).^2-(1-rho).^2), 0, T, 'ArrayValued', true)).^0.5;
        case 's'
            fun = @(beta, rho, tauC, T) (integral(@(x) 2*beta./T.*(1-x./T).*(rho*exp(-x./tauC)+1-rho).^2, 0, T, 'ArrayValued', true)).^0.5;
    end
    funs = @(beta, rho, tauC, T) arrayfun(fun, beta, rho, tauC, T);

    % Define the objective function for lsqcurvefit
    objfun = @(varFit,T) funs(varFit(1)*ones(length(T),1),varFit(2)*ones(length(T),1),varFit(3)*ones(length(T),1),T);

    % Fit the function using lsqcurvefit
    x0 = [0.5, 0.5, 1]; % Initial guess for [beta rho tauC]
    lb = [1, 0, 10e-3]; % lower boundary for [beta rho tauC]
    ub = [1, 1, 20]; % upper boundary for [beta rho tauC]
    options = optimoptions('lsqcurvefit','FunctionTolerance',1e-6,'Display','off');
    [varFit, resnorm, residual, exitflag, output] = lsqcurvefit(objfun, x0, T', K, lb, ub,options);
    fittingresults.varFit = varFit;
    fittingresults.KFit = objfun(varFit, T');
    fittingresults.R = (1-sum(abs(K-fittingresults.KFit).^2)./sum(abs(K-mean(K)).^2));
end