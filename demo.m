clear all
% system settings 
tau = 1:100; % time lag [ms]
Texp = 1:5:1000; % exposure time [ms]
beta = 1;
rho = 0.5;
tauC = 10; % ms 

% generate g2 & K data through analytical expression
fun_g2t = @(beta, rho, tauC) 1+beta*(rho^2*exp(-tau/tauC).^2+2*rho*(1-rho)*exp(-tau/tauC));
dg2t = fun_g2t(beta, rho, tauC);

fun_g2s = @(beta, rho, tauC) 1+beta*(rho*exp(-tau./tauC)+1-rho).^2;
dg2s = fun_g2s(beta, rho, tauC);

fun_Kt = @(beta, rho, tauC, T) ...
    (integral(@(x) 2*beta./T.*(1-x./T).*((rho*exp(-x./tauC)).^2+2*rho*(1-rho)*(exp(-x./tauC))), 0, T, 'ArrayValued', true)).^0.5;
for iT = 1: length(Texp)
    dKt(iT) = fun_Kt(beta, rho, tauC, Texp(iT));
end

fun_Ks = @(beta, rho, tauC, T) ...
    (integral(@(x) 2*beta./T.*(1-x./T).*(rho*exp(-x./tauC)+1-rho)^2, 0, T, 'ArrayValued', true)).^0.5;
for iT = 1: length(Texp)
    dKs(iT) = fun_Ks(beta, rho, tauC, Texp(iT));
end

% fitting 
fg2t = pixelfitg2(tau, dg2t', 't');
fg2s = pixelfitg2(tau, dg2s', 's');

fKt = pixelfitK(Texp, dKt', 't');
fKs = pixelfitK(Texp, dKs', 's');

% plot results
figure;
plot(tau, dg2t,'*b'); hold on; plot(tau, fg2t.g2Fit,'b')
hold on;
plot(tau, dg2s,'*r'); hold on; plot(tau, fg2s.g2Fit,'r')
xlabel('\tau [ms]');
ylabel('g_2(\tau)');
grid on;
legend({'g^t_2(\tau), analytical','g^t_2(\tau), fitted', 'g^s_2(\tau), analytical','g^s_2(\tau), fitted'})
set(gca,'FontSize',12);

figure;
plot(Texp, dKt,'*b'); hold on; plot(Texp, fKt.KFit,'b')
hold on;
plot(Texp, dKs,'*r'); hold on; plot(Texp, fKs.KFit,'r')
xlabel('T [ms]');
ylabel('K(T)');
grid on;
legend({'K^t(T), analytical','K^t(T), fitted', 'K^s(T), analytical','K^s(T), fitted'})
set(gca,'FontSize',12);