%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORWARD STARTING OPTION PRICING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The purpose of the project is to check if the analytical solution
% of the forward starting option is equal to the numerical approximation.

% The code proceeds as follow: 
% 1) Analytical and numerical approximation of GBM
% 2) Pricing a European call option numerically and in closed form
% 3) Pricing a Forward starting option numerically and in closed form

%% Housekeeping
clear;
clc;
rng("default")
addpath("functions")

%% Parameters
k = 0.2;
S0 = 10;
sigma = 0.8;
T = 1;
tf = 0.5;
r = 0.01;
sim = 1000;
discr_steps = 500;
dt = T/discr_steps;

%% 1) Analytical and numerical approximation of GBM

St_exact = zeros(2, sim) + S0;
for i=2:discr_steps
    St_exact(i,:) = GBM_approx(St_exact(i-1, :), r, sigma, dt, sim) ;
end

% Check if theoretical moments coincides with empirical one

%E[]
exp_value_theor = S0*exp(r*T);
exp_value_exact = mean(St_exact(end,:));

%Var()
std_theor = sqrt(exp(2*r*T)*(S0^2)*(exp((sigma^2)*T) - 1));
std_exact = std(St_exact(end,:));

figure;
plot(St_exact(:, :))
xline(tf*discr_steps)
ylabel('S(t)')
saveas(gcf,'images/MC_sim_GBM.jpg', 'jpg')

%% 2) Pricing a European call option numerically and in closed form

ST_exact = St_exact(end, :);
EC_price_numerical_exact = exp(-r*T)*mean(max(ST_exact-k*S0, 0));

d1 = (log(1/k) + (r + (sigma^2)/2)*T)/(sigma*sqrt(T));
d2 = (log(1/k) + (r - (sigma^2)/2)*T)/(sigma*sqrt(T));
EC_price_anal = S0*normcdf(d1) - exp(-r*T)*k*S0*normcdf(d2);

sim_iter = round(linspace(100, 10000, 200));
EC_price_numerical_exact_5 = zeros(length(sim_iter), 1);
EC_price_numerical_exact = zeros(length(sim_iter), 1);
EC_price_numerical_exact_95 = zeros(length(sim_iter), 1);
for i=1:length(sim_iter)
    i
    sim = sim_iter(i);
    St_exact = zeros(2, sim) + S0;
    for j=2:discr_steps
        St_exact(j,:) = GBM_approx(St_exact(j-1, :), r, sigma, dt, sim) ;
    end
    ST_exact = St_exact(end, :);
    tmp = max(ST_exact-k*S0, 0);
    EC_price_numerical_exact_5(i, 1) = exp(-r*T)*(mean(tmp) - 1.96*std(tmp)/sqrt(sim));
    EC_price_numerical_exact(i, 1) = exp(-r*T)*mean(tmp);
    EC_price_numerical_exact_95(i, 1) = exp(-r*T)*(mean(tmp) + 1.96*std(tmp)/sqrt(sim));
end

figure;
plot(sim_iter, EC_price_numerical_exact, ...
    sim_iter, EC_price_numerical_exact_5, "--r", ...
    sim_iter, EC_price_numerical_exact_95, "--r")
yline(EC_price_anal, '-')
ylabel('Call price')
xlabel('Number of simulations')
legend("Expected value", "95th CI")
saveas(gcf,'images/EC_approx.jpg', 'jpg')


%% 3) Pricing a Forward starting option numerically and in closed form

ST_exact = St_exact(end, :);
Stf_exact = St_exact(tf*T*discr_steps, :);
FSO_price_numerical_exact = exp(-r*T)*mean(max(ST_exact-k*Stf_exact, 0));

d1 = (log(1/k) + (r + (sigma^2)/2)*(T-tf))/(sigma*sqrt(T-tf));
d2 = (log(1/k) + (r - (sigma^2)/2)*(T-tf))/(sigma*sqrt(T-tf));
FSO_price_anal = S0*normcdf(d1) - exp(-r*(T-tf))*k*S0*normcdf(d2);

sim_iter = round(linspace(100, 10000, 200));
FSO_price_numerical_exact_5 = zeros(length(sim_iter), 1);
FSO_price_numerical_exact = zeros(length(sim_iter), 1);
FSO_price_numerical_exact_95 = zeros(length(sim_iter), 1);
for i=1:length(sim_iter)
    i
    sim = sim_iter(i);
    St_exact = zeros(2, sim) + S0;
    for j=2:discr_steps
        St_exact(j,:) = GBM_approx(St_exact(j-1, :), r, sigma, dt, sim) ;
    end
    ST_exact = St_exact(end, :);
    Stf_exact = St_exact(tf*T*discr_steps, :);
    tmp = max(ST_exact-k*Stf_exact, 0);
    FSO_price_numerical_exact_5(i, 1) = exp(-r*T)*(mean(tmp) - 1.96*std(tmp)/sqrt(sim));
    FSO_price_numerical_exact(i, 1) = exp(-r*T)*mean(tmp);
    FSO_price_numerical_exact_95(i, 1) = exp(-r*T)*(mean(tmp) + 1.96*std(tmp)/sqrt(sim));
end

figure;
plot(sim_iter, FSO_price_numerical_exact, ...
    sim_iter, FSO_price_numerical_exact_5, "--r", ...
    sim_iter, FSO_price_numerical_exact_95, "--r")
yline(FSO_price_anal, '-')
ylabel('Call price')
xlabel('Number of simulations')
legend("Expected value", "95th CI")
saveas(gcf,'images/FS_approx.jpg', 'jpg')



