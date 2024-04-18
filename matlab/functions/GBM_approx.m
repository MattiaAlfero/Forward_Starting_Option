function St1_exact = GBM_approx(St_exact, mu, sigma, dt, sim) 
    zt = normrnd(0, 1, 1, sim);
    St1_exact = St_exact.*exp((mu - (sigma^2)/2)*dt + sigma*sqrt(dt).*zt);
end

        