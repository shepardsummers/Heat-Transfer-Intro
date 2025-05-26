function g_eq_d2q9 = eqm_t_d2q9(Rho,T)
    % This function computes thermal pdf based on T (scalar) and Rho (scalar)
    
    w = [0 1/6 1/6 1/6 1/6 1/12 1/12 1/12 1/12];

    R = 8.314;
    
    g_eq_d2q9 = Rho*R*T*w';
end