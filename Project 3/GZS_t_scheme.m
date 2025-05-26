function g_new = GZS_t_scheme(x1,y1,x2,y2,R,x_circ,y_circ,T_f,T_ff,w,Rho,g,g_eq,gg,gg_eq,Tau)

    % Can be sped up by precomputing w*Rho
    T_w = 1;

    c_w = find_the_wall_point(x1,y1,x2,y2, R, x_circ, y_circ);
    delta = sqrt((c_w(1) - x2)^2 + (c_w(2) - y2)^2)/sqrt((x1-x2)^2 + (y1-y2)^2);

    T_b1 = ((delta-1)*T_f + T_w)/delta;

    if delta >= 0.75
        T_b = T_b1;
        g_a_eq = Rho*R*T_b*w;
        g_a_neq = g - g_eq;
    elseif delta < 0.75
        T_b2 = ((delta-1)*T_ff + 2*T_w)/(1 + delta);
        T_b = delta*T_b1 + (1-delta)*T_b2;
        g_a_eq = Rho*R*T_b*w;
        g_a_neq = g - g_eq + (1-delta)*(gg - gg_eq);
    end

    g_new = g_a_eq + (1 - 1/Tau)*g_a_neq;
end

