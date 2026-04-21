function F_D = drag(p, T, rho, A_f, u)
    if rho > 1e-3 && p > 0
        c = py.CoolProp.CoolProp.PropsSI('A','P',p,'T',T,'Air'); % Ambient speed of sound
        M = u/c;
        if M<=4
            C_D = interp1(0:0.25:4,...
                [0.06, 0.065, 0.07, 0.08, 0.11, 0.15, 0.155, 0.14, 0.13, 0.12, 0.115, 0.11, 0.105, 0.1, 0.097, 0.093, 0.09],...
                M); % From C_D vs M curve in M5
        else
            C_D = .2./M+.04; % Approximation of CD after M=4
        end
        F_D = 0.5*rho*u^2*A_f*C_D;
    else
        F_D = 0; % No drag if out of atmosphere
    end
end