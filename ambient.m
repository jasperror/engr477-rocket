% Ambient air conditions as a function of altitude
function [g, T, p, rho] = ambient(h)
    if h<11e3
        T = 15.04-.00649*h+273.15;
        p = 101.29*(T/288.08)^5.256;
    elseif h < 25e3
        T = -55.46+273.15;
        p = 22.65*exp(1.73-0.000157*h);
    else
        T = -131.21+0.00299*h+273.15;
        p = 2.488*(T/216.6)^-11.388;
    end
    p = p*1e3; % kPa to Pa
    rho = 1.2.*exp(-2.9e-5.*h.^1.15);
    g0 = 9.81; % Gravitational acceleration at surface, m/s^2
    r_e = 6371e3; % Earth radius, km
    g = g0*r_e^2/(r_e+h)^2; % Gravitational acceleration at h, m/s^2
end