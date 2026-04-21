function [T_o2, T_h2] = ceainlet(dmdt, ofr)
    Ru = 8.314;
    tank.o2.p = 6e5;
    tank.o2.T = 90;
    tank.o2.V = 30;
    tank.o2.R = Ru/32;
    tank.o2.dmdt = dmdt*ofr/(ofr+1);
    tank.o2.rho = py.CoolProp.CoolProp.PropsSI('D','T',tank.o2.T,'P',tank.o2.p,'Oxygen');
    tank.o2.fluid = 'Oxygen';
    tank.h2.p = 4e5;
    tank.h2.T = 20;
    tank.h2.V = 30;
    tank.h2.R = Ru/2;
    tank.h2.dmdt = dmdt-tank.o2.dmdt;
    tank.h2.rho = py.CoolProp.CoolProp.PropsSI('D','T',tank.h2.T,'P',tank.h2.p,'H2');
    tank.h2.fluid = 'H2';

    h2comp = comp(tank.h2, 25, 0.95, 1.3);
    o2comp = comp(tank.o2, 80/6, 0.95, 1.3);

    h2heat = h2comp;
    h2heat.T = h2heat.T + 300; % based off SSME, higher to account for larger nozzle
    % Assuming no change in pressure
    
    w_req = (h2comp.w*tank.h2.dmdt + o2comp.w*tank.o2.dmdt)/tank.h2.dmdt;
    turbh2 = turb2(h2heat, w_req, 0.8); % need to tune this to get 80 bar output
    
    T_o2 = o2comp.T;
    T_h2 = turbh2.T;
end

function out = comp(in, pr, eta, Ar)
    out.p = in.p*pr;
    out.R = in.R;
    out.dmdt = in.dmdt;
    out.fluid = in.fluid;
    h1 = py.CoolProp.CoolProp.PropsSI('H','T',in.T,'P',in.p,in.fluid);
    h01 = h1 + in.V^2/2;
    s1 = py.CoolProp.CoolProp.PropsSI('S','T',in.T,'P',in.p,in.fluid);

    tmp = in.T; i = 0; err = 1;
    while err > 0.001 && i < 1e4
        i = i + 1;
        out.rho = py.CoolProp.CoolProp.PropsSI('D','T',tmp,'P',out.p,in.fluid);
        out.V = in.rho/out.rho*Ar*in.V;
        h2s = py.CoolProp.CoolProp.PropsSI('H','P',out.p,'S',s1,in.fluid);
        h02s = h2s + out.V^2/2;
        h02 = (h02s-h01)/eta+h01;
        h2 = h02 - out.V^2/2;
        out.T = py.CoolProp.CoolProp.PropsSI('T','H',h2,'P',out.p,in.fluid);
        err = abs((out.T-tmp)/out.T);
        tmp = out.T;
    end
    
    cp = py.CoolProp.CoolProp.PropsSI('CPMASS','T',out.T,'P',out.p,in.fluid);
    cv = py.CoolProp.CoolProp.PropsSI('CVMASS','T',out.T,'P',out.p,in.fluid);
    out.gamma = cp/cv;
    out.h = h2;
    out.h0 = h02;
    out.M = out.V/sqrt(out.gamma*in.R*out.T);
    out.p0 = out.p*(1+(out.gamma-1)/2*out.M^2)^((out.gamma-1)/(out.gamma));
    out.T0 = out.T*(1+(out.gamma-1)/2*out.M^2);
    out.w = out.h-h1+0.5*(out.V^2-in.V^2);
end

function out = turb2(in, w_req, pr)
    out.dmdt = in.dmdt;
    out.R = in.R;
    out.p = in.p*pr;
    cp = py.CoolProp.CoolProp.PropsSI('CPMASS','T',in.T,'P',out.p,in.fluid);
    out.T = in.T - w_req/cp;
end
