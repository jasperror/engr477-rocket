function [m, dmdt, h, u_eq, F_D, du, a, u, E, t_b, t_ea, t_end, der, in, e, L, oxfl, T] = solver(...
    m_rocket,... Rocket empty mass, kg
    m_payload,... Payload mass, kg
    m_propellant,... Rocket total propellant mass, kg
    h_orbit,... Orbital altitude, m
    ceain,... CEA input data
    A_e,... Nozzle velocity area, m^2
    A_f,... Rocket frontal area, m^2
    flag... If you only want the initial cea data
    )
    
    m_zp = m_rocket + m_payload; % Zero propellent mass, kg
    m = m_zp + max([m_propellant,50]); % Total initial mass, kg
    dt = 0.1; % Timestep, s
    t_CEA = 100; % How often to run the CEA, s (VERY PERFORMANCE INTENSIVE!)
    h = 0;
    u = 0;
    a = 0;
    F_D = 0;
    e = 5*9.81;
    in = e;
    dmdt = 400;
    i = 1;
    Kp = 3;
    % Kp = 1;
    Kd = 1.5;
    % Kd = 0;
    Ki = 2.5;
    % Ki = 0;
    Kb = 0.1;
    r_e = 6371e3; % Earth radius, km
    m_e = 5.9722e24; % Earth mass, kg
    G = 6.6743e-11; % Gravitational constant, m^3/kg-s^2
    v_orbit = sqrt(G*m_e/(r_e+h_orbit));
    E_req = 1.05*(m_payload*ambient(h_orbit)*h_orbit+0.5*m_payload*v_orbit^2);
    turnrate = 10*dt; % 10 deg/s turn (ignoring inertia and probably a lot of orbital mechanics)
    initialcea = CEA('problem','rocket','eql','fac','acat',ceain.acat,'phi',ceain.phi,'p,bar',ceain.p,'supsonic(ae/at)',ceain.aeat,'reactants','fuel','H2','wt%',100,'t,k',300,'oxid','O2','O',2,'wt%',100,'t,k',300,'output','transport','short','mks','end');
    oxfl = initialcea.output.oxfl;

    while m(i)>m_zp && i < 1e5
        i = i + 1;
        if i == 2
            [T_o2, T_h2] = ceainlet(dmdt(i-1), oxfl);
            warning('off','all');
            cearesult = CEA('problem','rocket','froz','nfz',1,'fac','acat',ceain.acat,'phi',ceain.phi,'p,bar',ceain.p,'supsonic(ae/at)',ceain.aeat,'reactants','fuel','H2','wt%',100,'t,k',T_h2,'oxid','O2','O',2,'wt%',100,'t,k',T_o2,'output','transport','short','mks','end');
            warning('on','all');
            u_e = cearesult.output.froz.mach(4)*cearesult.output.froz.sonvel(4)*cosd(ceain.alpha); % Exit velocity, m/s
            p_e = cearesult.output.froz.pressure(4)*1e5; % Exit pressure (bar->Pa), Pa
            L = combustorsize(cearesult);
            % u_e = cearesult.output.eql.mach(4)*cearesult.output.eql.sonvel(4);
            % p_e = cearesult.output.eql.pressure(4)*1e5;
            if flag
                u_eq = NaN;
                F_D = NaN;
                du = NaN;
                E = NaN;
                t_b = NaN;
                t_ea = NaN;
                t_end = NaN;
                der = NaN;
                L = NaN;
                return
            end
        elseif mod(i,round(t_CEA/dt)) == 0
            [T_o2, T_h2] = ceainlet(dmdt(i-1), oxfl);
            warning('off','all');
            cearesult = CEA('problem','rocket','froz','nfz',1,'fac','acat',ceain.acat,'phi',ceain.phi,'p,bar',ceain.p,'supsonic(ae/at)',ceain.aeat,'reactants','fuel','H2','wt%',100,'t,k',T_h2,'oxid','O2','O',2,'wt%',100,'t,k',T_o2,'output','transport','short','mks','end');
            warning('on','all');
            u_e = cearesult.output.froz.mach(4)*cearesult.output.froz.sonvel(4)*cosd(ceain.alpha); % Exit velocity, m/s
            p_e = cearesult.output.froz.pressure(4)*1e5; % Exit pressure (bar->Pa), Pa
            % u_e = cearesult.output.eql.mach(4)*cearesult.output.eql.sonvel(4);
            % p_e = cearesult.output.eql.pressure(4)*1e5;
        end
        
        dmts = -dmdt(i-1)*dt; % Propellant mass flow rate, kg/timestep
        m(i) = m(i-1)+dmts; % Total rocket mass, kg
        if h(i-1) < h_orbit
            h(i) = h(i-1)+u(i-1)*dt; % Altitude, m
            theta(i) = orbit();
        else % If we have reached the required altitude, just accelerate
            theta(i) = min([90,theta(i-1)+turnrate]);
            h(i) = h(i-1)+u(i-1)*dt*cosd(theta(i))*dt;
        end

        [g, T_amb, p_amb, rho_amb] = ambient(h(i));
        u_eq(i) = u_e + A_e/dmdt(i-1)*(p_e-p_amb);
        % u_eq(i) = u_e + A_e/dmdt(i-1)*(p_e-p_amb);
        
        F_D(i) = drag(p_amb, T_amb, rho_amb, A_f, u(i-1));
        du(i) = -u_eq(i)*dmts/m(i) - F_D(i)/m(i)*dt - g*cosd(theta(i))*dt;
        if h(i) == 0 && du(i) < 0 % if we are stuck on the ground
            du(i) = 0; % let the controller add more POWER!!
        end
        u(i) = u(i-1)+du(i);
        a(i) = du(i)/dt;
        E(i) = m_payload*g*h(i) + 0.5*m_payload*u(i)^2;
        T(i) = du(i)*dt/m(i) + F_D(i) + m(i).*g.*cosd(theta(i));

        % PID of dm/dt
        e(i) = 9.81*5-a(i); % Error
        in(i) = in(i-1)+e(i)*dt; % Integral error
        try
            ebar(i) = mean([e(i),e(i-1),e(i-2)]);
            der(i) = mean([e(i)-e(i-1),e(i-1)-e(i-2)])/dt; % Derivative error
            der(i) = mean(der(i-5:i));
        catch exception
            der(i) = (e(i)-e(i-1))/dt;
            ebar(i) = e(i); %#ok<*AGROW>
        end
        dmdtraw = Kp*ebar(i) + Ki*in(i) + Kd*der(i); % Propellant mass flow rate (PID), kg/s
        dmdt(i) = dmdtraw;
        dmdt(i) = max([50, min([dmdtraw,3000])]);
        dmdt(i) = max([0.98*dmdt(i-1),min([1.02*dmdt(i-1),dmdt(i)])]);
        if dmdt(i) ~= dmdtraw
            in(i) = in(i-1)+Ki*e(i)*dt + Kb*(dmdt(i)-dmdtraw)*dt;
            % in(i) = in(i-1);
        end
        if m-dmdt(i)*dt < m_zp
            dmdt(i) = (m-m_zp)/dt;
        end
        % dmdt(i) = 200;
    end

    t_b = i*dt; % burnout time, s
    while F_D(i) > 0 % In case burn didnt reach edge of atmosphere, keep checking the du
        i = i+1;
        u_eq(i) = 0;
        dmdt(i) = 0;
        T(i) = 0;
        m(i) = m(i-1);
        if h(i-1) < h_orbit
            h(i) = h(i-1)+u(i-1)*dt; % Altitude, m
            theta(i) = orbit();
        else % If we have reached the required altitude, just accelerate
            theta(i) = min([90,theta(i-1)+turnrate]);
            h(i) = h(i-1)+u(i-1)*dt*cosd(theta(i))*dt;
        end
        [g, T_amb, p_amb, rho_amb] = ambient(h(i));
        F_D(i) = drag(p_amb, T_amb, rho_amb, A_f, u(i-1));
        du(i) = - F_D(i)/m(i)*dt - g*cosd(theta(i))*dt;
        u(i) = u(i-1)+du(i);
        a(i) = du(i)/dt;
        E(i) = m_payload*g*h(i) + 0.5*m_payload*u(i)^2;
    end
    t_ea = i*dt; % Exit atmosphere time, s
    j = i;
    m(i) = m_payload; % Ditch the rocket now that we are out of the atmosphere
    u(i) = u(i)+2; % Assuming the separation adds 2 m/s to the velocity for clearance.
    dt = 10*dt; % Take bigger steps
    while (h(i) < h_orbit || E(i) < E_req) && h(i)-h(i-1) > 0
        i = i + 1;
        u_eq(i) = 0;
        F_D(i) = 0;
        T(i) = 0;
        dmdt(i) = 0;
        m(i) = m(i-1);
        if h(i-1) < h_orbit
            h(i) = h(i-1)+u(i-1)*dt; % Altitude, m
            theta(i) = orbit();
        else % If we have reached the required altitude, just accelerate
            theta(i) = min([90,theta(i-1)+turnrate]);
            h(i) = h(i-1)+u(i-1)*dt*cosd(theta(i))*dt;
        end
        du(i) = - g*cosd(theta(i))*dt;
        u(i) = u(i-1)+du(i);
        a(i) = du(i)/dt;
        E(i) = m_payload*g*h(i) + 0.5*m_payload*u(i)^2;
    end
    t_end = t_ea+(i-j)*dt; % Reached (or failed to reach) orbit time, s
end
