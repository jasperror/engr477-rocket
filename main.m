% ENGR477: Aircraft Propulsion
% Project 2: Rockets
% By Jasper Palmer
clc; clear; close all;

% Rocket Masses
m_rocket = 5000; % Rocket mass, kg
m_payload = 300; % Payload mass, kg
m_empty = m_rocket + m_payload;
m_propellant = 37607;

% Orbit Considerations
h_orbit = 550e3; % Orbit altitude, m
r_e = 6371e3; % Earth radius, km
m_e = 5.9722e24; % Earth mass, kg
G = 6.6743e-11; % Gravitational constant, m^3/kg-s^2
v_orbit = sqrt(G*m_e/(r_e+h_orbit)); % Orbital velocity (assuming circular orbit, satelite mass negligible), m/s
E_req = 1.05*(m_payload*ambient(h_orbit)*h_orbit+0.5*m_payload*v_orbit^2); % Required energy to enter orbit (assuming a 5% additive for turning), J

% Rocket Dimensions
A_e = pi; % Nozzle outlet area, m^2
A_f = 1.5*pi; % Rocket frontal area, m^2

% CEA
ceain = struct('acat', 3, 'phi', 2, 'p', 80, 'aeat', 50,'alpha',14);

% Calculation
[m, dmdt, h, u_eq, F_D, du, a, u, E, t_b, t_ea, t_end, der, in, e, L, oxfl, T] = solver(m_rocket, m_payload, m_propellant, h_orbit, ceain, A_e, A_f, false);

% Plot
%%
close all
titles={'$m$ [kg]','$\dot{m}_p$ [kg/s]','$a$ [m/s$^2$]','$h$ [m]','$u$ [m/s]', '$E_{sat}$ [J]','$T$ [N]','$F_D$ [N]','$u_{eq}$ [m/s]'};
data = cat(1, m, dmdt, a, h, u, E, T, F_D, u_eq);
labels={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'};
figure('Position',[0 0 720 874])
tiledlayout
hold on
for i = 1:length(titles)
    nexttile
    plot([0.1:0.1:t_ea-0.1, t_ea:t_end], data(i,:))
    xlim([0, t_end])
    ylabel(titles(i),'Interpreter','latex')
    title(labels(i),'Interpreter','latex')
    xlabel('$t$ [s]','Interpreter','latex')
    if isequal(titles(i),{'$a$ [m/s$^2$]'})
        yline(9.81*5,'r')
        ylim([0, 9.81*5*1.5])
    elseif isequal(titles(i),{'$E_{sat}$ [J]'})
        yline(E_req,'g')
    elseif isequal(titles(i),{'$h$ [m]'})
        yline(h_orbit,'g')
    elseif isequal(titles(i),{'$u$ [m/s]'})
        yline(v_orbit,'g')
    elseif isequal(titles(i), {'$m$ [kg]'})
        yline(m_empty)
    end
end
hold off

figure
hold on
plot(0.1:0.1:t_b,e)
plot([0.1:0.1:t_ea-0.1, t_ea:t_end],a)
yline(5*9.81)
plot(0.1:0.1:t_b,der)
ylim([-40,80])
yyaxis right
plot([0.1:0.1:t_ea-0.1, t_ea:t_end],dmdt)
plot(0.1:0.1:t_b,in)
xlabel('t [s]','interpreter','latex')
grid on
legend('Error','Plant','','Derivative','Signal','Integral')
% close(2)

%% Solution Space for propellant mass
m_propellant = linspace(200, 100e3, 50);
E = 0;
for i = 1:length(m_propellant)
    [~, ~, ~, ~, ~, ~, ~, ~, tmpE, ~, ~, ~, ~, ~, ~, ~, ~] = solver(m_rocket, m_payload, m_propellant(i), h_orbit, ceain, A_e, A_f, false);
    E(i) = tmpE(end); %#ok<*SAGROW>
end

figure
plot(m_propellant,E)
yline(E_req,'g')

%% Energy Equation

function E = solverE(m_rocket, m_payload, m_propellant, h_orbit, ceain, A_e, A_f)
    [~, ~, ~, ~, ~, ~, ~, ~, E, ~, ~, ~, ~, ~, ~, ~, ~] = solver(m_rocket, m_payload, m_propellant, h_orbit, ceain, A_e, A_f, false);
    E = E(end);
end

%% Solution to propellant mass

func2 = @(m_propellant) solverE(m_rocket, m_payload, m_propellant, h_orbit, ceain, A_e, A_f)-E_req;
m_propellant = fzero(func2, 50e3)

%% Solution space for phi
phi = 1:0.5:3;
rho_t = 3000; % 2195 Al density, kg/m^3
S_t = 560e6; % 2195 Al yield strength, Pa
n = 1.6; % Safety factor
p_O2 = 6e5; % O2 storage pressure
p_H2 = 4e5; % H2 storage pressure
T_O2 = 90; % O2 storage temperature
T_H2 = 20; % H2 storage temperature
Ru = 8.314;
m_rocket_nt = 2500; % Mass of rocket minus storage tanks

m_prop = zeros(size(phi));
oxflv = m_prop; m_r = m_prop;
eps = 1e-2;
for i = 1:length(phi)
    ceain.phi = phi(i);
    err = 1;
    j = 0;
    m_propellant = 50e3; % Guess at required weight of propellant
    m_O2 = m_propellant*3/4;
    m_H2 = m_propellant/4;
    while err > eps && j < 20
        j = j + 1;
    
        V_O2 = m_O2/py.CoolProp.CoolProp.PropsSI('D','T',T_O2,'P',p_O2,'Oxygen');
        V_H2 = m_H2/py.CoolProp.CoolProp.PropsSI('D','T',T_H2,'P',p_H2,'H2');
        r_O2 = (V_O2*3/4/pi)^(1/3); % Assuming spherical tanks
        r_H2 = (V_H2*3/4/pi)^(1/3);
        
        t_O2 = p_O2*r_O2/(2*S_t/n-p_O2/2); % Tank wall thickness, m 
        t_H2 = p_H2*r_H2/(2*S_t/n-p_H2/2);
        m_tO2 = 4/3*pi*rho_t*((r_O2+t_O2)^3-r_O2^3); % Oxygen tank mass, kg
        m_tH2 = 4/3*pi*rho_t*((r_H2+t_H2)^3-r_H2^3); % Hydrogen tank mass, kg
        m_rocket = m_rocket_nt + m_tO2 + m_tH2;
        
        func2 = @(m_propellant) solverE(m_rocket, m_payload, m_propellant, h_orbit, ceain, A_e, A_f)-E_req;
        m_propellant = fzero(func2, m_propellant);
    
        [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, oxfl] = solver(m_rocket, m_payload, m_propellant, h_orbit, ceain, A_e, A_f, true);
        tmpo = m_propellant*oxfl/(oxfl+1);
        tmph = m_propellant-tmpo;
        err = sqrt(mean( [(tmpo-m_O2)/tmpo,(tmph-m_H2)/tmph].^2 ));
        m_O2 = tmpo;
        m_H2 = tmph;
        disp(['i=',num2str(i),', j=',num2str(j),', m_p=',num2str(round(m_propellant,1))])
    end
    m_prop(i) = m_propellant;
    oxflv(i) = oxfl;
    m_r(i) = m_rocket;
end
figure
plot(phi,m_prop)
ylabel('Propellant mass [kg]')
yyaxis right
plot(phi,m_r)
ylabel('Rocket mass [kg]')
xlabel("\phi")
ceain.phi = 2; % Optimal value from above analysis

%% Solution space for Ac/A*

acat = 2:0.2:4;
m_prop = zeros(size(acat));
for i = 1:length(acat)
    i
    ceain.acat = acat(i);
    
    func2 = @(m_propellant) solverE(m_rocket, m_payload, m_propellant, h_orbit, ceain, A_e, A_f)-E_req;
    m_prop(i) = fzero(func2, 50e3);
end
figure
plot(acat, m_prop)

%% Nozzle dimensions

alpha = 1; % Initial nozzle angle, deg
i = 1;
dx = 0.05;
cearesult = CEA('problem','rocket','froz','nfz',1,'fac','acat',3,'phi',2,'p,bar',80,'supsonic(ae/at)',100,'reactants','fuel','H2','wt%',100,'t,k',300.0,'oxid','O2','O',2,'wt%',100,'t,k',90,'output','transport','short','mks','end');
Cstar = cearesult.output.eql.cstar(2);
p2 = cearesult.output.eql.pressure(2)*1e5;
gamma = cearesult.output.eql.gamma(2);
M2 = cearesult.output.eql.mach(2);
p02 = p2*(1+M2^2*(gamma-1)/2)^(gamma/(gamma-1));
A = Cstar*550/p02; % Throat area, m^2
r = sqrt(A/pi);
x = 0;
Isp = 0;
m_nozzle = 0;

while i < 1e4
    i = i+1;
    if mod(i,10) == 0
        disp( ['i=', num2str(i), ', α=', num2str(alpha(i-1))] )
    end
    x(i) = dx*i;
    r(i) = r(i-1) + tand(alpha(i-1))*dx;
    A(i) = pi*r(i)^2;
    AeAt = A(i)/A(1);
    cearesult = CEA('problem','rocket','eql','fac','acat',3,'phi',2.5,'p,bar',80,'supsonic(ae/at)',AeAt,'reactants','fuel','H2','wt%',100,'t,k',300.0,'oxid','O2','O',2,'wt%',100,'t,k',200,'output','transport','short','mks','end');
    alpha(i) = atand(1/sqrt( cearesult.output.eql.mach(4)^2-1 )); % Nozzle angle
    Isp(i) = cearesult.output.eql.isp(4); % Isp
    p(i) = cearesult.output.eql.pressure(4)*1e5; % Pressure, Pa
    m_nozzle(i) = 2000*pi*( (r(i)+0.015)^2-r(i)^2 )*dx + m_nozzle(i-1);
    if abs(alpha(i)-alpha(i-1))/alpha(i) < 1e-3
        break
    end
end

figure('Position',[0 0 720 874])
tiledlayout(5,1)
data = [alpha; A; r; Isp; p];
titles = {'\alpha [deg]','A [m^2]','r [m]','I_{sp} [s]','p [Pa]'};
for i = 1:size(data,1)
    nexttile
    plot(x,data(i,:))
    xlim([0,max(x)])
    ylabel(titles(i))
    if isequal(titles(i), {'r [m]'})
        hold on
        plot(x,-data(i,:))
        yline(0,'k--')
        hold off
    end
end

%% Solution space for Ae/A*
% Need to run nozzle first

aeat = 50:10:120;
m_rocket_nn = 5e3; % mass of rocket

m_prop = zeros(size(aeat));
for i = 4:length(aeat)
    i
    if aeat(i) == 70
        m_prop(i) = NaN;
        continue % it does not like 70
    end
    j = find(A./A(1)>=aeat(i),1);
    m_rocket = m_rocket_nn + m_nozzle(j);

    ceain.aeat = aeat(i);
    ceain.alpha = interp1(A./A(1),alpha,aeat(i));
    func2 = @(m_propellant) solverE(m_rocket, m_payload, m_propellant, h_orbit, ceain, aeat(i)*A(1), max([A_f,aeat(i)*A(1)]))-E_req;
    m_prop(i) = fzero(func2, 50e3);
end
figure
plot(aeat, m_prop)
xlabel('Ae/A*')
ylabel('m_{propellant}')