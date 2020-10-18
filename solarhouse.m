% todo:

% fix bug where thicker walls makes house colder??


% How does the storage unit thickness affect the house’s thermal behavior?
% 
% The thickness of the storage unit affects the temperature range
% experienced by the floor and air in the house, but does not influence the
% average temperature over time. The thinner the absorber, the more
% fluctuation occurs on a daily basis, and the longer it would take for the
% house to reach a comfortable average temperature over time.
% 
% How does the insulation thickness affect the house’s thermal behavior?
% 
% The thickness of the insulation controls the asymptotic behavior of the
% system, and in turn the resultant temperature in the house, as well as the
% time it takes to get there. The thinner the insulation, the less time it
% takes to reach a steady state, and the lower the steady state temperature
% will be. Reversing this, with thicker walls, there's more insulation and
% therefore a higher resultant internal temperature, but it will take
% longer to achieve this state.
% 
% Why?
% 
% Once again, we can compare the the absorber to a capacitor and the
% insulation to a resistor, which means we've modeled something similar to
% an RC circuit. In varying the thicknesses of each of these
% parts, we're effectively changing the resultant air temperature
% (analogous to voltage, affected by resistance) and the smoothing on the
% air temperature (analogous to the cutoff frequency of an electronic
% filter, which affects voltage over time). 

clc, clear; % clean slate


% Plot temperature over time with optimized values
[t, dT, M] = housetemps(0.5, .1, 6, 2, 3, 2, 40);
fig1 = figure(1);
hold on;
grid on;
plot(t,dT,'-')
plot(t,M, '--')
h = legend("$T_{floor}$", "$T_{air}$", "$\overline{T}_{floor}$", "$\overline{T}_{air}$");
set(h,'interpreter','Latex','FontSize',12)
title("Temperatures in House & 4 hr Moving Average");
xlabel("Time (seconds)");
ylabel("Temperature (Celsius)");
hold off

% Plot temperature over time with optimized values over a single day
fig2 = figure(2);
hold on;
grid on;
plot(t,dT,'-')
plot(t,M, '-.')
h = legend("$T_{floor}$", "$T_{air}$", "$\overline{T}_{floor}$", "$\overline{T}_{air}$");
set(h,'interpreter','Latex','FontSize',12);
title("Temperatures in House & 4 hr Moving Average, Single Day");
xlabel("Time (seconds)");
ylabel("Temperature (Celsius)");
x_start = 993000;
axis([x_start x_start+86400 0 40])
hold off;

% make fancy 3d plot
figure(3);
hold on;
grid on;
absorb_thick_range = linspace(0.1, 1, 50); % thickness of thermal mass/absorber
insul_thick_range = linspace(0.2, 1.5, 50); % thickness of wall insulation
heats = zeros(50, 50);
for a = 1:50
    for i = 1:50
        [t, dT, M] = housetemps(absorb_thick_range(a), insul_thick_range(i), 6, 2, 3, 2, 30);
        heats(a, i) = M(end);
    end
    a
end
contour(absorb_thick_range, insul_thick_range, heats, 'ShowText', 'On')
title('Sweeping absorber and insulator thickness values')
xlabel('Absorber thickness (m)'), ylabel('Wall insulator thickness (m)')
legend('Degrees celsius')
hold off;

function [t, dT, M] = housetemps(absorber_thickness, insulation_thickness, ...
                              h_length, h_width, h_height, g_height, ...
                              num_days)

% standard values
% absorber_thickness = .5; % m
% insulation_thickness = .2; % m
% h_length = 6; % m
% h_width = 2; % m
% h_height = 3; % m
% g_height = 2; % m
% num_days = 40; % days

tspan = [0 86400*num_days]; % s

g_area = h_length * g_height; % m ^ 2
% g_thick = 0.005; % m

r_height = sqrt(h_width^2 + (h_height - g_height)^2); % m
r_area = h_length * r_height; % m ^ 2
% r_thick = 0.1; % m

bwa_area = h_length * h_height; % m ^ 2
swa_area = (g_height * h_width) + .5 * (h_width * (h_height - g_height)); % m ^ 2
wa_thick = insulation_thickness; % m

f_area = h_width * h_length; % m ^ 2
f_thick = absorber_thickness; % m

ins_net_area = bwa_area + 2 * swa_area + r_area + f_area; % m ^ 2

% material thermal properties
h_glass = .7; % W/m^2-K
h_in = 15; % W/m^2-K
h_out = 30; % W/m^2-K

k_wall = .4; % W/m-K

v_air = swa_area * h_width; % m^3
d_air = 1.225; %kg/m^3
m_air = d_air * v_air; %kg
c_air = 1012; %J/kg-K

v_abs = f_area * f_thick; % m^3
d_abs = 3000; % kg/m^3
m_abs = d_abs * v_abs; % kg

C_abs = 800 * m_abs; % J/K
C_a = c_air * m_air; % J/K

% partial resistances
R_abs_air = 1/(h_in * f_area);

R_air_wall = 1/(h_in * ins_net_area);
R_wall_wall = wa_thick/(k_wall * ins_net_area);

R_air_glass = 1/(h_in * g_area);
R_glass_glass = 1/(h_glass * g_area);

R_outer_air = 1/(h_out * ins_net_area);

% equivalent resistances
R_fa = R_abs_air; % total thermal resistance from floor to air
R_ao = + 1/(1/(R_air_wall + R_wall_wall)+1/(R_air_glass + R_glass_glass)) + R_outer_air; % total thermal resistance from air to outside

% f = @(t,T) [(1/C_abs)*(g_area*(-361*cos(pi*t/43200) + 224*cos(pi*t/21600) + 2210) - ((T(1)-T(2))/R_fa)); ...
%             (1/C_a)*((T(1)-T(2))/R_fa) - (T(2)-(6*sin((2*pi*t/86400)+3*pi/4)-3)/R_ao)];

f = @(t,T) [(1/C_abs)*(g_area*(-361*cos(pi*t/43200) + 224*cos(pi*t/21600) + 210) - ((T(1)-T(2))/R_fa)); ...
            (1/C_a)*(((T(1)-T(2))/R_fa) - (T(2) - -3)/R_ao)];

[t,dT] = ode45(f, tspan, [0 0]); % compute the ode

M = movmean(dT,3600*[4 0]); % 4 hour moving average

end