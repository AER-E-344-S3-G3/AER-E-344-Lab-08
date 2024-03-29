clear, clc, close all
%% Kyle Ostendorf Lab 8
%% Constants
x = [0:10,15:5:65];
n_x = length(x); % number of tests
taps = 1:38;
n_taps = length(taps); % number of taps
bad_taps = [1,2,21,22,35];  
good_taps = taps;
good_taps(bad_taps) = [];
d_m = .004; % [m]
y = taps * d_m;

mu = 1.8e-5; % [N*s/m^2]
rho = 1.225; % [kg/m^3]
IN0 = readmatrix("0in.csv");
chosen = 2:2:20; %% Data chosen to be graphed
%% Reading Cell array
total_p = zeros(n_x, n_taps);
static_p = zeros(n_x, 1);
for i = 1:length(x)

MAREIX = readmatrix(sprintf("%din.csv",x(i)));
total_p(i,:) = mean([MAREIX(:,2:17),MAREIX(:,35:50),MAREIX(:,68:73)],"omitmissing");
static_p(i,:) = mean(MAREIX(:,74));
end

%% Interpolating Data
interp_data = zeros(n_x, length(bad_taps));
i_total_p = zeros(n_x, n_taps);
for i = 1:n_x
good_data = total_p(i,good_taps);
interp_data(i,:) = interp1(good_taps,good_data,bad_taps,"linear","extrap");
i_total_p(i,:) = [interp_data(i,1:2),good_data(1:18),interp_data(i,3:4),good_data(19:30),interp_data(i,5),good_data(31:33)];
end

%% Velocity
Q = i_total_p - static_p;
U = sqrt(2*Q/rho);
U_inf = zeros(1,n_x);
U_Uinf = zeros(n_x, n_taps);
for i = 1:n_x
    U_inf(i) = mean(U(i,25:1:30));
    for j = 1:n_taps
        U_Uinf(i,j) = U(i,j)/U_inf(i);
    end
end

%% Delta Estimations
delta = zeros(1,n_x);
y_delta = zeros(n_x,n_taps);
for i = 1:n_x
    delta_loc = 1;
    u = 0.99*U_inf(i);
    while (U(i,delta_loc) < u)
        delta_loc = delta_loc +1;
    end
    delta(i) = d_m*(delta_loc-1);
    y_delta(i,:) = y / delta(i);
end

%% Plot y_delta vs. U_Uinf
for i  = chosen
    figure(i)
    plot(y_delta(i,:), U_Uinf(i,:))
    hold on
    plot([y_delta(i,1),y_delta(i,n_taps)],[.99,.99])

    fontname("Times New Roman");
    fontsize(12, "points");
    title_str = "Y/\delta vs. U/U_{inf} at "+ x(i) +" in";
    title(title_str);
    xlabel("Y/\delta");
    ylabel("U/U_{inf}");
    grid on;
    legend({'Data', '99% of U_{inf}'},"Location",'southeast')
end

%% Displacement and Momentum thickness
theta = zeros(1,n_x);
delta_star = zeros(1,n_x);
for i = 1:n_x
    for j = 1:n_taps-1
        f1 = U(i,j)/U_inf(i)*(1- U(i,j)/U_inf(i));
        f2 = U(i,j+1)/U_inf(i)*(1-U(i,j+1)/U_inf(i));
        theta(i) = theta(i) + d_m*((f1 + f2)/2);

        f3 = (1 - U(i,j)/U_inf(i));
        f4 = (1 - U(i,j+1)/U_inf(i));
        delta_star(i) = delta_star(i) + d_m*((f3 + f4)/2);
    end
end

%% Theoretical vs Experimental delta
V_inf = mean(U_inf); % [m/s] <-- this is for 10 Hz
Re_transition = 10^5; % []
m_to_in = 39.3700787;
x_transition = Re_transition * mu / (rho * V_inf); % [m]
x_transition = x_transition * m_to_in; % [in]
x_laminar = 0 : 0.1 : x_transition; % [in]
x_turbulent = x_transition : 0.1 : 70; % [in]
Re_laminar = rho * V_inf * (x_laminar * m_to_in^-1) / mu; % []
Re_turbulent = rho * V_inf * (x_turbulent * m_to_in^-1) / mu; % []
boundary_layer_laminar = 5.0 * x_laminar ./ sqrt(Re_laminar); % [in]
boundary_layer_turbulent = 0.37 * x_turbulent ./ Re_turbulent.^(1 / 5); % [in]

% Output
figure(n_taps+1);
yyaxis left
plot(x_laminar, boundary_layer_laminar*m_to_in^-1);
hold on;
plot(x_turbulent, boundary_layer_turbulent*m_to_in^-1);
hold on
plot(x(chosen),delta(chosen));
hold on
ylabel("\delta [m]");
yyaxis right
plot(x(chosen), theta(chosen));
ylabel("\theta");
hold off
title("Boundary Layer Thickness vs. Distance from LE");
xlabel("x [in]");
legend({'Laminar','Turbulent','Experimental \delta','Experimental \theta'},"Location",'southeast')
grid on;

%% Shear Stress Coefficient
cf = gradient(theta,x);
cf_rel = [(0.0583 ./ Re_laminar.^0.2), (0.0583 ./ Re_turbulent.^0.2)];
% Probably the wrong way to do cf 
% tau_w = zeros(n_x,n_taps);
% cf = zeros(n_x,n_taps);
% for i = 1:n_x
%     tau_w(i,:) = mu * gradient(U(i,:),y);
%     cf(i,:) = tau_w(i,:) / (1/2 * rho * U_inf(i)^2);
% end
figure(n_taps +2)
plot(x,cf)
hold on
plot([x_laminar, x_turbulent], cf_rel)
title("Local Shear Stress Coefficient vs. Distance from LE");
xlabel("x [in]");
ylabel("c_{f}")
legend({'Experimental','Emperical Relation'},"Location",'best')
grid on;

%% Drag Coefficient
Cd_rel = [(0.074 ./ Re_laminar.^0.2), (0.074 ./ Re_turbulent.^0.2)];
Cd = 2 * theta ./ x;
figure(n_taps +3)
plot(x,Cd)
hold on
plot([x_laminar, x_turbulent], Cd_rel)
title("Drag Coefficient vs. Distance from LE");
xlabel("x [in]");
ylabel("C_{d}")
legend({'Experimental','Emperical Relation'},"Location",'best')
grid on;

% gg ez need to add svg save to each plot and fix figure numbering
% and probably some calculations cause stuff be wack