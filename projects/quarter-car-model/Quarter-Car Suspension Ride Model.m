%% quarter_car.m
% 2-DOF quarter-car ride model (sprung + unsprung mass)
% States: x = [zs; vs; zu; vu]   (displacements in m, velocities in m/s)

clear; clc; close all;

%% 1.Parameters
p.ms = 70;      % sprung mass, quarter body [kg]
p.mu = 15;       % unsprung mass [kg]
p.ks = 18000;    % suspension stiffness at the wheel [N/m]
p.cs = 1200;     % damper coefficient at the wheel [N*s/m]
p.kt = 75000;   % tire vertical stiffness [N/m]

% p.ks = k_spring * MR^2;   p.cs = c_damper * MR^2;   (MR = motion ratio)

%% 2.natural frequencies and damping from eigenvalues
[A, ~] = state_matrices(p);
ev = eig(A);
ev = ev(imag(ev) > 0);                 % one of each complex pair
[~, idx] = sort(abs(ev));
ev = ev(idx);
fn   = abs(ev) / (2*pi);               % natural frequency [Hz]
zeta = -real(ev) ./ abs(ev);           % damping ratio [-]

keff = p.ks*p.kt/(p.ks + p.kt);
fn_body_est = sqrt(keff/p.ms)/(2*pi);
fn_hop_est  = sqrt((p.ks + p.kt)/p.mu)/(2*pi);

fprintf('--- Mode check ---\n');
fprintf('Body mode:  %.2f Hz, zeta = %.2f  (hand estimate %.2f Hz)\n', fn(1), zeta(1), fn_body_est);
fprintf('Wheel hop:  %.2f Hz, zeta = %.2f  (hand estimate %.2f Hz)\n', fn(2), zeta(2), fn_hop_est);
fprintf('Expect body ~2.5 Hz, wheel hop ~10-15 Hz.\n\n');

%% 3. Step bump response (5 cm step at t = 0.2 s)
h = 0.05;
road_step = @(t) h * (t >= 0.2);

[t, x, out] = simulate(p, road_step, [0 3]);

fprintf('--- Step response check ---\n');
fprintf('Final zs = %.4f m, zu = %.4f m (both should approach %.2f)\n', x(end,1), x(end,3), h);
fprintf('Final suspension travel = %.5f m (should approach 0)\n', out.travel(end));
fprintf('Peak body accel = %.2f m/s^2 (%.2f g)\n\n', max(abs(out.acc)), max(abs(out.acc))/9.81);

figure('Name','Step response');
subplot(4,1,1);
plot(t, road_step(t), 'k--', t, x(:,1), 'b', t, x(:,3), 'r', 'LineWidth', 1.2);
ylabel('Position [m]'); legend('road','body z_s','wheel z_u'); grid on;
title('5 cm step bump');
subplot(4,1,2);
plot(t, out.acc, 'LineWidth', 1.2);
ylabel('Body accel [m/s^2]'); grid on;
subplot(4,1,3);
plot(t, out.travel, 'LineWidth', 1.2);
ylabel('Susp. travel z_s - z_u [m]'); grid on;
subplot(4,1,4);
plot(t, out.tire_defl, 'LineWidth', 1.2);
ylabel('Tire defl. z_u - z_r [m]'); xlabel('Time [s]'); grid on;

%% 4. Sine road input (near body resonance)
A_sine = 0.02;                       % 2 cm amplitude
f_sine = fn(1);                      % excite at body natural frequency
road_sine = @(t) A_sine * sin(2*pi*f_sine*t);

[t2, x2, out2] = simulate(p, road_sine, [0 6]);

figure('Name','Sine road');
subplot(3,1,1);
plot(t2, road_sine(t2), 'k--', t2, x2(:,1), 'b', 'LineWidth', 1.2);
ylabel('Position [m]'); legend('road','body z_s'); grid on;
title(sprintf('Sine road, %.2f Hz, %.0f mm amplitude', f_sine, A_sine*1000));
subplot(3,1,2);
plot(t2, out2.acc, 'LineWidth', 1.2);
ylabel('Body accel [m/s^2]'); grid on;
subplot(3,1,3);
plot(t2, out2.tire_defl, 'LineWidth', 1.2);
ylabel('Tire defl. [m]'); xlabel('Time [s]'); grid on;

%% 5. Damping sweep: comfort vs road holding
% Comfort    -> RMS body acceleration (lower is better)
% Road hold  -> RMS dynamic tire deflection (lower is better)
cs_vals = linspace(200, 3000, 15);
rms_acc  = zeros(size(cs_vals));
rms_tire = zeros(size(cs_vals));

for i = 1:numel(cs_vals)
    q = p;
    q.cs = cs_vals(i);
    [ti, ~, oi] = simulate(q, road_step, [0 3]);
    rms_acc(i)  = sqrt(trapz(ti, oi.acc.^2)       / (ti(end)-ti(1)));
    rms_tire(i) = sqrt(trapz(ti, oi.tire_defl.^2) / (ti(end)-ti(1)));
end

figure('Name','Damping sweep');
subplot(1,2,1);
yyaxis left;  plot(cs_vals, rms_acc, 'o-');  ylabel('RMS body accel [m/s^2]');
yyaxis right; plot(cs_vals, rms_tire*1000, 's-'); ylabel('RMS tire defl. [mm]');
xlabel('c_s [N s/m]'); grid on; title('Effect of damping');
subplot(1,2,2);
plot(rms_acc, rms_tire*1000, 'o-'); hold on;
plot(rms_acc(1), rms_tire(1)*1000, 'g^', 'MarkerSize', 9, 'MarkerFaceColor', 'g');
plot(rms_acc(end), rms_tire(end)*1000, 'rv', 'MarkerSize', 9, 'MarkerFaceColor', 'r');
xlabel('RMS body accel [m/s^2]  (comfort)');
ylabel('RMS tire defl. [mm]  (road holding)');
legend('sweep','lowest c_s','highest c_s','Location','best');
grid on; title('Trade-off curve');

%% ===== Local functions =====

function [A, B] = state_matrices(p)
% x_dot = A*x + B*zr
A = [ 0,            1,            0,            0;
     -p.ks/p.ms,  -p.cs/p.ms,    p.ks/p.ms,    p.cs/p.ms;
      0,            0,            0,            1;
      p.ks/p.mu,   p.cs/p.mu,   -(p.ks+p.kt)/p.mu, -p.cs/p.mu];
B = [0; 0; 0; p.kt/p.mu];
end

function dx = qc_ode(t, x, p, road)
[A, B] = state_matrices(p);
dx = A*x + B*road(t);
end

function [t, x, out] = simulate(p, road, tspan)
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', 1e-3);
[t, x] = ode45(@(t,x) qc_ode(t, x, p, road), tspan, zeros(4,1), opts);

zr = arrayfun(road, t);
out.acc       = zeros(size(t));
for i = 1:numel(t)
    dxi = qc_ode(t(i), x(i,:).', p, road);
    out.acc(i) = dxi(2);                 % sprung mass acceleration
end
out.travel    = x(:,1) - x(:,3);         % suspension travel
out.tire_defl = x(:,3) - zr;             % tire deflection
end