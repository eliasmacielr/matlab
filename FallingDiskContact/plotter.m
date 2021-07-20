clearvars

R = 0.5;
m = 5;
I_A = 1/2*m*R^2;
I_T = 1/4*m*R^2;
g = 9.81;

alpha = 0.1; % concat this to file name

exp_n = 43; % Experiment number

%% Load contact result

load(strcat('res-contact2-h0.1-alpha',num2str(alpha),'.mat'))
T = tf - t0;

t = t0:h:tf;

% Get coordinates
X = q(1,:);
Y = q(2,:);
theta = q(3,:);
phi = q(4,:);
psi = q(5,:);

% Compute velocities from coordinates
Xdot = [y0(6),diff(X)/h];
Ydot = [y0(7),diff(Y)/h];
thetadot = [y0(8),diff(theta)/h];
phidot = [y0(9),diff(phi)/h];
psidot = [y0(10),diff(psi)/h];

E = 1/2*m*(Xdot.^2 + Ydot.^2 + R^2*sin(theta).*thetadot.^2) + ...
    1/2*(I_A*(psidot - phidot.*sin(theta)).^2 + ...
        I_T*(thetadot.^2 + phidot.^2.*(cos(theta).^2))) + ...
    m*g*R*cos(theta);


%% Load reference result

Ref = load(strcat('res-ode15i-h0.01-alpha',num2str(alpha),'.mat'));

X_ref = transpose(Ref.q(:,1));
Y_ref = transpose(Ref.q(:,2));
theta_ref = transpose(Ref.q(:,3));
phi_ref = transpose(Ref.q(:,4));
psi_ref = transpose(Ref.q(:,5));

h_ref = Ref.h;
t_ref = t0:h_ref:(max(size(X_ref))-1)*Ref.h;

% Compute velocities from coordinates
Xdot_ref = [Ref.y0(6),diff(X_ref)/h_ref];
Ydot_ref = [Ref.y0(7),diff(Y_ref)/h_ref];
thetadot_ref = [Ref.y0(8),diff(theta_ref)/h_ref];
phidot_ref = [Ref.y0(9),diff(phi_ref)/h_ref];
psidot_ref = [Ref.y0(10),diff(psi_ref)/h_ref];

E_ref = 1/2*m*(Xdot_ref.^2 + Ydot_ref.^2 + R^2*sin(theta_ref).*thetadot_ref.^2) + ...
    1/2*(I_A*(psidot_ref - phidot_ref.*sin(theta_ref)).^2 + ...
        I_T*(thetadot_ref.^2 + phidot_ref.^2.*(cos(theta_ref).^2))) + ...
    m*g*R*cos(theta_ref);

%% Plot graphics

marker_indices = T*h*10;

f = figure;
f.Position(3:4) = [680 420];
plot(t,X,'-o','LineWidth',.5,'MarkerIndices',1:marker_indices:length(X))
hold on
plot(t,Y,'-+','LineWidth',.5,'MarkerIndices',1:marker_indices:length(Y))
hold on
plot(t_ref,X_ref,':','LineWidth',1.5)
hold on
plot(t_ref,Y_ref,':','LineWidth',1.5)
xlim([t0 tf])
ylim('padded')
lgd = legend('$X$','$Y$','$X_{ref}$','$Y_{ref}$','Interpreter','latex','FontSize',16);
lgd.Location = 'northeastoutside';
xlabel('Tiempo (s)')
ylabel('$X,Y$ (m)', 'Interpreter', 'latex','FontSize',16)
ax = gca;
ax.FontSize = 12;
savefig(f,strcat(num2str(exp_n),'XY.fig'))
saveas(f,strcat(num2str(exp_n),'XY'),'eps')
close(f)

f = figure;
f.Position(3:4) = [680 420];
plot(t,theta,'-o','LineWidth',.5,'MarkerIndices',1:marker_indices:length(theta))
hold on
plot(t,phi,'-+','LineWidth',.5,'MarkerIndices',1:marker_indices:length(phi))
hold on
plot(t,psi,'-*','LineWidth',.5,'MarkerIndices',1:marker_indices:length(psi))
hold on
plot(t_ref,theta_ref,':','LineWidth',1.5)
hold on
plot(t_ref,phi_ref,':','LineWidth',1.5)
hold on
plot(t_ref,psi_ref,':','LineWidth',1.5)
ax = gca;
ax.FontSize = 12;
xlim([t0 tf])
ylim('padded')
lgd = legend('$\theta$','$\phi$','$\psi$','$\theta_{ref}$','$\phi_{ref}$','$\psi_{ref}$','Interpreter','latex','FontSize',16);
lgd.Location = 'northeastoutside';
xlabel('Tiempo (s)')
ylabel('$\theta,\phi,\psi$ (rad)', 'Interpreter', 'latex','FontSize',16)
savefig(f,strcat(num2str(exp_n),'thetaphipsi.fig'))
saveas(f,strcat(num2str(exp_n),'thetaphipsi'),'eps')
close(f)

f = figure;
f.Position(3:4) = [680 420];
plot(t,Xdot,'-o','LineWidth',.5,'MarkerIndices',1:marker_indices:length(X))
hold on
plot(t,Ydot,'-+','LineWidth',.5,'MarkerIndices',1:marker_indices:length(Y))
hold on
plot(t_ref,Xdot_ref,':','LineWidth',1.5)
hold on
plot(t_ref,Ydot_ref,':','LineWidth',1.5)
xlim([t0 tf])
ylim('padded')
lgd = legend('$\dot{X}$','$\dot{Y}$','$\dot{X}_{ref}$','$\dot{Y}_{ref}$','Interpreter','latex','FontSize',16);
lgd.Location = 'northeastoutside';
xlabel('Tiempo (s)')
ylabel('$\dot{X},\dot{Y}$ (m/s)', 'Interpreter', 'latex','FontSize',16)
ax = gca;
ax.FontSize = 12;
savefig(f,strcat(num2str(exp_n),'XdotYdot.fig'))
saveas(f,strcat(num2str(exp_n),'XdotYdot'),'eps')
close(f)

f = figure;
f.Position(3:4) = [680 420];
plot(t,thetadot,'-o','LineWidth',.5,'MarkerIndices',1:marker_indices:length(theta))
hold on
plot(t,phidot,'-+','LineWidth',.5,'MarkerIndices',1:marker_indices:length(phi))
hold on
plot(t,psidot,'-*','LineWidth',.5,'MarkerIndices',1:marker_indices:length(psi))
hold on
plot(t_ref,thetadot_ref,':','LineWidth',1.5)
hold on
plot(t_ref,phidot_ref,':','LineWidth',1.5)
hold on
plot(t_ref,psidot_ref,':','LineWidth',1.5)
ax = gca;
ax.FontSize = 12;
xlim([t0 tf])
ylim('padded')
lgd = legend('$\dot{\theta}$','$\dot{\phi}$','$\dot{\psi}$','$\dot{\theta}_{ref}$','$\dot{\phi}_{ref}$','$\dot{\psi}_{ref}$','Interpreter','latex','FontSize',16);
lgd.Location = 'northeastoutside';
xlabel('Tiempo (s)')
ylabel('$\dot{\theta},\dot{\phi},\dot{\psi}$ (rad/s)', 'Interpreter', 'latex','FontSize',16)
savefig(f,strcat(num2str(exp_n),'thetadotphidotpsidot.fig'))
saveas(f,strcat(num2str(exp_n),'thetadotphidotpsidot'),'eps')
close(f)

f = figure;
f.Position(3:4) = [680 420];
plot(t,E,'-o','LineWidth',.5,'MarkerIndices',1:marker_indices:length(E))
hold on
plot(t_ref,E_ref,':','LineWidth',1.5)
ax = gca;
ax.FontSize = 12;
xlim([t0 tf])
ylim('padded')
lgd = legend('Contacto (orden 2)','Referencia','FontSize',12);
lgd.Location = 'northeastoutside';
xlabel('Tiempo (s)')
ylabel('Energía mecánica total (J)')
savefig(f,strcat(num2str(exp_n),'energia.fig'))
saveas(f,strcat(num2str(exp_n),'energia'),'eps')
close(f)
