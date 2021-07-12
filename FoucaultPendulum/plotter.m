m = 28;
g = 9.83;
l = 67;

t0 = 0;
tf = 3600;
h = 0.1;
t = t0:h:tf;

alpha = 0.001;

formatSpec = '%f %f';
figHeight = 260;

sizeA = [2 Inf];

filename = strcat('foucault_pendulum-herglotz-alpha=',num2str(alpha));
fileID = fopen(strcat(filename,'.txt'),'r');
Cherglotz = fscanf(fileID,formatSpec,sizeA);

filename = strcat('foucault_pendulum-ld-alpha=',num2str(alpha));
fileID = fopen(strcat(filename,'.txt'),'r');
Cld = fscanf(fileID,formatSpec,sizeA);

filename = strcat('foucault_pendulum-newton-alpha=',num2str(alpha));
fileID = fopen(strcat(filename,'.txt'),'r');
Cnewton = fscanf(fileID,formatSpec,sizeA);

marker_indices = 20;

figure
plot(t,Cherglotz(1,:),'-o','LineWidth',.5,'MarkerIndices',1:marker_indices:length(Cld(1,:)))
hold on
plot(t,Cld(1,:),'-+','LineWidth',.5,'MarkerIndices',1:marker_indices:length(Cld(1,:)))
hold on
plot(t,Cnewton(1,:),':','LineWidth',1.5)
xlim([3500 3600])
ylim([-0.6 0.6])
legend('Contacto','LA','Referencia')
ax = gca;
ax.FontSize = 12;
xlabel('Tiempo (s)')
ylabel('$x$ (m)', 'Interpreter', 'latex', 'FontSize', 16)
% close

figure
plot(t,Cherglotz(2,:),'-o','LineWidth',.5,'MarkerIndices',1:marker_indices:length(Cld(1,:)))
hold on
plot(t,Cld(2,:),'-+','LineWidth',.5,'MarkerIndices',1:marker_indices:length(Cld(1,:)))
hold on
plot(t,Cnewton(2,:),':','LineWidth',1.5)
xlim([3500 3600])
ylim([-0.25 0.25])
legend('Contacto','LA','Referencia')
ax = gca;
ax.FontSize = 12;
xlabel('Tiempo (s)')
ylabel('$y$ (m)', 'Interpreter', 'latex', 'FontSize', 16)
% close

N = 100/h;
from = 35000;
d1 = zeros(1, N);
d2 = zeros(1, N);
for i = 1:N
    j = i + from;
    d1(i) = sqrt((Cherglotz(1,j) - Cnewton(1,j))^2 + ...
        (Cherglotz(2,j) - Cnewton(2,j))^2);
    d2(i) = sqrt((Cld(1,j) - Cnewton(1,j))^2 + ...
        (Cld(2,j) - Cnewton(2,j))^2);
end

figure
plot(((1:N)+from)/10,d1,'-o','MarkerIndices',1:marker_indices:N)
hold on
plot(((1:N)+from)/10,d2,'-+','MarkerIndices',1:marker_indices:N)
legend('Contacto','LDA')
ax = gca;
ax.FontSize = 12;
xlabel('Tiempo (s)')
ylabel('$\|q_{t_j} - q_{ref}(t_j)\|$','Interpreter','latex','FontSize', 16)

X = Cherglotz(1,from+1:end);
Y = Cherglotz(2,from+1:end);
Xdot = diff(Cherglotz(1,from:end))/h;
Ydot = diff(Cherglotz(2,from:end))/h;
Eherglotz = (1/2)*m*(Xdot.^2 + Ydot.^2) + ...
    (1/2)*m*(g/l)*(X.^2 + Y.^2);

X = Cld(1,from+1:end);
Y = Cld(2,from+1:end);
Xdot = diff(Cld(1,from:end))/h;
Ydot = diff(Cld(2,from:end))/h;
Eld = (1/2)*m*(Xdot.^2 + Ydot.^2) + (1/2)*m*(g/l)*(X.^2 + Y.^2);

X = Cnewton(1,from+1:end);
Y = Cnewton(2,from+1:end);
Xdot = diff(Cnewton(1,from:end))/h;
Ydot = diff(Cnewton(2,from:end))/h;
Enewton = (1/2)*m*(Xdot.^2 + Ydot.^2) + (1/2)*m*(g/l)*(X.^2 + Y.^2);

figure
plot(3500:h:3600,Eherglotz,'-o','LineWidth',.5,'MarkerIndices',1:marker_indices:N)
hold on
plot(3500:h:3600,Eld,'-+','LineWidth',.5,'MarkerIndices',1:marker_indices:N)
hold on
plot(3500:h:3600,Enewton,':','LineWidth',1.5)
legend('Contacto','LDA','Referencia')
ax = gca;
ax.FontSize = 12;
xlabel('Tiempo (s)')
ylabel('Energía mecánica total (J)')

% Trayectorias

figure
plot(Cherglotz(1,:),Cherglotz(2,:),'Color',[0 0.4470 0.7410])
axis equal
xlim([-.8 .8])
ylim([-.25 .25])
ax = gca;
ax.FontSize = 12;
xlabel('$x$ (m)', 'fontsize', 18, 'Interpreter', 'latex')
ylabel('$y$ (m)', 'fontsize', 18, 'Interpreter', 'latex')
fig = gcf;
fig.Position = [fig.Position(1:3),figHeight];

figure
plot(Cld(1,:),Cld(2,:),'Color',[0.8500 0.3250 0.0980])
axis equal
xlim([-.8 .8])
ylim([-.25 .25])
ax = gca;
ax.FontSize = 12;
xlabel('$x$ (m)', 'fontsize', 18, 'Interpreter', 'latex')
ylabel('$y$ (m)', 'fontsize', 18, 'Interpreter', 'latex')
fig = gcf;
fig.Position = [fig.Position(1:3),figHeight];

figure
plot(Cnewton(1,:),Cnewton(2,:),'Color',[0.9290 0.6940 0.1250])
axis equal
xlim([-.8 .8])
ylim([-.25 .25])
ax = gca;
ax.FontSize = 12;
xlabel('$x$ (m)', 'fontsize', 18, 'Interpreter', 'latex')
ylabel('$y$ (m)', 'fontsize', 18, 'Interpreter', 'latex')
fig = gcf;
fig.Position = [fig.Position(1:3),figHeight];
