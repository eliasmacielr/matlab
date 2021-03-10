m = 28;
g = 9.83;
l = 67;
h = 0.1;
alpha = 0.0001;
formatSpec = ['(' '%f' ', ' '%f' ')'];
figHeight = 260;

filename = strcat('foucault_pendulum-newton-alpha=',num2str(alpha));
fileID = fopen(strcat(filename,'.txt'),'r');
Cnewton = textscan(fileID,formatSpec);

figure
plot(Cnewton{1},Cnewton{2},'Color',[0.4660 0.6740 0.1880])
axis equal
xlim([-.8 .8])
ylim([-.25 .25])
xlabel('$x$ (m)', 'fontsize', 18, 'Interpreter', 'latex')
ylabel('$y$ (m)', 'fontsize', 18, 'Interpreter', 'latex')
fig = gcf;
fig.Position = [fig.Position(1:3),figHeight];
savefig(strcat(filename,'.fig'))
%close


filename = strcat('foucault_pendulum-herglotz-alpha=',num2str(alpha),...
    '-60min-h0.1');
fileID = fopen(strcat(filename,'.txt'),'r');
Cherglotz = textscan(fileID,formatSpec);

figure
plot(Cherglotz{1},Cherglotz{2},'Color',[0 0.4470 0.7410])
axis equal
xlim([-.8 .8])
ylim([-.25 .25])
xlabel('$x$ (m)', 'fontsize', 18, 'Interpreter', 'latex')
ylabel('$y$ (m)', 'fontsize', 18, 'Interpreter', 'latex')
fig = gcf;
fig.Position = [fig.Position(1:3),figHeight];
savefig(strcat(filename,'.fig'))
%close


filename = strcat('foucault_pendulum-ld-alpha=',num2str(alpha),...
    '-60min-h0.1');
fileID = fopen(strcat(filename,'.txt'),'r');
Cld = textscan(fileID,formatSpec);

figure
plot(Cld{1},Cld{2},'Color',[0.8500 0.3250 0.0980])
axis equal
xlim([-.8 .8])
ylim([-.25 .25])
xlabel('$x$ (m)', 'fontsize', 18, 'Interpreter', 'latex')
ylabel('$y$ (m)', 'fontsize', 18, 'Interpreter', 'latex')
fig = gcf;
fig.Position = [fig.Position(1:3),figHeight];
savefig(strcat(filename,'.fig'))
%close

N = 100/h;
from = 35000;
d1 = zeros(1, N);
d2 = zeros(1, N);
for i = 1:N
    j = i + from;
    d1(i) = sqrt((Cherglotz{1}(j) - Cnewton{1}(j))^2 + ...
        (Cherglotz{2}(j) - Cnewton{2}(j))^2);
    d2(i) = sqrt((Cld{1}(j) - Cnewton{1}(j))^2 + ...
        (Cld{2}(j) - Cnewton{2}(j))^2);
end

% figure
% plot((1:N)+from,d1,(1:N)+from,d2)
% legend('Contacto','LDA')
% xlabel('$t$ (s)','fontsize',18,'Interpreter','latex')
% ylabel('$\|q_{t_j} - q_{ref}(t_j)\|$','fontsize',18,'Interpreter','latex')
% savefig(strcat('Error relativo alpha=',num2str(alpha),'.fig'))
% close
% 
% X = Cnewton{1}(from+1:end);
% Y = Cnewton{2}(from+1:end);
% Xdot = diff(Cnewton{1}(from:end))/h;
% Ydot = diff(Cnewton{2}(from:end))/h;
% Enewton = (1/2)*m*(Xdot.^2 + Ydot.^2) + (1/2)*m*(g/l)*(X.^2 + Y.^2);
% 
% X = Cherglotz{1}(from+1:end);
% Y = Cherglotz{2}(from+1:end);
% Xdot = diff(Cherglotz{1}(from:end))/h;
% Ydot = diff(Cherglotz{2}(from:end))/h;
% Eherglotz = (1/2)*m*(Xdot.^2 + Ydot.^2) + (1/2)*m*(g/l)*(X.^2 + Y.^2);
% 
% X = Cld{1}(from+1:end);
% Y = Cld{2}(from+1:end);
% Xdot = diff(Cld{1}(from:end))/h;
% Ydot = diff(Cld{2}(from:end))/h;
% Eld = (1/2)*m*(Xdot.^2 + Ydot.^2) + (1/2)*m*(g/l)*(X.^2 + Y.^2);
% 
% figure
% plot((1:N)+from,Eherglotz,'Color',[0 0.4470 0.7410],'Marker','s')
% hold on
% plot((1:N)+from,Eld,'Color',[0.8500 0.3250 0.0980],'Marker','*')
% hold on
% plot((1:N)+from,Enewton,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
% legend('Contacto','LDA','Referencia')
% xlabel('Tiempo (s)','fontsize',12)
% ylabel('Energía mecánica total (J)','fontsize',12)
% savefig(strcat('Energia alpha=',num2str(alpha),'.fig'))
% close
