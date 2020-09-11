%  Daniel Kawano, Rose-Hulman Institute of Technology

%  Adapted by Elias Maciel to the case of the disk described in "A Disk
%  Rolling on a Horizontal Surface Without Slip", Paris and Zhang, 2002.
%  Last modified:  Sep 10, 2020

function animate_rolling_disk(X,Y,theta,phi,psi,R,h)

%  Create a circle for the disk:
angle = linspace(0, 2*pi, 36)';
circX = R*cos(angle);
circZ = R*sin(angle);

for k = 1:length(psi)
    R1 = [cos(phi(k)), sin(phi(k)), 0;              %  3-2-1 set of
          -sin(phi(k)), cos(phi(k)), 0;             %  Euler angles
          0, 0, 1];
    R2 = [cos(theta(k)), 0, -sin(theta(k));
          0, 1, 0;
          sin(theta(k)), 0, cos(theta(k))];
    R3 = [1, 0, 0;
          0, cos(psi(k)), sin(psi(k));
          0, -sin(psi(k)), cos(psi(k))];

    Q = R3*R2*R1;
    % Reference configuration
    e1 = Q*[0; -1; 0];                       %  e1
    e2 = Q*[1; 0; 0];                        %  e2
    e3 = Q*[0; 0; 1];                        %  e3

    % Disk
    xcirc(:,k) = X(k)                 + e1(2)*circX + e3(2)*circZ;  %  m
    ycirc(:,k) = Y(k)                 + e1(1)*circX + e3(1)*circZ;  %  m
    zcirc(:,k) = R*abs(cos(theta(k))) + e1(3)*circX + e3(3)*circZ;  %  m

    % Path
    xP(k,1) = X(k);                                 %  m
    yP(k,1) = Y(k);                                 %  m
    zP(k,1) = 0;                                    %  m

    % Point
    xA(k,1) = X(k)                 - R*e3(2);       %  m
    yA(k,1) = Y(k)                 - R*e3(1);       %  m
    zA(k,1) = R*abs(cos(theta(k))) - R*e3(3);       %  m
end

%  Set up the figure window:

figure
set(gcf, 'color', 'w')
plot3(X(1), Y(1), R*abs(cos(theta(1))));
xlabel('\itX (m)')
set(gca, 'xdir', 'reverse')
ylabel('\itY (m)')
set(gca, 'ydir', 'reverse')
zlabel('\itZ (m)            ', 'rotation', 0)
axis equal
xlim([min(X)-R, max(X)+R])
ylim([min(Y)-R, max(Y)+R])
zlim([0, 1.2*(2*R)])
grid on

%  Trace out the path of the instantaneous contact point P, and orient the 
%  disk appropriately.  Also, highlight a material point A on the disk's 
%  periphery:

path = line('xdata', xP(1:1), 'ydata', yP(1:1), 'zdata', zP(1:1), ...
    'color', 'b', 'linewidth', 2);
disk = patch('xdata', xcirc(:,1), 'ydata', ycirc(:,1), ...
    'zdata', zcirc(:,1), 'facecolor', [1, 0.8, 0.6], 'linewidth', 2);
pointA = line('xdata', xA(1), 'ydata', yA(1), 'zdata', zA(1), ...
    'marker', 'o', 'color', 'R', 'markerfacecolor', 'R', 'linewidth', 3);

%  Animate the disk's motion by updating the figure with its current 
%  location and orientation:

% pause

animation = VideoWriter('rolling-disk-contact.avi');
animation.FrameRate = 1/h;
open(animation);

for k = 1:length(xP)
    set(path, 'xdata', xP(1:k), 'ydata', yP(1:k), 'zdata', zP(1:k));
    set(disk, 'xdata', xcirc(:,k), 'ydata', ycirc(:,k), ...
        'zdata', zcirc(:,k));
    set(pointA, 'xdata', xA(k), 'ydata', yA(k), 'zdata', zA(k));
    drawnow   
    writeVideo(animation, getframe(gcf));
end 

close(animation);
