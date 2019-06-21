% https://stackoverflow.com/questions/24409045/finite-difference-time-domain-ftdt-method-for-1d-em-wave

%FDTD Yee algorithm to solve coupled EM wave equations
clear
clc


G=50;           %Specify size of the grid
f=10^3;         %choose initial frequency of wave in hertz
e=1;            %specify permitivity and permeability (normalised condition)
u=1;            
Nt=150;         %time steps
E0=1;           %Electric Field initial amplitude

%Specify step sizes using corruant condition
c=3*10^8;
dx=0.01;
dt=dx/2*c;

%make constant terms
c1=-dt/(dx*e);
c2=-dt/(dx*u);

%create vgector place holders
Ex=zeros(1,G);
Hy=zeros(1,G);

%create updating loop
M=moviein(Nt);
for t=1:Nt

    % Spatial Ex

    for k=2:G-1
        Ex(k)=Ex(k)+c1*(Hy(k)-Hy(k-1));
    end
    Ex(G)=0; %PEC boundary condition
    %E Source at LHS boundary 
    Ex(1)=E0*sin(2*pi*f*t);
    %Spatial Hy
    for n=1:G-1
        Hy(n)=Hy(n)+c2*(Ex(n)-Ex(n+1));
    end
    Hy(G)=0; %PMC boundary condition

plot(Ex);
M(:,t) = getframe;
end
movie(M,1);
