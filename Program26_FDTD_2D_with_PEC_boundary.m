%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written for Course :- Computational Electromagnetics, Fall 2011
%                       Department of Electrical Engineering
%                       Indian Institute of Technology Madras
%                       Chennai - 600036, India
%
% Authors            :- Sathya Swaroop Ganta, B.Tech., M.Tech. Electrical Engg.
%                       Kayatri, M.S. Engineering Design
%                       Pankaj, M.S. Electrical Engg.
%                       Sumanthra Chaudhuri, M.S. Electrical Engg.
%                       Projesh Basu, M.S. Electrical Engg.
%                       Nikhil Kumar CS, M.S. Electrical Engg.
%
% Instructor :- Ananth Krishnan
%               Assistant Professor
%               Department of Electrical Engineering
%               Indian Institute of Technology Madras
%
% Any correspondance regarding this program may be addressed to
% Prof. Ananth Krishnan at 'computational.em.at.iit.madras@gmail.com'
%
% Copyright/Licensing :- For educational and research purposes only. No
% part of this program may be used for any financial benefit of any kind
% without the consent of the instructor of the course.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "2D FDTD solution for a boundary condition of Perfect Electric Conductor"
% 
% Objective of the program is to solve for the Maxwell's equation for a TM 
% wave containing the xy-plane polarized magnetic field having components Hy
% and Hx and z-polarized electric field Ez. The fields are updated at every 
% timestep, in a space, where all physical parameters of free space are not
% normalized to 1 but are given real and known values. The update is done 
% using standard update equations obtained from the difference form of Maxwell's 
% curl equations with very low electric and magnetic conductivities of 4x10^(-4) 
% units incorporated in them. The field points are defined in a grid described 
% by Yee's algorithm. The H fields are defined at every half coordinate of 
% spacesteps. More precisely, the Hx part is defined at every half y-coordinate 
% and full x-coordinate and the Hy part is defined at every half x-coordinate 
% and full y-coordinate and E fields i.e the Ez part is defined at every full 
% x and full y-coordinate points.Also here, the space-step length is taken 
% as 1 micron instead of 1 unit in unitless domain assumed in previous programs. 
% Also, the time update is done using Leapfrog time-stepping. Here, H-fields
% i.e. Hx and Hy are updated every half time-step and E fields i.e Ez are 
% updated every full time-step. This is shown by two alternating vector updates 
% spanning only a part of spatial grid where the wave, starting from source, 
% has reached at that particular time instant avoiding field updates at all 
% points in the grid which is unnecessary at that time instant. These spatial 
% updates are inside the main for-loop for time update, spanning the entire 
% time grid. Also, here, the matrices used as multiplication factors for update 
% equations are initialized before the loop starts to avoid repeated calculation 
% of the same in every loop iteration, a minor attempt at optimization. The 
% boundary condition here is perfect electric boundary in the sense that the 
% boundary grid points have zero electric field values irrespective of the 
% influence of external fields.
%
% A source of electric field is defined at the center of the spatial domain, 
% which is a hard source, in that it does not change its value due to 
% interference from external fields i.e in other words, the source is a 
% perfect electric conductor. The form of the source can be varied using 
% the variables sine, gaussian and impulse. The source is available in four 
% standard forms- Unit-time step, Impulse, Gausian and Sinusoidal forms. 
% The color scaled plot of Ez field over the entire spatial domain is shown 
% at every time step. The simulation can be ended by closing this plot window 
% or by waiting till all the time step updates are completed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Clearing variables in memory and Matlab command screen
clear all;
clc;
% Grid Dimension in x (xdim) and y (ydim) directions
xdim=200;
ydim=200;
%Total no of time steps
time_tot=350;
%Position of the source (center of the domain)
xsource=100;
ysource=100;
%Courant stability factor
S=1/(2^0.5);
% Parameters of free space (permittivity and permeability and speed of
% light) are all not 1 and are given real values
epsilon0=(1/(36*pi))*1e-9;
mu0=4*pi*1e-7;
c=3e+8;
% Spatial grid step length (spatial grid step= 1 micron and can be changed)
delta=1e-6;
% Temporal grid step obtained using Courant condition
deltat=S*delta/c;
% Initialization of field matrices
Ez=zeros(xdim,ydim);
Hy=zeros(xdim,ydim);
Hx=zeros(xdim,ydim);
% Initialization of permittivity and permeability matrices
epsilon=epsilon0*ones(xdim,ydim);
mu=mu0*ones(xdim,ydim);
% Initializing electric and magnetic conductivity matrices
sigma=4e-4*ones(xdim,ydim);
sigma_star=4e-4*ones(xdim,ydim);
%Choice of nature of source
gaussian=0;
sine=0;
% The user can give a frequency of his choice for sinusoidal (if sine=1 above) waves in Hz 
frequency=1.5e+13;
impulse=0;
%Choose any one as 1 and rest as 0. Default (when all are 0): Unit time step
%Multiplication factor matrices for H matrix update to avoid being calculated many times 
%in the time update loop so as to increase computation speed
A=((mu-0.5*deltat*sigma_star)./(mu+0.5*deltat*sigma_star)); 
B=(deltat/delta)./(mu+0.5*deltat*sigma_star);
                          
%Multiplication factor matrices for E matrix update to avoid being calculated many times 
%in the time update loop so as to increase computation speed                          
C=((epsilon-0.5*deltat*sigma)./(epsilon+0.5*deltat*sigma)); 
D=(deltat/delta)./(epsilon+0.5*deltat*sigma);                     
% Update loop begins
for n=1:1:time_tot
    
    %if source is impulse or unit-time step 
    if gaussian==0 && sine==0 && n==1
        Ez(xsource,ysource)=1;
    end
    
    % Setting time dependent boundaries to update only relevant parts of the 
    % vector where the wave has reached to avoid unnecessary updates.
    if n<xsource-2
        n1=xsource-n-1;
    else
        n1=1;
    end
    if n<xdim-1-xsource
        n2=xsource+n;
    else
        n2=xdim-1;
    end
    if n<ysource-2
        n11=ysource-n-1;
    else
        n11=1;
    end
    if n<ydim-1-ysource
        n21=ysource+n;
    else
        n21=ydim-1;
    end
    
    %Vector update instead of for-loop for Hy and Hx fields
    Hy(n1:n2,n11:n21)=A(n1:n2,n11:n21).*Hy(n1:n2,n11:n21)+B(n1:n2,n11:n21).*(Ez(n1+1:n2+1,n11:n21)-Ez(n1:n2,n11:n21));
    Hx(n1:n2,n11:n21)=A(n1:n2,n11:n21).*Hx(n1:n2,n11:n21)-B(n1:n2,n11:n21).*(Ez(n1:n2,n11+1:n21+1)-Ez(n1:n2,n11:n21));
    
    %Vector update instead of for-loop for Ez field
    Ez(n1+1:n2+1,n11+1:n21+1)=C(n1+1:n2+1,n11+1:n21+1).*Ez(n1+1:n2+1,n11+1:n21+1)+D(n1+1:n2+1,n11+1:n21+1).*(Hy(n1+1:n2+1,n11+1:n21+1)-Hy(n1:n2,n11+1:n21+1)-Hx(n1+1:n2+1,n11+1:n21+1)+Hx(n1+1:n2+1,n11:n21));
    
    % Perfect Electric Conductor boundary condition
    Ez(1:xdim,1)=0;
    Ez(1:xdim,ydim)=0;
    Ez(1,1:ydim)=0;
    Ez(xdim,1:ydim)=0;
    
    % Source conditions
    if impulse==0
        % If unit-time step
        if gaussian==0 && sine==0
            Ez(xsource,ysource)=1;
        end
        %if sine
        if sine==1
            tstart=1;
            N_lambda=c/(frequency*delta);
            Ez(xsource,ysource)=sin(((2*pi*(c/(delta*N_lambda))*(n-tstart)*deltat)));
        end
        %if gaussian
        if gaussian==1
            if n<=42
                Ez(xsource,ysource)=(10-15*cos(n*pi/20)+6*cos(2*n*pi/20)-cos(3*n*pi/20))/32;
            else
                Ez(xsource,ysource)=0;
            end
        end
    else
        %if impulse
        Ez(xsource,ysource)=0;
    end
    
    %Movie type colour scaled image plot of Ez
    imagesc(delta*(1:1:xdim)*1e+6,(1e+6*delta*(1:1:ydim))',Ez',[-1,1]);colorbar;
    title(['\fontsize{20}Colour-scaled image plot of Ez in a spatial domain with PEC boundary and at time = ',num2str(round(n*deltat*1e+15)),' fs']); 
    xlabel('x (in um)','FontSize',20);
    ylabel('y (in um)','FontSize',20);
    set(gca,'FontSize',20);
    getframe;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF PROGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
