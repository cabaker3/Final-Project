%Godunov Scheme..............
clear all
clc
 
%Constants:
L = 1;
g = 1.4;
dx = 1/100;
dt = 1/500;
i = 1;
IM = L/dx;
 
%Lax Initial Conditions:
for x=0:dx:1
    if x<=0.5
        Qi(:,i) = [0.445; 0.311; 8.928];
        %Qi(:,i) = [1; 0; 1];
    elseif x>=0.5 
        Qi(:,i) = [0.5; 0; 1.4275];
        %Qi(:,i) = [0.125; 0; 0.1];
    end
    i = i+1;
end
Qold(:,:) = Qi(:,:);
Qnew(:,:) = Qold(:,:);
%Initial Flow Properties:
%Density.
for i = 1:IM+1
    rhoi(i,:) = Qi(1,i); 
end
%Velocity.
for i = 1:IM+1
    ui(i,:) = Qi(2,i) / rhoi(i,:);
end
%Total Energy.
for i = 1:IM+1
    eti(i,:) = Qi(3,i) / rhoi(i,:);
end
%Pressure, from the equation of state.
for i = 1:IM+1
    pi(i,:) = (g-1).*(rhoi(i,:) .* eti(i,:) - 0.5 .* rhoi(i,:).*ui(i,:).^2) ;
end
 
%Speed of Sound
for i = 1:IM+1
    a(i,:) = sqrt(g*pi(i,:)/rhoi(i,:));
end
 
%Intial E Matrix:
for i = 1:IM+1
    E(:,i) = [rhoi(i,:)*ui(i,:); rhoi(i,:)*(ui(i,:)^2) + pi(i,:);...
        eti(i,:)*rhoi(i,:)*ui(i,:) + pi(i,:)*ui(i,:)];
end
 
%Eigenvalues
    for i = 1:IM+1
        eigen(:,i) = [ui(i,:); ui(i,:) + a(i,:); ui(i,:) - a(i,:)]; 
    end
 
%Alpha
alpha = max(abs(eigen))'; 
 
k=1;
 
for t = 0:dt:0.16
%.......................................................................................
%Flux
for j = 1:3
    for i = 1:IM
        F(j,i) = (0.5)*(E(j,i) + E(j,i+1)) - (0.5)*(abs(alpha(i,:)))*(Qold(j,i+1) - Qnew(j,i));
    end
end
 
%Qnew
for j = 1:3
    for i = 2:IM
        Qn1(j,i) = Qold(j,i) - (dt/dx)*(F(j,i) - F(j,i-1));
    end
end
 
Qn1(:,1) = Qi(:,1);
Qn1(:,101) = Qi(:,101);
Qnew = Qn1;
Qold = Qnew;
 
%Density.
for i = 1:IM+1
    rho(i,:) = Qnew(1,i); 
end
%Velocity.
for i = 1:IM+1
    u(i,:) = Qnew(2,i) / rho(i,:);
end
%Total Energy.
for i = 1:IM+1
    et(i,:) = Qnew(3,i) / rho(i,:);
end
%Pressure, from the equation of state.
for i = 1:IM+1
    p(i,:) = (g-1).*(rho(i,:) .* et(i,:) - 0.5 .* rho(i,:).*u(i,:).^2) ;
end
 
%Speed of Sound
for i = 1:IM+1
    a(i,:) = sqrt(g*p(i,:)/rho(i,:));
end
 
%Intial E Matrix:
for i = 1:IM+1
    E(:,i) = [rho(i,:)*u(i,:); rho(i,:)*(u(i,:)^2) + p(i,:);...
        et(i,:)*rho(i,:)*u(i,:) + p(i,:)*u(i,:)];
end
 
%Eigenvalues
    for i = 1:IM+1
        eigen(:,i) = [u(i,:); u(i,:) + a(i,:); u(i,:) - a(i,:)]; 
    end
 
%Alpha
alpha = max(abs(eigen))'; 
%.................................................................................
um(:,k) = u(:,1);
rhom(:,k) = rho(:,1);
pm(:,k) = p(:,1);
etm(:,k) = et(:,1);
k = k+1;
end
 
figure(1)
plot(0:dx:1,ui,'k',0:dx:1,um(:,21),'k--sq',0:dx:1,um(:,41),'k--^',0:dx:1,um(:,81),'k--v')
xlabel('X-Domain')
ylabel('Velcoity')
title('Shock Tube Velocity Profile, CFL = 0.8')
 
figure(2)
plot(0:dx:1,rhoi,'k',0:dx:1,rhom(:,21),'k--sq',0:dx:1,rhom(:,41),'k--^',0:dx:1,rhom(:,81),'k--v')
xlabel('X-Domain')
ylabel('Density')
title('Shock Tube Density Profile, CFL = 0.8')
 
figure(3)
plot(0:dx:1,pi,'k',0:dx:1,pm(:,21),'k--sq',0:dx:1,pm(:,41),'k--^',0:dx:1,pm(:,81),'k--v')
xlabel('X-Domain')
ylabel('Pressure')
title('Shock Tube Pressure Profile, CFL = 0.8')
