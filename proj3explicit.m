clear all
clc;

dt=0.05;
Time=5;
gamma=0.1;
N=50;
w=Time/dt;
a=1;
dx=1/N;
R=a*dt/dx;
Z=gamma*dt/(dx^2);


u=zeros(w+1,N+2);
u(1,1)=0.1;

u(:,N+2)=1;

for T=2:w+1
   t=dt*(T-1);
    u(T,1)=0.1*(1+sin(6*t));
  for i=3:N 
      
   unew(i)=u(T-1,i)-R*(u(T-1,i+1)-2*u(T-1,i)+u(T-1,i-1))+Z*(u(T-1,i+1)-2*u(T-1,i)+u(T-1,i-1));
   u(T,i)=unew(i);
  end
   u(T,2)=u(T-1,2)-R*(u(T-1,3)-2*u(T-1,2)+u(T-1,1))+Z*(u(T-1,3)-2*u(T-1,2)+u(T-1,1));
   u(T,N+1)=u(T-1,N)-R*(u(T-1,N+1)-2*u(T-1,N)+u(T-1,N-1))+Z*(u(T-1,N+1)-2*u(T-1,N)+u(T-1,N-1));
  
end
 

