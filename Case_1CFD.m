%Case 1 CFD code
%Givens
clc;
clear;
L=0.02;
Ta=100;
Tb=200;
k=0.5;
S=1000000;
n=5;
dx=L/n;
a = zeros(n);
Error=0;
AbsError=0;
Sa = S*dx +(2*k*Ta/dx);



%a(i)x(i) = b(i)*x(i+1) + c(i)*x(i-1) + d(i);

% generating tri diagonal matrix
for i=1:n
    
        if(i==1)
            a(i,1)= (3*k/dx);
            a(i,2)= (-1*k/dx);
            for j=3:n
                a(i,j)=0;
            end
            d(i,1)= Sa;
        end
        
        if(i==n)
            a(i,n)= (3*k/dx);
            a(i,n-1)= (-1*k/dx);
            for j=1:n-2
                a(i,j)=0;
            end
            d(i,1)= S*dx + (2*k*Tb/dx);
        end
        
        %Everything besides the boundary values
        if(i>1 && i<n)
            for j=1:n
                a(i,j)=0;
            end
            a(i,i-1)= (-1*k/dx);
            a(i,i)= (2*k/dx);
            a(i,i+1)= (-1*k/dx);
            d(i,1)= S*dx;
        end
            
end


%Thomas algorithm - Gauss Elimination
for i=1:n
    if(i==1)
        a(i,2)=a(i,2)/a(i,1);
        d(1,1)=d(1,1)/a(1,1);
        a(1,1)=1;
    end
     if(i==n)
        a(i,i)=a(i,i)/a(i,i-1);
        d(i,1)=d(i,1)/a(i,i-1);
        a(i,i-1)=1;       
        d(i,i-1)=0;
        a(i,i)=a(i,i)-a(i-1,i);
        d(i,1)=d(i,1)-d(i-1,1);
        
        d(i,1)=d(i,1)/a(i,i);
        a(i,i)=1;
        
    end
    
    if(i>1 && i<n)
        a(i,i)=a(i,i)/a(i,i-1);
        a(i,i+1)=a(i,i+1)/a(i,i-1);
        d(i,1)=d(i,1)/a(i,i-1);
        a(i,i-1)=1;
        
        a(i,i-1)=0;
        a(i,i)=a(i,i)-a(i-1,i);
        a(i,i+1)=a(i,i+1)-a(i-1,i+1);
        d(i,1)=d(i,1)-d(i-1,1);
        
        
        a(i,i+1)=a(i,i+1)/a(i,i);
        d(i,1)=d(i,1)/a(i,i);
        a(i,i)= 1;
    end
    
         
end


% Calculation of temperature distribution and error
for i=n:-1:1
    
    if(i==n)
        T(n+1)=d(i);      
        
        T(1)=Ta; T(n+2)=Tb;        
        Tempactual(1)=Ta; Tempactual(n+2)=Tb;        
        u(1)=0; u(n+2)=L;        
        Errper(1)=0; Errper(n+2)=0;
    end
    
    if(i<n)
        T(i+1)=d(i)-(T(i+2)*a(i,i+1));      
    end    
    u(i+1)=dx/2 + dx*(i-1);
    
    Tempactual(i+1) = Ta + ((Tb-Ta)/L + S*(L-u(i+1))/(2*k)) * u(i+1);
    Error=T(i+1) - Tempactual(i+1);
    AbsError=AbsError+Error^2
    Errper(i+1)=100 * Error/Tempactual(i+1);
    
end

%Output
u'
Tempactual'
T'
Errper'
k=sqrt(AbsError/n)


% Plot
plot(u,T,'ro');
hold on;
plot(u,Tempactual, 'g');
hold on;
plot(u,Errper, 'b');

