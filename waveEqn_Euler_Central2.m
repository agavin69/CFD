clear all;
clc;

format long;

% switch on frame saving
saveframes=0;

% solve the advection equation: du/dt = c*du/dt
% numerical method: RK4 in time, central-2 differencing in space

% number of data points
M = 320;
CFL=0.05;
% periodic interval dimensions (i.e. x \in [a,b))
a = 0; b =  2;

% data point spacing
dx = (b-a)/M;

% coordinates of data points ("nodes")
x = linspace(a, b-dx, M);

% start time
StartTime = 0;
EndTime = 8;

% dt chosen small to make space errors dominant
dt = CFL*dx;

% choose order of RK integrator
s = 4;

% compute number of time steps
Nsteps = ceil((EndTime-StartTime)/dt);

% modify dt to end exactly at EndTime
dt = (EndTime-StartTime)/Nsteps;

% initial condition 
u = exactsolution(x);

t = StartTime; % set time variable
fram = 1;
for n=1:Nsteps % begin time stepping
    
    k1 = u;
    
    up2 = [k1(3:M),k1(1:2)];
    up1 = [k1(2:M),k1(1)];   % k1_{m+1}
    um1 = [k1(M),k1(1:M-1)]; % k1_{m-1}
    um2 = [k1(M-1),k1(M),k1(1:M-2)];
    fa = 1/12*(up2-8*up1+8*um1-um2)/dx;
    k2  = u + 0.5*dt*fa;     % compute k2 using dt/2

    up2 = [k2(3:M),k2(1:2)];
    up1 = [k2(2:M),k2(1)];   % k2_{m+1}
    um1 = [k2(M),k2(1:M-1)]; % k2_{m-1}
    um2 = [k2(M-1),k2(M),k2(1:M-2)];
    fb = 1/12*(up2-8*up1+8*um1-um2)/dx;
    k3  = u + 0.5*dt*fb;     % compute k3 using dt/2
    
    up2 = [k3(3:M),k3(1:2)];
    up1 = [k3(2:M),k3(1)];   % k3_{m+1}
    um1 = [k3(M),k3(1:M-1)]; % k3_{m-1}
    um2 = [k3(M-1),k3(M),k3(1:M-2)];
    fc = 1/12*(up2-8*up1+8*um1-um2)/dx;
    k4  = u + dt*fc;         % compute k4 using dt 
    
    up2 = [k4(3:M),k4(1:2)];
    up1 = [k4(2:M),k4(1)];   % k4_{m+1}
    um1 = [k4(M),k4(1:M-1)]; % k4_{m-1}
    um2 = [k4(M-1),k4(M),k4(1:M-2)];
    fd = 1/12*(up2-8*up1+8*um1-um2)/dx;
    
    %Compute final value
    utilde = u + (dt/6)*(fa + 2*fb + 2*fc + fd);
        
    u = utilde; % finish RK step
    t  = t+dt;   % update time
    
    if(mod(n,40)==0 |n==Nsteps) % selective plotting
    
        plot(x, u, 'r-*', 'LineWidth', 1); % plot numerical solution
        hold on;

        xplot = linspace(a, b, 1000); % plot exact solution at lots of points
        xmod = mod(xplot+t,b); % coordinate for exact pulse
        uexact = exactsolution(xmod);
        plot(xplot, uexact, 'k-', 'LineWidth', 1);hold off;
        axis([0 2 -0.5 1.5]);
        legend(sprintf('Numerical solution (M=%d,t=%g)', M,t), 'Exact solution');
        title('CENTRAL DIFFERENCE - ORDER 4');
        drawnow; pause(0.02);
        if(saveframes == 1)        % save frame to file
            if(fram<10)            fname = sprintf('anim_000%d.ppm', fram);
            elseif(fram<100)   fname = sprintf('anim_00%d.ppm', fram);
            elseif(fram<1000) fname = sprintf('anim_0%d.ppm', fram);
            end
            print('-dppm', fname);
            fram = fram+1;
        end
    end
end

% to make the animated gif: I used ImageMagick's convert
% under cygwin:   convert *ppm deltaplus.gif

% Compute and display the max error
xmod = mod(x+t,b);
uexact    = exactsolution(xmod);
finalerror = abs(uexact-u);
L2=sqrt(sum((uexact(1:M)-u(1:M)).^2));
max_error = max(abs(finalerror));
disp(max_error);