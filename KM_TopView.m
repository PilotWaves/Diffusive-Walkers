tic
warning('off','MATLAB:hg:willberemoved')
set(0,'DefaultAxesFontSize',20)

T_max = 1300;

%%% Setting up geometry for illustration purposes
theta_geo = [0:2*pi/100:2*pi];
x1 = 1.25*cos(theta_geo);
y1 = 1.25*sin(theta_geo);
x2 = 0.75*cos(theta_geo);
y2 = 0.75*sin(theta_geo);

%%% Dimensional parameters in mm, kg, s

g = 9810;
gammaF = 4.29*g; % dim bath forcing
gamma = 1.047*gammaF;
eps = gamma/gammaF - 1; % reduced acceleration

a = 0.5/sqrt(1.047 - 1); % dimensional amplitude at gamma = 1.047gammaF
A = a*sqrt(eps); % dimensional amplitude

f = 80; % dim frequency
Tf = 2/f; % faraday period

rho = 9.5e-7; % dim oil density
R = 0.28; % dim drop radius
sigma = 0.0206; % dim surface tension

lam_f = 4.75; % dim wavelength
%mu = 1.84e-8; % dimensional viscosity
varphi = pi/2;

%nu = 6*pi*mu/(rho*R^2);

%%% Dimensionless parameters

varphi = f*varphi;
%nu = 2*nu/f;
G = 16*pi^2;
A = A*(2*pi*f)^2/g;
gamma = gamma/g;
Romega = 15; %R/lam_f

%%%  Initializations
y = (A+1)*1.1;
x(1) = 0.02;
vx = 0;
vy = 0;

dt = 0.001;

%%% Initial drop

%%% Illustration purposes
writerObj = VideoWriter('TopView.avi');
writerObj.FrameRate = 75;
open(writerObj);

for t = dt:dt:4 % Discretization of intermediate projectile motion
   
    %%% Path of flight
    eta = y - G*t^2/2 + vy*t;
    Psi = A*sin(pi/2 + 2*pi*t)*sin(pi/2 + 2*pi*x(1)) + gamma*sin(pi/2 + 2*pi*(2*t-varphi));
    
    if eta < Psi + R/(g/(2*pi*f)^2) % Just past impact
       
        %%% Bisection method to find time of intersection
        
        t_low = t;
        t_high = t-dt;
        
        for i = 1:10000
            
            t_impact = (t_low+t_high)/2;
            eta = y - G*t_impact^2/2 + vy*t_impact;
            Psi = A*sin(pi/2 + 2*pi*t_impact)*sin(pi/2 + 2*pi*x(1)) + gamma*sin(pi/2 + 2*pi*(2*t_impact-varphi));
        
            if abs(eta - (Psi + R/(g/(2*pi*f)^2))) < 0.0001
                break
            else if eta < Psi + R/(g/(2*pi*f)^2)
                    t_low = t_impact;
                else
                    t_high = t_impact;
                end
            end
        end
        
        Psi_plot = A*sin(pi/2 + 2*pi*t)*sin(pi/2 + 2*pi*x_plot) + gamma*sin(pi/2 + 2*pi*(2*t-varphi));
       
        %%% Derivatives
    dPsidx = -A*2*pi*sin(pi/2 + 2*pi*t_impact)*sin(2*pi*x(1));
    dPsidt = -A*2*pi*sin(2*pi*t_impact)*sin(pi/2 + 2*pi*x(1)) - 4*pi*gamma*sin(2*pi*(2*t_impact-varphi));
    detadt = vy - G*t_impact;
    
    %%% Angles and Velocities
    theta_in = -pi/2;
    phi = atan(dPsidx*(g/(2*pi*f)^2)/lam_f);
    alpha_in = (theta_in + pi)-phi;
    u_in = abs((detadt - dPsidt)*(g/(2*pi*f)^2)/Tf);
    ut = -u_in*sin(pi/2 + alpha_in);
    un = u_in*sin(alpha_in);
    
        % Coefficient of restitution
            We = rho*R*un^2/sigma;
            We = log10(We);
            
            CRn = -0.077877342024652*We + 0.219769978437467;
            CRt = -0.004433030858978*We^4 - 0.045354932949122*We^3 + ...
                    -0.190613483794666*We^2 - 0.394841202070383*We + ...
                    0.620909597750850;
                
                if CRn < 0.2
                   CRn = 0.2;
                else if CRn > 0.35
                    CRn = 0.35;
                    end
                end
                
                if CRt < 0.2
                   CRt = 0.2;
                else if CRt > 0.78
                        CRt = 0.78;
                    end
                end
           %%%%%%%%%%%%%%%%%%%%%%%
           
      if ut ~= 0
           alpha_out = atan2(CRn*un,CRt*ut);
      else
           alpha_out = pi/2;
      end
                
      theta_out = alpha_out + phi;
      u_out = sqrt((CRn*un)^2 + (CRt*ut)^2);
      
      %%% Iterating velocity and coordinates
      vx = u_out*sin(pi/2 + theta_out)*Tf/lam_f;
      vy = u_out*sin(theta_out)*Tf/(g/(2*pi*f)^2);
      y = Psi+R/(g/(2*pi*f)^2);
      
      break
      
    else
      
    end
        
end

count = 1;
            h = plot(x1,y1,'k',x2,y2,'k',cos(x(count)/Romega),sin(x(count)/Romega),'r.');
            set(h(1),'LineWidth',2);
            set(h(2),'LineWidth',2);
            set(h(3),'MarkerSize',20);
            axis([-1.5 1.5 -1.5 1.5],'square');
            axis off
            frame = getframe(gcf);
            writeVideo(writerObj,frame);

%%% Subsequent drops

Trapped = 0;
T=t_impact;

while T <= T_max
for t = dt:dt:4  % Discretization of intermediate projectile motion
   
    T = T + dt;
    
    %%% Path of flight
    xi = x(count) + vx*t;
    eta = y - G*t^2/2 + vy*t;
    Psi = A*sin(pi/2 + 2*pi*T)*sin(pi/2 + 2*pi*xi) + gamma*sin(pi/2 + 2*pi*(2*T-varphi));

    
    if eta > Psi + R/(g/(2*pi*f)^2)
        Trapped = 0;
        
                
      %%%  Makes video
%            if rem(t,0.1) == 0
%            h = plot(cos(x(1:count)/Romega),sin(x(1:count)/Romega),'o',x1,y1,'k',x2,y2,'k',cos(xi/Romega),sin(xi/Romega),'r.');
%            set(h(1),'MarkerSize',10);
%            set(h(2),'LineWidth',2);
%            set(h(3),'LineWidth',2);
%            set(h(4),'MarkerSize',20);
%            axis([-1.5 1.5 -1.5 1.5],'square');
%            axis off
%            frame = getframe(gcf);
%            writeVideo(writerObj,frame);
%            end
    
    else %%% Just past impact
        
        if Trapped == 1
            
            %%% Chatering fix
            dxidt = vx;
            dPsidx = -2*pi*A*sin(pi/2 + 2*pi*T)*sin(2*pi*xi);
            dPsidt = -2*pi*A*sin(2*pi*T)*sin(pi/2 + 2*pi*xi) - 4*pi*gamma*sin(2*pi*(2*T-varphi)) + dPsidx*dxidt;

            y = Psi + R/(g/(2*pi*f)^2);
            vy = dPsidt;
            
            break
            
        else
       
            if eta == Psi + R/(g/(2*pi*f)^2)
                t_impact = t;
            else
        t_low = t;
        t_high = t-dt;
        T_prev = T-t;
        
        for i = 1:10000
        
            t_impact = (t_low+t_high)/2;
            T = T_prev + t_impact;
            xi = x(count) + vx*t_impact;
            eta = y - G*t_impact^2/2 + vy*t_impact;
            Psi = A*sin(pi/2 + 2*pi*T)*sin(pi/2 + 2*pi*xi) + gamma*sin(pi/2 + 2*pi*(2*T-varphi));
            
            if abs(eta - (Psi + R/(g/(2*pi*f)^2))) < 0.0001
                break
            else if eta < Psi + R/(g/(2*pi*f)^2)
                    t_low = t_impact;
                else
                    t_high = t_impact;
                end
            end
        end
            end
            
            %%% Derivatives
            dxidt = vx;
            detadt = vy - G*t_impact;
            dPsidx = -2*pi*A*sin(pi/2 + 2*pi*T)*sin(2*pi*xi);
            dPsidt = -2*pi*A*sin(2*pi*T)*sin(pi/2 + 2*pi*xi) - 4*pi*gamma*sin(2*pi*(2*T-varphi)) + dPsidx*dxidt;
            
            %%% Rectifying angles at impact; i.e., finding theta

            theta_in = atan2(detadt*(g/(2*pi*f)^2)/Tf,dxidt*lam_f/Tf);
            
            phi = atan(dPsidx*(g/(2*pi*f)^2)/lam_f);
            alpha_in = (theta_in + pi) - phi;
            u_in = sqrt((dxidt*lam_f/Tf)^2 + ((detadt - dPsidt)*(g/(2*pi*f)^2)/Tf)^2);
            ut = -u_in*sin(pi/2 + alpha_in);
            un = u_in*sin(alpha_in);
            
            % Coefficient of restitution
            We = rho*R*un^2/sigma;
            We = log10(We);
            
            CRn = -0.077877342024652*We + 0.219769978437467;
            CRt = -0.004433030858978*We^4 - 0.045354932949122*We^3 + ...
                    -0.190613483794666*We^2 - 0.394841202070383*We + ...
                    0.620909597750850;
                
                if CRn < 0.2
                   CRn = 0.2;
                else if CRn > 0.35
                    CRn = 0.35;
                    end
                end
                
                if CRt < 0.2
                   CRt = 0.2;
                else if CRt > 0.78
                        CRt = 0.78;
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%
                
            alpha_out = atan2(CRn*un,CRt*ut);
                 
            theta_out = alpha_out + phi;
            u_out = sqrt((CRn*un)^2 + (CRt*ut)^2);
            
            %%% Iterating velocity and coordinates
            count = count + 1;
            vx = u_out*sin(pi/2 + theta_out)*Tf/lam_f;
            if abs(vx) < 1e-15
               vx = 0; 
            end
            vy = u_out*sin(theta_out)*Tf/(g/(2*pi*f)^2);
            x(count) = xi;
            y = Psi + R/(g/(2*pi*f)^2);
            
            FA = 0.45;
            if vy < dPsidt
                %vy = dPsidt;
                Trapped = 1;
                FA = 0.4;
            end
            
      %%%  Makes video
           
                    h = plot(cos(x(1:count-1)/Romega),sin(x(1:count-1)/Romega),'o',x1,y1,'k',x2,y2,'k',cos(x(count)/Romega),sin(x(count)/Romega),'r.');
                    set(h(1),'MarkerSize',10);
                    set(h(2),'LineWidth',2);
                    set(h(3),'LineWidth',2);
                    set(h(4),'MarkerSize',20);
                    axis([-1.5 1.5 -1.5 1.5],'square');
                    axis off
                    frame = getframe(gcf);
                    writeVideo(writerObj,frame);
           
                        break
            
        end
            
    end
    
end
end

close(writerObj);
toc