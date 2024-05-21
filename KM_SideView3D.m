tic
warning('off','MATLAB:hg:willberemoved')
set(0,'DefaultAxesFontSize',20)

T_max = 11;
%T_max = 22;

%%% Dimensional parameters in mm, kg, s

g = 9810;
gammaF = 4.29*g; % dim bath forcing
gamma = 1.047*gammaF;
%gamma = 1.0055*gammaF;
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

%%%  Initializations
y = (A+1)*1.1;
x0 = -0.15;
x = x0;
vx = 0;
vy = 0;

dt = 0.002;

%%% Initial drop

%%% Illustration purposes
writerObj = VideoWriter('SideView3D.avi');
writerObj.FrameRate = 100;
%writerObj.FrameRate = 50;
open(writerObj);

M = 1;
strat = [-100:100/M:0]';
[X_sphere, Z_sphere, Y_sphere]=sphere;
X_sphere = X_sphere*R/(lam_f/1.5);
Y_sphere = Y_sphere*R/(g/(2*pi*f)^2)/1.5; % gamma = 1.047*gammaF;
Z_sphere = Z_sphere*R/lam_f;
x_plot = [x0 - 1.5:0.01:x0 + 2.5];
[X,Z] = meshgrid(x_plot,[-0.8 0.8]);

A_axis = 0.5*(2*pi*f)^2/g;


for t = dt:dt:4 % Discretization of intermediate projectile motion
   
    %%% Path of flight
    eta = y - G*t^2/2 + vy*t;
    Psi = A*sin(pi/2 + 2*pi*t)*sin(pi/2 + 2*pi*x) + gamma*sin(pi/2 + 2*pi*(2*t-varphi));
        
    if eta < Psi + R/(g/(2*pi*f)^2) % Just past impact
       
        %%% Bisection method to find time of intersection
        
        t_low = t;
        t_high = t-dt;
        
        for i = 1:10000
            
            t_impact = (t_low+t_high)/2;
            eta = y - G*t_impact^2/2 + vy*t_impact;
            Psi = A*sin(pi/2 + 2*pi*t_impact)*sin(pi/2 + 2*pi*x) + gamma*sin(pi/2 + 2*pi*(2*t_impact-varphi));
        
            if abs(eta - (Psi + R/(g/(2*pi*f)^2))) < 0.0001
                break
            else if eta < Psi + R/(g/(2*pi*f)^2)
                    t_low = t_impact;
                else
                    t_high = t_impact;
                end
            end
        end
               
        %%% Derivatives
    dPsidx = -A*2*pi*sin(pi/2 + 2*pi*t_impact)*sin(2*pi*x);
    dPsidt = -A*2*pi*sin(2*pi*t_impact)*sin(pi/2 + 2*pi*x) - 4*pi*gamma*sin(2*pi*(2*t_impact-varphi));
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
      
      %%%  Makes video
                       
            Psi_plot = A*sin(pi/2 + 2*pi*t_impact)*sin(pi/2 + 2*pi*X) + gamma*sin(pi/2 + 2*pi*(2*t_impact-varphi));
            blah = strat*ones(1,length(x_plot))+Psi_plot;
            S0 = surf(X_sphere+x, Z_sphere, Y_sphere-y+2*Psi,'edgecolor','none');
            hold on
            S1 = surf(X,Z,Psi_plot,'edgecolor','none');
            alpha(S1,0.9)
            hold on
            S2 = surf(X_sphere+x, Z_sphere, Y_sphere+y,'edgecolor','none');
            alpha(S2,0.5)
            hold on
            axis([x-1.5 x+2.5 -0.8 0.8 -5*A_axis 5*A_axis])
            colormap([0.9 0.9 0.1])
            view(0,60)
            shading interp
            lightangle(0,45)
            S1.FaceLighting = 'gouraud';
            S1.AmbientStrength = 0.3;
            S1.DiffuseStrength = 1;
            S1.SpecularStrength = 0.5;
            S1.SpecularExponent = 25;
            S1.BackFaceLighting = 'unlit';
            S2.FaceLighting = 'gouraud';
            S2.AmbientStrength = 0.3;
            S2.DiffuseStrength = 0.5;
            S2.SpecularStrength = 1;
            S2.SpecularExponent = 25;
            S2.SpecularColorReflectance = 1;
            S2.BackFaceLighting = 'unlit';
            h = fill3([X,fliplr(X)],[Z,fliplr(Z)],[blah(1,:),fliplr(blah(2,:))],[0.9,0.9,0], 'edgecolor','none');
            alpha(h,0.5)
            hold off
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
      
      break
      
    else
        
            if rem(t,0.004) == 0
            Psi_plot = A*sin(pi/2 + 2*pi*t)*sin(pi/2 + 2*pi*X) + gamma*sin(pi/2 + 2*pi*(2*t-varphi));
            blah = strat*ones(1,length(x_plot))+Psi_plot;
            S0 = surf(X_sphere+x, Z_sphere, Y_sphere-eta+2*Psi,'edgecolor','none');
            hold on
            S1 = surf(X,Z,Psi_plot,'edgecolor','none');
            alpha(S1,0.9)
            hold on
            S2 = surf(X_sphere+x, Z_sphere, Y_sphere+eta,'edgecolor','none');
            alpha(S2,0.5)
            hold on
            axis([x-1.5 x+2.5 -0.8 0.8 -5*A_axis 5*A_axis])
            colormap([0.9 0.9 0.1])
            view(0,60)
            shading interp
            lightangle(0,45)
            S1.FaceLighting = 'gouraud';
            S1.AmbientStrength = 0.3;
            S1.DiffuseStrength = 1;
            S1.SpecularStrength = 0.5;
            S1.SpecularExponent = 25;
            S1.BackFaceLighting = 'unlit';
            S2.FaceLighting = 'gouraud';
            S2.AmbientStrength = 0.3;
            S2.DiffuseStrength = 0.5;
            S2.SpecularStrength = 1;
            S2.SpecularExponent = 25;
            S2.SpecularColorReflectance = 1;
            S2.BackFaceLighting = 'unlit';
            h = fill3([X,fliplr(X)],[Z,fliplr(Z)],[blah(1,:),fliplr(blah(2,:))],[0.9,0.9,0], 'edgecolor','none');
            alpha(h,0.5)
            hold off
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
            end
      
    end
        
end

%%% Subsequent drops

Trapped = 0;
T=t_impact;

while T <= T_max
for t = dt:dt:4  % Discretization of intermediate projectile motion
   
    T = T + dt;
    
    %%% Path of flight
    xi = x + vx*t;
    eta = y - G*t^2/2 + vy*t;
    Psi = A*sin(pi/2 + 2*pi*T)*sin(pi/2 + 2*pi*xi) + gamma*sin(pi/2 + 2*pi*(2*T-varphi));
        
    if eta > Psi + R/(g/(2*pi*f)^2)
        Trapped = 0;
                    
            if rem(t,0.004) == 0
            Psi_plot = A*sin(pi/2 + 2*pi*T)*sin(pi/2 + 2*pi*X) + gamma*sin(pi/2 + 2*pi*(2*T-varphi));
            blah = strat*ones(1,length(x_plot))+Psi_plot;
            S0 = surf(X_sphere+xi, Z_sphere, Y_sphere-eta+2*Psi,'edgecolor','none');
            hold on
            S1 = surf(X,Z,Psi_plot,'edgecolor','none');
            alpha(S1,0.9)
            hold on
            S2 = surf(X_sphere+xi, Z_sphere, Y_sphere+eta,'edgecolor','none');
            alpha(S2,0.5)
            hold on
            axis([x0-1.5 x0+2.5 -0.8 0.8 -5*A_axis 5*A_axis])
            colormap([0.9 0.9 0.1])
            view(0,60)
            shading interp
            lightangle(0,45)
            S1.FaceLighting = 'gouraud';
            S1.AmbientStrength = 0.3;
            S1.DiffuseStrength = 1;
            S1.SpecularStrength = 0.5;
            S1.SpecularExponent = 25;
            S1.BackFaceLighting = 'unlit';
            S2.FaceLighting = 'gouraud';
            S2.AmbientStrength = 0.3;
            S2.DiffuseStrength = 0.5;
            S2.SpecularStrength = 1;
            S2.SpecularExponent = 25;
            S2.SpecularColorReflectance = 1;
            S2.BackFaceLighting = 'unlit';
            h = fill3([X,fliplr(X)],[Z,fliplr(Z)],[blah(1,:),fliplr(blah(2,:))],[0.9,0.9,0], 'edgecolor','none');
            alpha(h,0.5)
            hold off
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
            end
    
    else %%% Just past impact
        
        if Trapped == 1
            
            %%% Chatering fix
            dxidt = vx;
            dPsidx = -2*pi*A*sin(pi/2 + 2*pi*T)*sin(2*pi*xi);
            dPsidt = -2*pi*A*sin(2*pi*T)*sin(pi/2 + 2*pi*xi) - 4*pi*gamma*sin(2*pi*(2*T-varphi)) + dPsidx*dxidt;

            
            y = Psi + R/(g/(2*pi*f)^2);
            vy = dPsidt;
                        
            if rem(t,0.004) == 0
            Psi_plot = A*sin(pi/2 + 2*pi*T)*sin(pi/2 + 2*pi*X) + gamma*sin(pi/2 + 2*pi*(2*T-varphi));
            blah = strat*ones(1,length(x_plot))+Psi_plot;
            S0 = surf(X_sphere+xi, Z_sphere, Y_sphere-y+2*Psi,'edgecolor','none');
            hold on
            S1 = surf(X,Z,Psi_plot,'edgecolor','none');
            alpha(S1,0.9)
            hold on
            S2 = surf(X_sphere+xi, Z_sphere, Y_sphere+y,'edgecolor','none');
            alpha(S2,0.5)
            hold on
            axis([x0-1.5 x0+2.5 -0.8 0.8 -5*A_axis 5*A_axis])
            colormap([0.9 0.9 0.1])
            view(0,60)
            shading interp
            lightangle(0,45)
            S1.FaceLighting = 'gouraud';
            S1.AmbientStrength = 0.3;
            S1.DiffuseStrength = 1;
            S1.SpecularStrength = 0.5;
            S1.SpecularExponent = 25;
            S1.BackFaceLighting = 'unlit';
            S2.FaceLighting = 'gouraud';
            S2.AmbientStrength = 0.3;
            S2.DiffuseStrength = 0.5;
            S2.SpecularStrength = 1;
            S2.SpecularExponent = 25;
            S2.SpecularColorReflectance = 1;
            S2.BackFaceLighting = 'unlit';
            h = fill3([X,fliplr(X)],[Z,fliplr(Z)],[blah(1,:),fliplr(blah(2,:))],[0.9,0.9,0], 'edgecolor','none');
            alpha(h,0.5)
            hold off
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
            end
            
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
            xi = x + vx*t_impact;
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
            
            x_plot = [x0 - 1:0.01:x0 + 3];
            
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
            vx = u_out*sin(pi/2 + theta_out)*Tf/lam_f;
            if abs(vx) < 1e-15
               vx = 0; 
            end
            vy = u_out*sin(theta_out)*Tf/(g/(2*pi*f)^2);
            x = xi;
            y = Psi + R/(g/(2*pi*f)^2);
            
            FA = 0.45;
            if vy < dPsidt
                %vy = dPsidt;
                Trapped = 1;
            end
                        
            Psi_plot = A*sin(pi/2 + 2*pi*T)*sin(pi/2 + 2*pi*X) + gamma*sin(pi/2 + 2*pi*(2*T-varphi));
            blah = strat*ones(1,length(x_plot))+Psi_plot;
            S0 = surf(X_sphere+x, Z_sphere, Y_sphere-y+2*Psi,'edgecolor','none');
            hold on
            S1 = surf(X,Z,Psi_plot,'edgecolor','none');
            alpha(S1,0.9)
            hold on
            S2 = surf(X_sphere+x, Z_sphere, Y_sphere+y,'edgecolor','none');
            alpha(S2,0.5)
            hold on
            axis([x0-1.5 x0+2.5 -0.8 0.8 -5*A_axis 5*A_axis])
            colormap([0.9 0.9 0.1])
            view(0,60)
            shading interp
            lightangle(0,45)
            S1.FaceLighting = 'gouraud';
            S1.AmbientStrength = 0.3;
            S1.DiffuseStrength = 1;
            S1.SpecularStrength = 0.5;
            S1.SpecularExponent = 25;
            S1.BackFaceLighting = 'unlit';
            S2.FaceLighting = 'gouraud';
            S2.AmbientStrength = 0.3;
            S2.DiffuseStrength = 0.5;
            S2.SpecularStrength = 1;
            S2.SpecularExponent = 25;
            S2.SpecularColorReflectance = 1;
            S2.BackFaceLighting = 'unlit';
            h = fill3([X,fliplr(X)],[Z,fliplr(Z)],[blah(1,:),fliplr(blah(2,:))],[0.9,0.9,0], 'edgecolor','none');
            alpha(h,0.5)
            hold off
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
           
                        break
            
        end
            
    end
    
end
end

close(writerObj);
toc