function sol = ell_shearcrack(t,theta,phi,S,N,taus)
% sol = ell_shearcrack(t,theta,phi,S,N,taus)
%
% This function computes the normalised 3 components of acceleration and
% stress rates produced by a dynamically growing elliptical crack located
% in the plane (x,y) with center (0,0,0), for a stress drop in tau_zx,
% at an observation point given by angular positions (theta,phi) from the
% crack center.
%
% The solution is given in Richards, IJSS, 1973 (corrected from a few
% typos).
%
% accelerations are normalised by: (8 cs^2 S N a)/(pi cd^3 R)
% stress rates are normalised by: (16 mu cs^2 S N a)/(pi cd^4 R)
% time is normalised by: cd/R
%
% where: cs     is shear wave speed
%        cd     is P wave speed
%        S  is rupture speed along the x-direction
%        N     is rupture speed along the y-direction
%        a      is the slip rate at the crack center
%        R      is the distance between observation point and crack center
%        mu     is shear modulus
%
% input:
%   t:      an array of times (typically linspace(0,tmax))
%   theta:  angle between normal to crack (z-axis) and observation point
%   phi:    third polar cordinate (i.e., the other angle ! phi=0 is along x-coordinate)
%   S:      rupture speed along x-axis normalised by P wave speed
%   N:      rupture speed along y-axis normalised by P wave speed
%   taus:   Vp/Vs ratio
%
% output:
%   sol:    a structure with fields:
%       t:     array of time (same as input "t")
%       theta: same as input
%       phi:   same as input
%       S:     same as input
%       N:     same as input
%       taus:  same as input
%       tR:    time when rupture passes nearby
%       uj:    acceleration along axis j (j=x,y,z)
%       zj:    stress rate component tau_{zj} (j=x,y,z)
%
% EXAMPLE USE:
% 
% t = linspace(0,3,500);
% theta = pi/2 + 0.02; %(a point very close to the rupture plane)
% phi = -pi/48; %(mixed mode II/III)
% S = 0.6;
% N = 0.8;
% taus = 0.8/0.6;
% sol = ell_shearcrack(t,theta,phi,S,N,taus)
% 
% plot(sol.t,sol.ux,sol.t,sol.uy,sol.t,sol.uz,sol.tR,0,'o');
% legend('{\itx} component','{\ity} component','{\itz} component','rupture arrival');
% xlabel('time {\itt}/({\itR}/{\itc}_d)');
% ylabel('acceleration {\itu}_{\ittt}/((8 {\itc}_s^2 {\itS N a})/(\pi {\itc}_d^3 {\itR}))');
% xlim([0.5 2.5]);
% title('\phi=-\pi/48,  \theta=\pi/2 + 0.02,  S=0.6,  N=0.8');
%
% CODED BY N. BRANTUT (nicolas.brantut@normalesup.org). Last edited: 5 Nov. 2015.

%% Some definitions (Notations same as in Richards, 1973)

D = S.^2 .*(cos(phi)).^2 + N.^2 .*(sin(phi)).^2;
F = S.^2 .*(sin(phi)).^2 + N.^2 .*(cos(phi)).^2;

Delta = (S.^2 - N.^2)*cos(phi)*sin(phi);

G = D.*(sin(theta)).^2 .*(F + D.*(cos(theta).^2)) - S.^2.*N.^2;

%critical values of w on the Cagniard path in complex w-plane
w0d = (D.*(cos(theta)).^2 .*(1-D.*(sin(theta).^2))./G)^(1/2);
w0s = (D.*(cos(theta)).^2 .*(1-D.*taus.^2.*(sin(theta).^2))./G)^(1/2);

%corresponding times
tau0d = (w0d.^2 .*S.^2.*N.^2 + D).^(1/2) ./(sin(theta)*D);
tau0s = (w0s.^2 .*S.^2.*N.^2 + D).^(1/2) ./(sin(theta)*D);

% Singularity of the integrand
q_SN = @(w) (w.*Delta + 1i*(w.^2.*S^2*N^2 + D).^(1/2))./D;

% Cagniard paths
q_d = @(tau,w) (tau.^2 - w.^2 -1).^(1/2).*cos(theta) + 1i.*tau.*sin(theta);
q_s = @(tau,w) (tau.^2 - w.^2 -1).^(1/2).*cos(theta) + 1i.*tau.*sin(theta);

%asymptotics for wd(tau) (to check)
wd_asymp = @(tau) tau*D./( (D^2+(Delta+1i*S*N)^2)^(1/2)*cos(theta) - 1i*(Delta+1i*S*N)*sin(theta));

%% Set the solution structure with the input parameters (for future use)

sol.t = t;
sol.theta = theta;
sol.phi = phi;
sol.S = S;
sol.N = N;
sol.taus = taus;

% get the time at which rupture tip is close to observation point
% i.e., when tip passes at the position which is the projection of the 
% observation point on the (xy) plane
x = sin(theta)*cos(phi);
y = sin(theta)*sin(phi);

sol.tR = sqrt((x/S)^2 + (y/N)^2);

%% computation of the solution

%initialise wd and ws
wd = 1+1i;
ws = 1+1i;

%loop for all times

for k = 1:length(t)
    
    % normalised time limits for integration
    Td = (t(k).^2 - 1).^(1/2);
    Ts = (t(k).^2 - taus^2).^(1/2);
    
    % initialise values of each component
    Id_ux = 0;
    Id_uy = 0;
    Id_uz = 0;
    Id_zx = 0;
    Id_zy = 0;
    Id_zz = 0;
    
    Is_ux = 0;
    Is_uy = 0;
    Is_uz = 0;
    Is_zx = 0;
    Is_zy = 0;
    Is_zz = 0;
    
    Rd_ux = 0;
    Rd_uy = 0;
    Rd_uz = 0;
    Rd_zx = 0;
    Rd_zy = 0;
    Rd_zz = 0;
    
    Rs_ux = 0;
    Rs_uy = 0;
    Rs_uz = 0;
    Rs_zx = 0;
    Rs_zy = 0;
    Rs_zz = 0;
    
    if G>0
        if t(k)>1
            if t(k)>tau0d
            
                % find the complex w_d value (cagniard path)
                [wd] = findwq(t(k),'d',wd);
                
                % compute the corresponding q_snd and q_d values
                qsnd = q_SN(wd);
                qd   = q_d(t(k),wd);
                

                %compute the residual terms
                Rd_ux = Rd_ux_calc(wd,qsnd);
                Rd_uy = Rd_uy_calc(wd,qsnd);
                Rd_uz = Rd_uz_calc(wd,qsnd);
                Rd_zx = Rd_zx_calc(wd,qsnd);
                Rd_zy = Rd_zy_calc(wd,qsnd);
                Rd_zz = Rd_zz_calc(wd,qsnd);
                
                %compute the integral terms
                Id_ux = quadgk(@(w) integrand_Id_ux(w, t(k)) - Sd_ux(w,t(k),wd,qd), 0,Td) + int_Sd_ux(Td,t(k),wd,qd);
                Id_uy = quadgk(@(w) integrand_Id_uy(w, t(k)) - Sd_uy(w,t(k),wd,qd), 0,Td) + int_Sd_uy(Td,t(k),wd,qd);
                Id_uz = quadgk(@(w) integrand_Id_uz(w, t(k)) - Sd_uz(w,t(k),wd,qd), 0,Td) + int_Sd_uz(Td,t(k),wd,qd);
                Id_zx = quadgk(@(w) integrand_Id_zx(w, t(k)) - Sd_zx(w,t(k),wd,qd), 0,Td) + int_Sd_zx(Td,t(k),wd,qd);
                Id_zy = quadgk(@(w) integrand_Id_zy(w, t(k)) - Sd_zy(w,t(k),wd,qd), 0,Td) + int_Sd_zy(Td,t(k),wd,qd);
                Id_zz = quadgk(@(w) integrand_Id_zz(w, t(k)) - Sd_zz(w,t(k),wd,qd), 0,Td) + int_Sd_zz(Td,t(k),wd,qd);

            else
                %other wise just compute the integral terms (if integrand
                %is not singular)
                Id_ux = quadgk(@(w) integrand_Id_ux(w, t(k)), 0,Td);
                Id_uy = quadgk(@(w) integrand_Id_uy(w, t(k)), 0,Td);
                Id_uz = quadgk(@(w) integrand_Id_uz(w, t(k)), 0,Td);
                Id_zx = quadgk(@(w) integrand_Id_zx(w, t(k)), 0,Td);
                Id_zy = quadgk(@(w) integrand_Id_zy(w, t(k)), 0,Td);
                Id_zz = quadgk(@(w) integrand_Id_zz(w, t(k)), 0,Td);
            end
        end
    
        if t(k)>taus
            if t(k)>tau0s
            
                % same shit for the S components....
                [ws] = findwq(t(k),'s',ws);
                
                qsns = q_SN(ws);
                qs   = q_s(t(k),ws);

                Rs_ux = Rs_ux_calc(ws,qsns);
                Rs_uy = Rs_uy_calc(ws,qsns);
                Rs_uz = Rs_uz_calc(ws,qsns);
                Rs_zx = Rs_zx_calc(ws,qsns);
                Rs_zy = Rs_zy_calc(ws,qsns);
                Rs_zz = Rs_zz_calc(ws,qsns);
                
                Is_ux = quadgk(@(w) integrand_Is_ux(w, t(k)) - Ss_ux(w,t(k),ws,qs), 0,Ts) + int_Ss_ux(Ts,t(k),ws,qs);
                Is_uy = quadgk(@(w) integrand_Is_uy(w, t(k)) - Ss_uy(w,t(k),ws,qs), 0,Ts) + int_Ss_uy(Ts,t(k),ws,qs);
                Is_uz = quadgk(@(w) integrand_Is_uz(w, t(k)) - Ss_uz(w,t(k),ws,qs), 0,Ts) + int_Ss_uz(Ts,t(k),ws,qs);
                Is_zx = quadgk(@(w) integrand_Is_zx(w, t(k)) - Ss_zx(w,t(k),ws,qs), 0,Ts) + int_Ss_zx(Ts,t(k),ws,qs);
                Is_zy = quadgk(@(w) integrand_Is_zy(w, t(k)) - Ss_zy(w,t(k),ws,qs), 0,Ts) + int_Ss_zy(Ts,t(k),ws,qs);
                Is_zz = quadgk(@(w) integrand_Is_zz(w, t(k)) - Ss_zz(w,t(k),ws,qs), 0,Ts) + int_Ss_zz(Ts,t(k),ws,qs);

            else
                Is_ux = quadgk(@(w) integrand_Is_ux(w, t(k)), 0,Ts);
                Is_uy = quadgk(@(w) integrand_Is_uy(w, t(k)), 0,Ts);
                Is_uz = quadgk(@(w) integrand_Is_uz(w, t(k)), 0,Ts);
                Is_zx = quadgk(@(w) integrand_Is_zx(w, t(k)), 0,Ts);
                Is_zy = quadgk(@(w) integrand_Is_zy(w, t(k)), 0,Ts);
                Is_zz = quadgk(@(w) integrand_Is_zz(w, t(k)), 0,Ts);
            end
        end
    else

        if t(k)>1
            Id_ux = quadgk(@(w) integrand_Id_ux(w, t(k)), 0,Td);
            Id_uy = quadgk(@(w) integrand_Id_uy(w, t(k)), 0,Td);
            Id_uz = quadgk(@(w) integrand_Id_uz(w, t(k)), 0,Td);
            Id_zx = quadgk(@(w) integrand_Id_zx(w, t(k)), 0,Td);
            Id_zy = quadgk(@(w) integrand_Id_zy(w, t(k)), 0,Td);
            Id_zz = quadgk(@(w) integrand_Id_zz(w, t(k)), 0,Td);
        end
    
        if t(k)>taus
            Is_ux = quadgk(@(w) integrand_Is_ux(w, t(k)), 0,Ts);
            Is_uy = quadgk(@(w) integrand_Is_uy(w, t(k)), 0,Ts);
            Is_uz = quadgk(@(w) integrand_Is_uz(w, t(k)), 0,Ts);
            Is_zx = quadgk(@(w) integrand_Is_zx(w, t(k)), 0,Ts);
            Is_zy = quadgk(@(w) integrand_Is_zy(w, t(k)), 0,Ts);
            Is_zz = quadgk(@(w) integrand_Is_zz(w, t(k)), 0,Ts);
        end
    end
    

    %% assign solution to structure
    sol.ux(k) = Id_ux + Is_ux + Rd_ux + Rs_ux;
    sol.uy(k) = Id_uy + Is_uy + Rd_uy + Rs_uy;
    sol.uz(k) = Id_uz + Is_uz + Rd_uz + Rs_uz;
    sol.zx(k) = Id_zx + Is_zx + Rd_zx + Rs_zx;
    sol.zy(k) = Id_zy + Is_zy + Rd_zy + Rs_zy;
    sol.zz(k) = Id_zz + Is_zz + Rd_zz + Rs_zz;
    
  
end

%% Function to compute w-plane Cagniard path
    function [wa] = findwq(tau,alpha,w0)
        if strcmp(alpha,'d')
            taua = 1;
        elseif strcmp(alpha,'s')
            taua = taus;
        else
            disp('Wrong wave type in findwq');
            return
        end
        
        Z      = @(w) -1i.*q_SN(w).*sin(theta) + (q_SN(w).^2 + w.^2 + taua^2).^(1/2).*cos(theta);
        dqsndw = @(w) (Delta + 1i*w.*S^2*N^2*(w.^2*S^2*N^2 + D).^(-1/2))/D;
        dZdw   = @(w) -1i*sin(theta)*dqsndw(w) + cos(theta)*(q_SN(w).^2 + w.^2 + taua^2).^(-1/2).*(q_SN(w).*dqsndw(w) + w);

        %find wa(tau)
        tol = 1e-8;
        imax = 200;
        counter = 0;
        wa = w0;
        residual = tau - Z(wa);
        while abs(residual)>tol
            wa = wa + residual/dZdw(wa);
            residual = tau - Z(wa);
            counter = counter + 1;
            if counter >imax
                disp('Max. Nmber of iteration reached !')
                return
            end
        end
        
        %{
        rew = linspace(2,5);
        imw = linspace(-.5,.5);
        for k1=1:length(rew)
            for k2=1:length(imw)
                check(k1,k2) = tau - Z(complex(rew(k1),imw(k2)));
            end
        end
        
        pcolor(rew,imw,abs(check'));
        axis equal
        hold on;
        plot(wa,'o');
        plot(w0s,0,'ks');
        plot(w0d,0,'rs');
        pause;
        hold off
        
        %}
        
        %wa = complex(abs(real(wa)),abs(imag(wa)));
        %disp(abs(tau-Z(wa)));
        
    end

%% integral terms for ux

    function res = integrand_Id_ux(w,tau)
        
        tauwd = 1.*(w.^2 + 1).^(1/2);
        q = (tau.^2 - tauwd.^2).^(1/2).*cos(theta) + 1i.*tau.*sin(theta);
        E = 1 + q.^2.*D + w.^2.*F;
        O = -2.*q.*w.*(S.^2 - N.^2).*cos(phi).*sin(phi);
        dqddt = (cos(theta).*tau.*(tau.^2-tauwd.^2).^(-1/2) + 1i.*sin(theta));
        
        Md = -q.^2.*cos(phi).^2 - w.^2.*sin(phi).^2;
        Nd = 2.*q.*w.*cos(phi).*sin(phi);
        
        res = real( ((E.^2 + O.^2).*Md - 2*E.*O.*Nd)./((E.^2 - O.^2).^2 ) .* dqddt);
        
    end

    function res = Sd_ux(w,tau,wa,qa)
        alpha = 'd';
        
        Va  = -a4(wa,qa)^2;
        Wap = -2*a4(wa,qa)*a8(a3(wa,alpha,tau));
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa./(w - wa).^2 + Ya./(w - wa));
    
    end

    function res = int_Sd_ux(Ta,tau,wa,qa)
        alpha = 'd';
        
        Va  = -a4(wa,qa)^2;
        Wap = -2*a4(wa,qa)*a8(a3(wa,alpha,tau));
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa*Ta./(wa*(wa - Ta)) + Ya*(log(abs(1-Ta/wa)) + 1i*angle(1-Ta/wa) ));
    
    end

    function res = integrand_Is_ux(w,tau)
        
        tauws = 1.*(w.^2 + taus.^2).^(1/2);
        q = (tau.^2 - tauws.^2).^(1/2).*cos(theta) + 1i.*tau.*sin(theta);
        E = 1 + q.^2.*D + w.^2.*F;
        O = -2.*q.*w.*(S.^2 - N.^2).*cos(phi).*sin(phi);
        dqsdt = (cos(theta).*tau.*(tau.^2-tauws.^2).^(-1/2) + 1i.*sin(theta));
        
        Ms = q.^2.*cos(phi).^2 + w.^2.*sin(phi).^2 + .5.*taus.^2;
        Ns = -2.*q.*w.*cos(phi).*sin(phi);
        
        res = real( ((E.^2 + O.^2).*Ms - 2*E.*O.*Ns)./((E.^2 - O.^2).^2 ) .* dqsdt);
        
    end

    function res = Ss_ux(w,tau,wa,qa)
        alpha = 's';
        
        Va  = a4(wa,qa)^2 + .5*taus^2;
        Wap = 2*a4(wa,qa)*a8(a3(wa,alpha,tau));
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa./(w - wa).^2 + Ya./(w - wa));
    
    end

    function res = int_Ss_ux(Ta,tau,wa,qa)
        alpha = 's';
        
        Va  = a4(wa,qa)^2 + .5*taus^2;
        Wap = 2*a4(wa,qa)*a8(a3(wa,alpha,tau));
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa*Ta./(wa*(wa - Ta)) + Ya*(log(abs(1-Ta/wa)) + 1i*angle(1-Ta/wa)));
    
    end

%% integral terms for uy

    function res = integrand_Id_uy(w,tau)
        
        tauwd = 1.*(w.^2 + 1).^(1/2);
        q = (tau.^2 - tauwd.^2).^(1/2).*cos(theta) + 1i.*tau.*sin(theta);
        E = 1 + q.^2.*D + w.^2.*F;
        O = -2.*q.*w.*(S.^2 - N.^2).*cos(phi).*sin(phi);
        dqddt = (cos(theta).*tau.*(tau.^2-tauwd.^2).^(-1/2) + 1i.*sin(theta));
        
        Md = (-q.^2 + w.^2).*cos(phi).*sin(phi);
        Nd = -q.*w.*(cos(phi).^2 - sin(phi).^2);
        
        res = real( ((E.^2 + O.^2).*Md - 2*E.*O.*Nd)./((E.^2 - O.^2).^2 ) .* dqddt);
        
    end

    function res = Sd_uy(w,tau,wa,qa)
        alpha = 'd';
        
        Va  = -a5(wa,qa);
        Wap = -a9(wa,qa,a3(wa,alpha,tau));
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa./(w - wa).^2 + Ya./(w - wa));
    
    end

    function res = int_Sd_uy(Ta,tau,wa,qa)
        alpha = 'd';
        
        Va  = -a5(wa,qa);
        Wap = -a9(wa,qa,a3(wa,alpha,tau));
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa*Ta./(wa*(wa - Ta)) + Ya*(log(abs(1-Ta/wa)) + 1i*angle(1-Ta/wa) ));
    
    end

    function res = integrand_Is_uy(w,tau)
        
        tauws = 1.*(w.^2 + taus.^2).^(1/2);
        q = (tau.^2 - tauws.^2).^(1/2).*cos(theta) + 1i.*tau.*sin(theta);
        E = 1 + q.^2.*D + w.^2.*F;
        O = -2.*q.*w.*(S.^2 - N.^2).*cos(phi).*sin(phi);
        dqsdt = (cos(theta).*tau.*(tau.^2-tauws.^2).^(-1/2) + 1i.*sin(theta));
        
        Ns = q.*w.*(cos(phi).^2 - sin(phi).^2);
        Ms = (q.^2-w.^2).*cos(phi).*sin(phi);
        
        res = real( ((E.^2 + O.^2).*Ms - 2*E.*O.*Ns)./((E.^2 - O.^2).^2 ) .* dqsdt);
        
    end

    function res = Ss_uy(w,tau,wa,qa)
        alpha = 's';
        
        Va  = a5(wa,qa);
        Wap = a9(wa,qa,a3(wa,alpha,tau));
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa./(w - wa).^2 + Ya./(w - wa));
    
    end

    function res = int_Ss_uy(Ta,tau,wa,qa)
        alpha = 's';
        
        Va  = a5(wa,qa);
        Wap = a9(wa,qa,a3(wa,alpha,tau));
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa*Ta./(wa*(wa - Ta)) + Ya*(log(abs(1-Ta/wa)) + 1i*angle(1-Ta/wa)));
    
    end

%% integral terms for uz

    function res = integrand_Id_uz(w,tau)
        
        tauwd = 1.*(w.^2 + 1).^(1/2);
        q = (tau.^2 - tauwd.^2).^(1/2).*cos(theta) + 1i.*tau.*sin(theta);
        md = (q.^2+w.^2+1).^(1/2);
        E = 1 + q.^2.*D + w.^2.*F;
        O = -2.*q.*w.*(S.^2 - N.^2).*cos(phi).*sin(phi);
        dqddt = (cos(theta).*tau.*(tau.^2-tauwd.^2).^(-1/2) + 1i.*sin(theta));
        
        Md = -1i*q.*md.*cos(phi);
        Nd = 1i*w.*md.*sin(phi);
        
        res = real( ((E.^2 + O.^2).*Md - 2*E.*O.*Nd)./((E.^2 - O.^2).^2 ) .* dqddt);
    end

    function res = Sd_uz(w,tau,wa,qa)
        alpha = 'd';
        
        md = (qa^2 + wa^2 + 1)^(1/2);
        aa3 = a3(wa,alpha,tau);
        
        Va  = -1i*a4(wa,qa)^2 * md;
        Wap = -1i*(a8(aa3)*md + a4(wa,qa)*a10(wa,qa,aa3)/md);
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa./(w - wa).^2 + Ya./(w - wa));
    
    end

    function res = int_Sd_uz(Ta,tau,wa,qa)
        alpha = 'd';
        
        md = (qa^2 + wa^2 + 1)^(1/2);
        aa3 = a3(wa,alpha,tau);
        
        Va  = -1i*a4(wa,qa)^2 * md;
        Wap = -1i*(a8(aa3)*md + a4(wa,qa)*a10(wa,qa,aa3)/md);
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa*Ta./(wa*(wa - Ta)) + Ya*(log(abs(1-Ta/wa)) + 1i*angle(1-Ta/wa)));
    
    end

    function res = integrand_Is_uz(w,tau)
        
        tauws = 1.*(w.^2 + taus.^2).^(1/2);
        q = (tau.^2 - tauws.^2).^(1/2).*cos(theta) + 1i.*tau.*sin(theta);
        ms = (q.^2 + w.^2+ taus^2).^(1/2);
        E = 1 + q.^2.*D + w.^2.*F;
        O = -2.*q.*w.*(S.^2 - N.^2).*cos(phi).*sin(phi);
        dqsdt = (cos(theta).*tau.*(tau.^2-tauws.^2).^(-1/2) + 1i.*sin(theta));
        
        Ms = 1i*q.*(ms - .5*taus^2./ms)*cos(phi);
        Ns = -1i*w.*(ms - .5*taus^2./ms)*sin(phi);
        
        res = real( ((E.^2 + O.^2).*Ms - 2*E.*O.*Ns)./((E.^2 - O.^2).^2 ) .* dqsdt);
    end

    function res = Ss_uz(w,tau,wa,qa)
        alpha = 's';
        
        aa3 = a3(wa,alpha,tau);
        
        Va  = 1i*a4(wa,qa)*a71(wa,qa);
        Wap = 1i*(a8(aa3)*a71(wa,qa) + a4(wa,qa)*a10(wa,qa,aa3)*a73(wa,qa)/qa);
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa./(w - wa).^2 + Ya./(w - wa));
    
    end

    function res = int_Ss_uz(Ta,tau,wa,qa)
        alpha = 's';
        
        aa3 = a3(wa,alpha,tau);
        
        Va  = 1i*a4(wa,qa)*a71(wa,qa);
        Wap = 1i*(a8(aa3)*a71(wa,qa) + a4(wa,qa)*a10(wa,qa,aa3)*a73(wa,qa)/qa);
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa*Ta./(wa*(wa - Ta)) + Ya*(log(abs(1-Ta/wa)) + 1i*angle(1-Ta/wa)));
    
    end

%% integral terms for tau_zx

    function res = integrand_Id_zx(w,tau)
        
        tauwd = 1.*(w.^2 + 1).^(1/2);
        q = (tau.^2 - tauwd.^2).^(1/2).*cos(theta) + 1i.*tau.*sin(theta);
        md = (q.^2+w.^2+1).^(1/2);
        E = 1 + q.^2.*D + w.^2.*F;
        O = -2.*q.*w.*(S.^2 - N.^2).*cos(phi).*sin(phi);
        dqddt = (cos(theta).*tau.*(tau.^2-tauwd.^2).^(-1/2) + 1i.*sin(theta));
        
        Md = md.*(q.^2*cos(phi)^2 + w.^2*sin(phi));
        Nd = -2*q.*w.*md.*cos(phi)*sin(phi);
        
        res = real( ((E.^2 + O.^2).*Md - 2*E.*O.*Nd)./((E.^2 - O.^2).^2 ) .* dqddt);
    end

    function res = Sd_zx(w,tau,wa,qa)
        alpha = 'd';
        
        md = (qa^2 + wa^2 + 1)^(1/2);
        aa3 = a3(wa,alpha,tau);
        
        Va  = a4(wa,qa)^2*md;
        Wap = 2*a4(wa,qa)*a8(aa3)*md + a4(wa,qa)^2*a10(wa,qa,aa3)/md;
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa./(w - wa).^2 + Ya./(w - wa));
    
    end

    function res = int_Sd_zx(Ta,tau,wa,qa)
        alpha = 'd';
        
        md = (qa^2 + wa^2 + 1)^(1/2);
        aa3 = a3(wa,alpha,tau);
        
        Va  = a4(wa,qa)^2*md;
        Wap = 2*a4(wa,qa)*a8(aa3)*md + a4(wa,qa)^2*a10(wa,qa,aa3)/md;
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa*Ta./(wa*(wa - Ta)) + Ya*(log(abs(1-Ta/wa)) + 1i*angle(1-Ta/wa)));
    
    end

    function res = integrand_Is_zx(w,tau)
        
        tauws = 1.*(w.^2 + taus.^2).^(1/2);
        q = (tau.^2 - tauws.^2).^(1/2).*cos(theta) + 1i.*tau.*sin(theta);
        ms = (q.^2 + w.^2+ taus^2).^(1/2);
        E = 1 + q.^2.*D + w.^2.*F;
        O = -2.*q.*w.*(S.^2 - N.^2).*cos(phi).*sin(phi);
        dqsdt = (cos(theta).*tau.*(tau.^2-tauws.^2).^(-1/2) + 1i.*sin(theta));
        
        Ms = (-q.^2*cos(phi)^2 - w.^2*sin(phi)^2).*(ms - .25*taus^2./ms)-.25*taus^2.*ms;
        Ns = 2*q.*w.*(ms - .25*taus^2./ms)*cos(phi)*sin(phi);
        
        res = real( ((E.^2 + O.^2).*Ms - 2*E.*O.*Ns)./((E.^2 - O.^2).^2 ) .* dqsdt);
    end

    function res = Ss_zx(w,tau,wa,qa)
        alpha = 's';
        
        aa3 = a3(wa,alpha,tau);
        
        Va  = -a4(wa,qa)^2*a72(wa,qa) - .25*taus^2*a7(wa,qa);
        Wap = -2*a4(wa,qa)*a8(aa3)*a72(wa,qa) - (a4(wa,qa)^2*a74(wa,qa)/qa + .25*taus^2/a7(wa,qa))*a10(wa,qa,aa3);
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa./(w - wa).^2 + Ya./(w - wa));
    
    end

    function res = int_Ss_zx(Ta,tau,wa,qa)
        alpha = 's';
        
        aa3 = a3(wa,alpha,tau);
        
        Va  = -a4(wa,qa)^2*a72(wa,qa) - .25*taus^2*a7(wa,qa);
        Wap = -2*a4(wa,qa)*a8(aa3)*a72(wa,qa) - (a4(wa,qa)^2*a74(wa,qa)/qa + .25*taus^2/a7(wa,qa))*a10(wa,qa,aa3);
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa*Ta./(wa*(wa - Ta)) + Ya*(log(abs(1-Ta/wa)) + 1i*angle(1-Ta/wa)));
    
    end

%% integral terms for tau_zy

    function res = integrand_Id_zy(w,tau)
        
        tauwd = 1.*(w.^2 + 1).^(1/2);
        q = (tau.^2 - tauwd.^2).^(1/2).*cos(theta) + 1i.*tau.*sin(theta);
        md = (q.^2+w.^2+1).^(1/2);
        E = 1 + q.^2.*D + w.^2.*F;
        O = -2.*q.*w.*(S.^2 - N.^2).*cos(phi).*sin(phi);
        dqddt = (cos(theta).*tau.*(tau.^2-tauwd.^2).^(-1/2) + 1i.*sin(theta));
        
        Md = md.*(q.^2-w.^2)*cos(phi)*sin(phi);
        Nd = q.*w.*md.*(cos(phi)^2 - sin(phi)^2);
        
        res = real( ((E.^2 + O.^2).*Md - 2*E.*O.*Nd)./((E.^2 - O.^2).^2 ) .* dqddt);
    end

    function res = Sd_zy(w,tau,wa,qa)
        alpha = 'd';
        
        md = (qa^2 + wa^2 + 1)^(1/2);
        aa3 = a3(wa,alpha,tau);
        
        Va  = a5(wa,qa)*md;
        Wap = a9(wa,qa,aa3)*md + a5(wa,qa)*a10(wa,qa,aa3)/md;
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa./(w - wa).^2 + Ya./(w - wa));
    
    end

    function res = int_Sd_zy(Ta,tau,wa,qa)
        alpha = 'd';
        
        md = (qa^2 + wa^2 + 1)^(1/2);
        aa3 = a3(wa,alpha,tau);
        
        Va  = a5(wa,qa)*md;
        Wap = a9(wa,qa,aa3)*md + a5(wa,qa)*a10(wa,qa,aa3)/md;
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa*Ta./(wa*(wa - Ta)) + Ya*(log(abs(1-Ta/wa)) + 1i*angle(1-Ta/wa)));
    
    end

    function res = integrand_Is_zy(w,tau)
        
        tauws = 1.*(w.^2 + taus.^2).^(1/2);
        q = (tau.^2 - tauws.^2).^(1/2).*cos(theta) + 1i.*tau.*sin(theta);
        ms = (q.^2 + w.^2+ taus^2).^(1/2);
        E = 1 + q.^2.*D + w.^2.*F;
        O = -2.*q.*w.*(S.^2 - N.^2).*cos(phi).*sin(phi);
        dqsdt = (cos(theta).*tau.*(tau.^2-tauws.^2).^(-1/2) + 1i.*sin(theta));
        
        Ms = (-q.^2+w.^2).*(ms - .25*taus^2./ms)*cos(phi)*sin(phi);
        Ns = -q.*w.*(ms - .25*taus^2./ms).*(cos(phi)^2 - sin(phi)^2);
        
        res = real( ((E.^2 + O.^2).*Ms - 2*E.*O.*Ns)./((E.^2 - O.^2).^2 ) .* dqsdt);
    end

    function res = Ss_zy(w,tau,wa,qa)
        alpha = 's';
        
        aa3 = a3(wa,alpha,tau);
        
        Va  = -a5(wa,qa)*a72(wa,qa);
        Wap = -a9(wa,qa,aa3)*a72(wa,qa) - a5(wa,qa)*a10(wa,qa,aa3)*a74(wa,qa)/qa;
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa./(w - wa).^2 + Ya./(w - wa));
    
    end

    function res = int_Ss_zy(Ta,tau,wa,qa)
        alpha = 's';
        
        aa3 = a3(wa,alpha,tau);
        
        Va  = -a5(wa,qa)*a72(wa,qa);
        Wap = -a9(wa,qa,aa3)*a72(wa,qa) - a5(wa,qa)*a10(wa,qa,aa3)*a74(wa,qa)/qa;
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa*Ta./(wa*(wa - Ta)) + Ya*(log(abs(1-Ta/wa)) + 1i*angle(1-Ta/wa)));
    
    end

%% integral terms for tau_zz

    function res = integrand_Id_zz(w,tau)
        
        tauwd = 1.*(w.^2 + 1).^(1/2);
        q = (tau.^2 - tauwd.^2).^(1/2).*cos(theta) + 1i.*tau.*sin(theta);
        md = (q.^2+w.^2+1).^(1/2);
        E = 1 + q.^2.*D + w.^2.*F;
        O = -2.*q.*w.*(S.^2 - N.^2).*cos(phi).*sin(phi);
        dqddt = (cos(theta).*tau.*(tau.^2-tauwd.^2).^(-1/2) + 1i.*sin(theta));
        
        Md = 1i*(q.^2+w.^2+.5*taus^2).*q*cos(phi);
        Nd = -1i*(q.^2+w.^2+.5*taus^2).*w*sin(phi);
        
        res = real( ((E.^2 + O.^2).*Md - 2*E.*O.*Nd)./((E.^2 - O.^2).^2 ) .* dqddt);
    end

    function res = Sd_zz(w,tau,wa,qa)
        alpha = 'd';
        
        md = (qa^2 + wa^2 + 1)^(1/2);
        aa3 = a3(wa,alpha,tau);
        
        Va  = 1i*a4(wa,qa)*a13(wa,qa);
        Wap = 1i*(a8(aa3)*a13(wa,qa) + 2*a4(wa,qa)*a10(wa,qa,aa3));
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa./(w - wa).^2 + Ya./(w - wa));
    
    end

    function res = int_Sd_zz(Ta,tau,wa,qa)
        alpha = 'd';
        
        md = (qa^2 + wa^2 + 1)^(1/2);
        aa3 = a3(wa,alpha,tau);
        
        Va  = 1i*a4(wa,qa)*a13(wa,qa);
        Wap = 1i*(a8(aa3)*a13(wa,qa) + 2*a4(wa,qa)*a10(wa,qa,aa3));
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa*Ta./(wa*(wa - Ta)) + Ya*(log(abs(1-Ta/wa)) + 1i*angle(1-Ta/wa)));
    
    end

    function res = integrand_Is_zz(w,tau)
        
        tauws = 1.*(w.^2 + taus.^2).^(1/2);
        q = (tau.^2 - tauws.^2).^(1/2).*cos(theta) + 1i.*tau.*sin(theta);
        ms = (q.^2 + w.^2+ taus^2).^(1/2);
        E = 1 + q.^2.*D + w.^2.*F;
        O = -2.*q.*w.*(S.^2 - N.^2).*cos(phi).*sin(phi);
        dqsdt = (cos(theta).*tau.*(tau.^2-tauws.^2).^(-1/2) + 1i.*sin(theta));
        
        Ms = -1i*(q.^2+w.^2+.5*taus^2).*q*cos(phi);
        Ns = 1i*(q.^2+w.^2+.5*taus^2).*w*sin(phi);
        
        res = real( ((E.^2 + O.^2).*Ms - 2*E.*O.*Ns)./((E.^2 - O.^2).^2 ) .* dqsdt);
    end

    function res = Ss_zz(w,tau,wa,qa)
        alpha = 's';
        
        aa3 = a3(wa,alpha,tau);
        
        Va  = -1i*a4(wa,qa)*a13(wa,qa);
        Wap = -1i*(a8(aa3)*a13(wa,qa) + 2*a4(wa,qa)*a10(wa,qa,aa3));
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa./(w - wa).^2 + Ya./(w - wa));
    
    end

    function res = int_Ss_zz(Ta,tau,wa,qa)
        alpha = 's';
        
        aa3 = a3(wa,alpha,tau);
        
        Va  = -1i*a4(wa,qa)*a13(wa,qa);
        Wap = -1i*(a8(aa3)*a13(wa,qa) + 2*a4(wa,qa)*a10(wa,qa,aa3));
        
        Xa  = a18(wa,alpha,tau)*Va./(8*a14(wa,qa,alpha,tau)^2);
        Ya  = ( (a19(wa,alpha,tau)-a15(wa,qa,alpha,tau)*a18(wa,alpha,tau))*Va + a18(wa,alpha,tau)*Wap )./(8*a14(wa,qa,alpha,tau)^2);
    
        res = real( Xa*Ta./(wa*(wa - Ta)) + Ya*(log(abs(1-Ta/wa)) + 1i*angle(1-Ta/wa)));
    
    end

%% residual terms for ux

    function res = Rd_ux_calc(w,q)
        alpha = 'd';
        
        aa1  = a1(w);
        aa2  = a2(w,q,alpha);
        aa11 = a11(w,q,alpha,0);
        aa12 = a12(w,q,alpha,0);
        
        Ua = -2*a4(w,q)*cos(phi);
        Va = -a4(w,q)^2;
        Wa = -2*a4(w,q)*a8(a3(w,alpha,0));
        
        res = 2*pi*real( 1i*(-( Ua + 1i*Va*D/aa1 )*aa2/(8*aa1^2) + aa11*aa2*Wa + aa12*Va ) );
        
    end

    function res = Rs_ux_calc(w,q)
        alpha = 's';
        
        aa1  = a1(w);
        aa2  = a2(w,q,alpha);
        aa11 = a11(w,q,alpha,0);
        aa12 = a12(w,q,alpha,0);
        
        Ua = 2*a4(w,q)*cos(phi);
        Va = a4(w,q)^2 + .5*taus^2;
        Wa = 2*a4(w,q)*a8(a3(w,alpha,0));
        
        res = 2*pi*real( 1i*(-( Ua + 1i*Va*D/aa1 )*aa2/(8*aa1^2) + aa11*aa2*Wa + aa12*Va ) );
        
    end

%% residual terms for uy

    function res = Rd_uy_calc(w,q)
        alpha = 'd';
        
        aa1  = a1(w);
        aa2  = a2(w,q,alpha);
        aa11 = a11(w,q,alpha,0);
        aa12 = a12(w,q,alpha,0);
        
        Ua = -a6(w,q);
        Va = -a5(w,q);
        Wa = -a9(w,q,a3(w,q,0));
        
        res = 2*pi*real( 1i*(-( Ua + 1i*Va*D/aa1 )*aa2/(8*aa1^2) + aa11*aa2*Wa + aa12*Va ) );
        
    end

    function res = Rs_uy_calc(w,q)
        alpha = 's';
        
        aa1  = a1(w);
        aa2  = a2(w,q,alpha);
        aa11 = a11(w,q,alpha,0);
        aa12 = a12(w,q,alpha,0);
        
        Ua = a6(w,q);
        Va = a5(w,q);
        Wa = a9(w,q,a3(w,q,0));
        
        res = 2*pi*real( 1i*(-( Ua + 1i*Va*D/aa1 )*aa2/(8*aa1^2) + aa11*aa2*Wa + aa12*Va ) );
        
    end

%% residual terms for uz

    function res = Rd_uz_calc(w,q)
        md = (q^2+w^2+1).^(1/2);
        alpha = 'd';
        
        aa1  = a1(w);
        aa2  = a2(w,q,alpha);
        aa11 = a11(w,q,alpha,0);
        aa12 = a12(w,q,alpha,0);
        
        Ua = -1i*(md*cos(phi) + q*a4(w,q)/md);
        Va = -1i*a4(w,q)*md;
        Wa = -1i*(a8(a3(w,alpha,0))*md + a4(w,q)*a10(w,q,a3(w,alpha,0))/md);
        
        res = 2*pi*real( 1i*(-( Ua + 1i*Va*D/aa1 )*aa2/(8*aa1^2) + aa11*aa2*Wa + aa12*Va ) );
        
    end

    function res = Rs_uz_calc(w,q)
        alpha = 's';
        
        aa1  = a1(w);
        aa2  = a2(w,q,alpha);
        aa11 = a11(w,q,alpha,0);
        aa12 = a12(w,q,alpha,0);
        
        Ua = 1i*(a71(w,q)*cos(phi) + a4(w,q)*a73(w,q));
        Va = 1i*a4(w,q)*a71(w,q);
        Wa = 1i*(a8(a3(w,alpha,0))*a71(w,q) + a4(w,q)*a10(w,q,a3(w,alpha,0))*a73(w,q)/q);
        
        res = 2*pi*real( 1i*(-( Ua + 1i*Va*D/aa1 )*aa2/(8*aa1^2) + aa11*aa2*Wa + aa12*Va ) );
        
    end

%% residual terms for tau_zx

    function res = Rd_zx_calc(w,q)
        alpha = 'd';
        
        aa1  = a1(w);
        aa2  = a2(w,q,alpha);
        aa11 = a11(w,q,alpha,0);
        aa12 = a12(w,q,alpha,0);
        
        md = (q^2+w^2+1)^(1/2);
        
        Ua = 2*a4(w,q)*md*cos(phi) + q*a4(w,q)^2/md;
        Va = a4(w,q)^2*md;
        Wa = 2*a4(w,q)*a8(a3(w,q,0))*md + a4(w,q)^2*a10(w,q,a3(w,q,0))/md;
        
        res = 2*pi*real( 1i*(-( Ua + 1i*Va*D/aa1 )*aa2/(8*aa1^2) + aa11*aa2*Wa + aa12*Va ) );
        
    end

    function res = Rs_zx_calc(w,q)
        alpha = 's';
        
        aa1  = a1(w);
        aa2  = a2(w,q,alpha);
        aa11 = a11(w,q,alpha,0);
        aa12 = a12(w,q,alpha,0);
        
        Ua = -2*a4(w,q)*a72(w,q)*cos(phi) - a4(w,q)^2*a74(w,q) - .25*taus^2*q/a7(w,q);
        Va = -a4(w,q)^2*a72(w,q) - .25*taus^2*a7(w,q);
        Wa = -2*a4(w,q)*a8(a3(w,q,0))*a72(w,q) - (a4(w,q)^2*a74(w,q)/q + .25*taus^2/a7(w,q))*a10(w,q,a3(w,q,0));
        
        res = 2*pi*real( 1i*(-( Ua + 1i*Va*D/aa1 )*aa2/(8*aa1^2) + aa11*aa2*Wa + aa12*Va ) );
        
    end

%% residual terms for tau_zy

    function res = Rd_zy_calc(w,q)
        alpha = 'd';
        
        aa1  = a1(w);
        aa2  = a2(w,q,alpha);
        aa11 = a11(w,q,alpha,0);
        aa12 = a12(w,q,alpha,0);
        
        md = (q^2+w^2+1)^(1/2);
        
        Ua = a6(w,q)*md + q*a5(w,q)/md;
        Va = a5(w,q)*md;
        Wa = a9(w,q,a3(w,q,0))*md + a5(w,q)*a10(w,q,a3(w,q,0))/md;
        
        res = 2*pi*real( 1i*(-( Ua + 1i*Va*D/aa1 )*aa2/(8*aa1^2) + aa11*aa2*Wa + aa12*Va ) );
        
    end

    function res = Rs_zy_calc(w,q)
        alpha = 's';
        
        aa1  = a1(w);
        aa2  = a2(w,q,alpha);
        aa11 = a11(w,q,alpha,0);
        aa12 = a12(w,q,alpha,0);
        
        Ua = -a6(w,q)*a72(w,q) - a5(w,q)*a74(w,q);
        Va = -a5(w,q)*a72(w,q);
        Wa = -a9(w,q,a3(w,q,0))*a72(w,q) - a5(w,q)*a10(w,q,a3(w,q,0))*a74(w,q)/q;
        
        res = 2*pi*real( 1i*(-( Ua + 1i*Va*D/aa1 )*aa2/(8*aa1^2) + aa11*aa2*Wa + aa12*Va ) );
        
    end

%% residual terms for tau_zz

    function res = Rd_zz_calc(w,q)
        alpha = 'd';
        
        aa1  = a1(w);
        aa2  = a2(w,q,alpha);
        aa11 = a11(w,q,alpha,0);
        aa12 = a12(w,q,alpha,0);
        
        Ua = 1i*(a13(w,q)*cos(phi) + 2*q*a4(w,q));
        Va = 1i*a4(w,q)*a13(w,q);
        Wa = 1i*(a8(a3(w,q,0))*a13(w,q) + 2*a4(w,q)*a10(w,q,a3(w,q,0)));
        
        res = 2*pi*real( 1i*(-( Ua + 1i*Va*D/aa1 )*aa2/(8*aa1^2) + aa11*aa2*Wa + aa12*Va ) );
        
    end

    function res = Rs_zz_calc(w,q)
        alpha = 's';
        
        aa1  = a1(w);
        aa2  = a2(w,q,alpha);
        aa11 = a11(w,q,alpha,0);
        aa12 = a12(w,q,alpha,0);
        
        Ua = -1i*(a13(w,q)*cos(phi) + 2*q*a4(w,q));
        Va = -1i*a4(w,q)*a13(w,q);
        Wa = -1i*(a8(a3(w,q,0))*a13(w,q) + 2*a4(w,q)*a10(w,q,a3(w,q,0)));
        
        res = 2*pi*real( 1i*(-( Ua + 1i*Va*D/aa1 )*aa2/(8*aa1^2) + aa11*aa2*Wa + aa12*Va ) );
        
    end

%% now all the a's ...

    function res = a1(w)
        res = (w.^2.*S^2*N^2 + D).^(1/2);
    end

    function res = a2(w,q,alpha)
        tmp = ( Delta + 1i*w*S^2*N^2.*(w^2*S^2*N^2+D)^(-1/2) )/D;
        if strcmp(alpha,'d')
            md = (q^2 + w^2 + 1)^(1/2);
            res = 1./( -1i*sin(theta)*tmp + md^(-1)*cos(theta)*(q*tmp + w) );
        elseif strcmp(alpha,'s')
            ms = (q^2 + w^2 + taus^2)^(1/2);
            res = 1./( -1i*sin(theta)*tmp + ms^(-1)*cos(theta)*(q*tmp + w) );
        else
            disp('Wrong wave type in a2');
            return
        end
    end

    function res = a3(w,alpha,tau)
        if tau==0
            res = ( Delta + 1i*w*S^2*N^2.*(w^2*S^2*N^2+D)^(-1/2) )/D;
        elseif strcmp(alpha,'d')
            tauwd = 1.*(w.^2 + 1).^(1/2);
            res = -w*(tau^2 - tauwd.^2).^(-1/2).*cos(theta);
        elseif strcmp(alpha,'s')
            tauws = 1.*(w.^2 + taus.^2).^(1/2);
            res = -w*(tau^2 - tauws.^2).^(-1/2).*cos(theta);
        else
            disp('Wrong wave type in a3');
            return
        end
    end

    function res = a4(w,q)
        if (S^2 - N^2)*sin(2*phi) <0
            res = q*cos(phi) - w*sin(phi);
        else
            res = q*cos(phi) + w*sin(phi);
        end
    end

    function res = a5(w,q)
        if (S^2 - N^2)*sin(2*phi) <0
            res = (q^2 - w^2)*cos(phi)*sin(phi) + q*w*(cos(phi)^2 - sin(phi)^2);
        else
            res = (q^2 - w^2)*cos(phi)*sin(phi) - q*w*(cos(phi)^2 - sin(phi)^2);
        end
    end

    function res = a6(w,q)
        if (S^2 - N^2)*sin(2*phi) <0
            res = 2*q*cos(phi)*sin(phi) + w*(cos(phi)^2 - sin(phi)^2);
        else
            res = 2*q*cos(phi)*sin(phi) - w*(cos(phi)^2 - sin(phi)^2);            
        end
    end

    function res = a7(w,q)
        res = (q^2+w^2+taus^2)^(1/2);
    end

    function res = a71(w,q)
        ms = (q^2 + w^2 + taus^2)^(1/2);
        res = ms - .5*taus^2/ms;
    end

    function res = a72(w,q)
        ms = (q^2 + w^2 + taus^2)^(1/2);
        res = ms - .25*taus^2/ms;
    end

    function res = a73(w,q)
        ms = (q^2 + w^2 + taus^2)^(1/2);
        res = q/ms + .5*q*taus^2/ms^3;
    end

    function res = a74(w,q)
        ms = (q^2 + w^2 + taus^2)^(1/2);
        res = q/ms + .25*q*taus^2/ms^3;
    end

    function res = a8(aa3)
        if (S^2 - N^2)*sin(2*phi) <0
            res = aa3*cos(phi)- sin(phi);
        else
            res = aa3*cos(phi)+ sin(phi);            
        end
    end

    function res = a9(w,q,aa3)
        if (S^2 - N^2)*sin(2*phi) <0
            res = 2*(q*aa3 - w)*cos(phi)*sin(phi) + ...
                (w*aa3 + q)*(cos(phi)^2-sin(phi)^2);
        else
            res = 2*(q*aa3 - w)*cos(phi)*sin(phi) - ...
                (w*aa3 + q)*(cos(phi)^2-sin(phi)^2);         
        end
    end

    function res = a10(w,q,aa3)
        res = aa3*q + w;
    end

    function res = a11(w,q,alpha,tau)
        if strcmp(alpha,'d')
            md = (q^2 + w^2 + 1)^(1/2);
            res = (8*a1(w)^2*(a3(w,alpha,tau) + w*cos(theta)/(q*cos(theta) - 1i*md*sin(theta))))^(-1);
        elseif strcmp(alpha,'s')
            ms = (q^2 + w^2 + taus^2)^(1/2);
            res = (8*a1(w)^2*(a3(w,alpha,tau) + w*cos(theta)/(q*cos(theta) - 1i*ms*sin(theta))))^(-1);
        else
            disp('Wrong wave type in a11');
            return
        end
    end

    function res = a12(w,q,alpha,tau)
        aa3  = a3(w,alpha,tau);
        aa11 = a11(w,q,alpha,tau);
        aa1  = a1(w);
        aa2  = a2(w,q,alpha);
        
        da1dt = w*S^2*N^2*aa2/aa1;
        da3dt = 1i*S^2*N^2*aa2/aa1 *(1 - w^2*S^2*N^2/aa1^2)/D;
        
        if strcmp(alpha,'d')
            ma  = (q^2+w^2+1)^(1/2);
        elseif strcmp(alpha,'s')
            ma  = (q^2+w^2+taus^2)^(1/2);
        else
            disp('Wrong wave type in a11');
            return
        end
        
        dmadt = aa2*(q*aa3 + w)/ma;

        res = -aa11^2*( 16*aa1*da1dt*(aa3 + w*cos(theta)./(q*cos(theta) - 1i*ma*sin(theta))) + ...
            8*aa1^2*( da3dt + aa2*cos(theta)./(q*cos(theta) - 1i*ma*sin(theta)) - w*cos(theta)./(q*cos(theta) - 1i*ma*sin(theta))^2 *(aa2*aa3*cos(theta) - 1i*sin(theta)*dmadt)));
    end

    function res = a13(w,q)
        res = q^2 + w^2 + .5*taus^2;
    end

    function res = a14(w,q,alpha,tau)
        if strcmp(alpha,'d')
            tauwa = 1.*(w.^2 + 1).^(1/2);
        elseif strcmp(alpha,'s')
            tauwa = 1.*(w.^2 + taus^2).^(1/2);
        else
            disp('Wrong wave type in a14');
        end
        
        dqadw = -w*(tau^2 - tauwa^2)^(-1/2) *cos(theta);
        
        res = w*F - q*Delta +1i*a1(w)*dqadw;
    end

    function res = a15(w,q,alpha,tau)
        if strcmp(alpha,'d')
            tauwa = 1.*(w.^2 + 1).^(1/2);
        elseif strcmp(alpha,'s')
            tauwa = 1.*(w.^2 + taus^2).^(1/2);
        else
            disp('Wrong wave type in a15');
        end
        
        dqadw = -w*(tau^2 - tauwa^2)^(-1/2) *cos(theta);
        d2qadw2 = -cos(theta)*(tau^2 - tauwa^2)^(-1/2)*(1 + w^2./(tau^2-tauwa^2));
        
        aa17 = a17(w,q);
        aa14 = a14(w,q,alpha,tau);
        aa1  = a1(w);
        
        res = F*aa17 - dqadw*(Delta + q*S^2*N^2*aa17)*aa17 + (1i*aa1*d2qadw2 - dqadw^2*S^2*N^2*aa17^2)./aa14;
        
    end

    function res = a16(w,q,alpha,tau)
        if strcmp(alpha,'d')
            tauwa = 1.*(w.^2 + 1).^(1/2);
        elseif strcmp(alpha,'s')
            tauwa = 1.*(w.^2 + taus^2).^(1/2);
        else
            disp('Wrong wave type in a16');
        end
        
        dqadw = -w*(tau^2 - tauwa^2)^(-1/2) *cos(theta);
        
        res = 1./w + F./(q*Delta) + dqadw*(1./q + D./(w*D));
    end

    function res = a17(w,q)
        res = 1./(w*F - q*Delta);
    end

    function res = a18(w,alpha,tau)
        if strcmp(alpha,'d')
            taua = 1;
        elseif strcmp(alpha,'s')
            taua = taus;
        else
            disp('Wrong wave type in a18');
        end
        
        res = tau*(tau^2 - w^2 - taua^2)^(-1/2)*cos(theta) + 1i*sin(theta);
        
    end

    function res = a19(w,alpha,tau)
        if strcmp(alpha,'d')
            taua = 1;
        elseif strcmp(alpha,'s')
            taua = taus;
        else
            disp('Wrong wave type in a18');
        end
        
        res = w*tau*(tau^2 - w^2 - taua^2)^(-3/2)*cos(theta);
        
    end


end

