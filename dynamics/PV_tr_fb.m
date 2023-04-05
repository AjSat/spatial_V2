%PV implementation for kinematic trees

function  [qdd, nu, xd_fb, Xs, Vs, As, LA, IA, KA, a_ee] = PV_tr_fb( model, x_fb, q, qd, tau, f_ext, K_con, k_con, Soft)

import casadi.*;

if strcmp(class(q), 'casadi.MX')
    cs = MX;
    csX = @MX;
else
    cs = SX;
    csX = @SX;
end

% FDab  Forward Dynamics via Articulated-Body Algorithm
% FDab(model,q,qd,tau,f_ext,grav_accn)  calculates the forward dynamics of
% a kinematic tree via the articulated-body algorithm.  q, qd and tau are
% vectors of joint position, velocity and force variables; and the return
% value is a vector of joint acceleration variables.  f_ext is an optional
% argument specifying the external forces acting on the bodies.  It can be
% omitted if there are no external forces.  The format of f_ext is
% explained in the source code of apply_external_forces.

%A_con and m_con define the motion constraints on the end-effector.

q = [0; q];
qd = [0; qd];
tau = [0; tau];

n = size(q, 1);
a_grav = get_gravity(model);

prop_a_grav = true;

qdd = cs.zeros(model.NB,1);

F_cstr0 = [];

nu_val_indices = {};
nu_val_counter = 0;

qn = x_fb(1:4);				% unit quaternion fixed-->f.b.
r = x_fb(5:7);				% position of f.b. origin
Xup{1} = plux( rq(qn), r );		% xform fixed --> f.b. coords

vfb = x_fb(8:end);
v{1} = Xup{1} * vfb;			% f.b. vel in f.b. coords

IA{1} = model.I{1};
pA{1} = crf(v{1}) * model.I{1} * v{1};
KA{1} = [];
LA{1} = [];
lA{1} = [];

if prop_a_grav
    a_grav_links{1} = Xup{1}*a_grav;
end


for i = 2:model.NB
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
    vJ = S{i}*qd(i);
    Xup{i} = XJ * model.Xtree{i};
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    c{i} = crm(v{i}) * vJ;
    IA{i} = model.I{i};
    pA{i} = crf(v{i}) * (model.I{i} * v{i});
    
    if prop_a_grav
            a_grav_links{i} = Xup{i}*a_grav_links{model.parent(i)};
    end
    
    if nargin > 8 && ~isempty(Soft{i})
        IA{i} = IA{i} + Soft{i}.Ki'*Soft{i}.Ri*Soft{i}.Ki;
        if strcmp(class(cs), 'casadi.SX')
            IA{i} = casadi_symmetric(IA{i});
        end
        if prop_a_grav
            Soft{i}.ki = Soft{i}.ki - Soft{i}.Ki*a_grav_links{i};
        end
        pA{i} = pA{i} - Soft{i}.Ki'*Soft{i}.Ri*Soft{i}.ki;
    end
    
end

nu_vals = [];
counter = 0;
for i = 1:model.NB
    KA{i} = K_con{i};
    m_i = size(K_con{i}, 1);
    LA{i} = csX(m_i, m_i);
    lA{i} = -k_con{i};
    
    if m_i > 0
        counter = counter + size(LA{1}, 1);
        nu_vals = [nu_vals, counter];
        
    end
    if m_i > 0 && prop_a_grav
        lA{i} = lA{i} + KA{i}*a_grav_links{i};
    end
    
end

parent = model.parent;
    for i = 1:length(parent)
        if parent(i) == 0
            Xa{i} = Xup{i};
        else
            Xa{i} = Xup{i} * Xa{parent(i)};
        end
    end

if nargin >= 6 && size(f_ext,1) >0
    pA = apply_external_forces( model.parent, Xup, pA, f_ext );
end

%PV
% if length(EEs) > 0
%     
%     for i = 1:length(EEs)
%         a_grav_ee = Xa{EEs{i}} * (-a_grav);
%         F_cstr{EEs{i}} = A_con{i};
%         nu_val_indices{EEs{i}} = nu_val_counter + 1 : nu_val_counter + size(F_cstr{EEs{i}}, 2);
%         nu_val_counter = nu_val_counter + size(F_cstr{EEs{i}}, 2);
%         E_cstr{EEs{i}} = zeros (size(A_con{i},2), size(A_con{i},2)); %csX(size(A_con,2), size(A_con,2));
%         E_conf{EEs{i}} = F_cstr{EEs{i}}' * (-a_grav_ee) - A_con{i}'*Xa{EEs{i}}*[m_con{i}];
%     end
% end

for i = model.NB:-1:2
    U{i} = IA{i} * S{i};
    d{i} = S{i}' * U{i};
    u{i} = tau(i) - S{i}'*pA{i};
    
    Ia = IA{i} - U{i}/d{i}*U{i}';
    if strcmp(class(cs), 'casadi.SX')
        Ia = casadi_symmetric(Ia);
    end
    pa = pA{i} + Ia*c{i} + U{i} * u{i}/d{i};

        IA{model.parent(i)} = IA{model.parent(i)} + Xup{i}' * Ia * Xup{i};
        if strcmp(class(cs), 'casadi.SX')
            IA{model.parent(i)} = casadi_symmetric(IA{model.parent(i)});
        end
        pA{model.parent(i)} = pA{model.parent(i)} + Xup{i}' * pa;
        
        
        if size( KA{i}, 1) > 0
            %PV          
            FS = KA{i}*S{i};
            KA{model.parent(i)} = [KA{model.parent(i)};  (KA{i} - FS/d{i}*U{i}')*Xup{i}];
            zero_off_diag_LM = csX(size(LA{model.parent(i)}, 1), size(KA{i}, 1));
            LA{model.parent(i)} = [LA{model.parent(i)}, zero_off_diag_LM; zero_off_diag_LM', LA{i} + FS/d{i}*(FS')];
            if strcmp(class(cs), 'casadi.SX')
                LA{model.parent(i)} = casadi_symmetric(LA{model.parent(i)});
            end
            lA{model.parent(i)} = [lA{model.parent(i)}; lA{i} + KA{i} * (c{i} + S{i}/d{i}*(tau(i) - S{i}'*(pA{i} + IA{i}*c{i})))]; 
        end
    
end

if isempty(lA{1})
    if strcmp(class(cs), 'casadi.MX')
        a{1} = -IA{1}\pA{1};
    else
        Lchol = cholesky(IA{1});
        a{1} = -back_sub(Lchol', forward_sub(Lchol,pA{1})); 
        nu = 0;
%         a_fun = Function('f_rob_dyn', {x_fb, q(2:end), qd(2:end), tau(2:end)}, {a{1}});
%         a_fun.n_instructions
    end

else
    if strcmp(class(cs), 'casadi.MX')
        OSIM = inv(LA{1} + cs.eye(size(LA{1},1))*1e-6, 'ldl');
    else
%         OSIM = inv(LA{1} + cs.eye(size(LA{1},1))*1e-6);
          Lchol_osim = cholesky(LA{1} + cs.eye(size(LA{1},1))*1e-6*0);
          
          % uncomment only for fast OSIM chol
%           Lambda_b = cholesky(LA{1});
%           Lambda_b = casadi_symmetric(Lambda_b);
%           Lcrossk = back_sub(Lambda_b', forward_sub(Lambda_b, KA{1}));
%           H_inv = inv(IA{1}); H_inv = casadi_symmetric(H_inv);
%           full_osim = Lambda_b + (Lcrossk * casadi_symmetric( H_inv + KA{1}'*Lcrossk))*Lcrossk';
%           full_osim = casadi_symmetric(full_osim);
%           invosim_fun = Function('f_osim', {q(2:end)}, {Lambda_b, casadi_symmetric(inv(casadi_symmetric(H_inv + KA{1}'*Lcrossk))), Lcrossk});
%           invosim_fun.n_instructions  
%           invosim_fun.generate(strcat(strcat('fast_osim_go_18D_', num2str(size(LA{1},1))), '.c'), struct('with_header', true))
%           
          

% uncomment only to compute OSIM chol
%           full_osim = LA{1} + KA{1}*inv(IA{1})*KA{1}';
%           full_osim = casadi_symmetric(full_osim);
%           %invosim_fun = Function('f_osim', {q(2:end)}, {casadi_symmetric(inv(full_osim))});
%           invosim_fun = Function('f_osim', {q(2:end)}, {cholesky(full_osim)});
%           invosim_fun.n_instructions
%           invosim_fun.generate(strcat(strcat('osiminv_atlas_18D_', num2str(size(LA{1},1))), '.c'), struct('with_header', true))
    end
%     FB_inertia = IA{1} + KA{1}'*OSIM*KA{1};
    FB_inertia = IA{1} + KA{1}'*back_sub(Lchol_osim', forward_sub(Lchol_osim, KA{1}));
%     b = -pA{1} - KA{1}'*OSIM*lA{1}; % - A_con'*Xa{n}*[m_con];%Xa{n}(1:3, 1:3)*[m_con];%m_con; %m_con needs to be in body frame. Otherwise use Xa{n}*[m_con];
    b = -pA{1} - KA{1}'*back_sub(Lchol_osim', forward_sub(Lchol_osim,lA{1}));
    if strcmp(class(cs), 'casadi.SX')
        FB_inertia = casadi_symmetric(FB_inertia);
        Lchol = cholesky(FB_inertia);
        a{1} = back_sub(Lchol', forward_sub(Lchol,b)); 
%         a_fun = Function('f_rob_dyn', {x_fb, q(2:end), qd(2:end), tau(2:end)}, {a{1}});
%         a_fun.n_instructions
    else
        a{1} = solve(FB_inertia, b, 'ldl');
    end
    
      nu = back_sub(Lchol_osim', forward_sub(Lchol_osim,(KA{1}*a{1} + lA{1})));
%     nu = OSIM*(KA{1}*a{1} + lA{1}); 

end

con_torque = csX(model.NB, 1);

lambda{1} = nu;
for i = 2:model.NB

    a{i} = Xup{i} * a{model.parent(i)} + c{i};
    if size(KA{i}, 2) > 0
        mi = size(KA{i}, 1);
        lambda{i} = lambda{model.parent(i)}(end - mi + 1 : end);
        lambda{model.parent(i)} = lambda{model.parent(i)}(1 : end - mi);
        con_torque(i) = -S{i}' * KA{i}' * lambda{i};
    else
        con_torque(i) = 0;
    end
    qdd(i,1) = (u{i} - U{i}'*a{i} + con_torque(i))/d{i};
    a{i} = a{i} + S{i}*qdd(i);
end

a_ee = Xa{model.NB}\a{model.NB};

a1_fb = Xup{1}\a{1} + a_grav;
xd_fb = [x_fb(8:end); a1_fb];

%In world frame final segment has acceleration
%final_seg = Xa{n} *a{n};

%tau_apply = SX.zeros(model.NB, 1);
% for i = 1:model.NB
%     tau_apply(i) = con_torque(i) + tau(i);
% end



tau_apply = cs.zeros(model.NB, 1);
% p_con_T{model.NB} = zeros(6,1);
% tau_apply(model.NB) = con_torque{model.NB} + tau(model.NB);
% for i = model.NB-1:-1:1
%     p_con_T{i} = Xup{i+1}' * (p_con_T{i+1} + U{i+1} * (-S{i+1}'*p_con_T{i+1} + con_torque{i+1}) /d{i+1});
%     con_torque(i) = con_torque(i) + S{i}' * p_con_T{i};
%     tau_apply(i) = con_torque(i) + tau(i);
% end

As = collect_terms(a, n);
Vs = collect_terms(v, n);
Xs = collect_terms(Xa, n);

ee_acc_f_ext_fun = 0;
ee_acc_X_fun = 0;
q_acc_X_fun = 0;
tau_X_fun = 0;

end

function x = collect_terms (a, n)
    
    x = [];
    for i = 1:n
        x = [x; a{i}];
    end
end

% if nargin >= 5 && size(f_ext,1) >0
%     ee_acc_f_ext_fun = jacobian(a{n}, f_ext{n});
% 
% end

% if size(KA{1}, 1) > 0
%     q_acc_X_fun = jacobian([qdd;con_torque; nu], [tau; q]);
%     tau_X_fun = q_acc_X_fun(n+1:end,:);
%     q_acc_X_fun = q_acc_X_fun(1:n,:);
% end



