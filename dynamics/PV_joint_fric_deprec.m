function  [qdd, nu, tau_apply, a, Xa, v, ee_acc_f_ext_fun, q_acc_X_fun, tau_X_fun] = PV( model, q, qd, tau, f_ext, A_con, m_con )

import casadi.*;

% Popov-Vereschagin solver with joint friction assuming that the joint
% friction is independent of the constraint forces at the joints

%A_con and m_con define the motion constraints on the end-effector.

a_grav = get_gravity(model);

qdd = SX.zeros(model.NB,1);
%qdd = zeros(model.NB,1);

for i = 1:model.NB
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
    vJ = S{i}*qd(i);
    Xup{i} = XJ * model.Xtree{i};
    if model.parent(i) == 0
        v{i} = vJ;
        c{i} = zeros(size(a_grav));		% spatial or planar zero vector
    else
        v{i} = Xup{i}*v{model.parent(i)} + vJ;
        c{i} = crm(v{i}) * vJ;
    end
    IA{i} = model.I{i};
    pA{i} = crf(v{i}) * (model.I{i} * v{i});
    
end

parent = model.parent;
    for i = 1:length(parent)
        if parent(i) == 0
            Xa{i} = Xup{i};
        else
            Xa{i} = Xup{i} * Xa{parent(i)};
        end
    end

if nargin >= 5 && size(f_ext,1) >0
    pA = apply_external_forces( model.parent, Xup, pA, f_ext );
end

%PV
if size(A_con, 2) > 0
    
    
    a_grav_ee = Xa{7} * (-a_grav);
    F_cstr{model.NB} = A_con;
    E_cstr{model.NB} = zeros (size(A_con,2), size(A_con,2)); %SX(size(A_con,2), size(A_con,2));
    E_conf{model.NB} = F_cstr{model.NB}' * (-a_grav_ee);
end

for i = model.NB:-1:1
    U{i} = IA{i} * S{i};
    d{i} = S{i}' * U{i};
    u{i} = tau(i) - S{i}'*pA{i};
    eqd = exp(qd(i)*10);
    meqd = exp(-qd(i)*10);
    u{i} = u{i} - sign(qd(i))*model.joint_frictions(i);
    
    Ia = IA{i} - U{i}/d{i}*U{i}';
    pa = pA{i} + Ia*c{i} + U{i} * u{i}/d{i};
    if model.parent(i) ~= 0
        IA{model.parent(i)} = IA{model.parent(i)} + Xup{i}' * Ia * Xup{i};
        pA{model.parent(i)} = pA{model.parent(i)} + Xup{i}' * pa;
        
        
        if size(A_con, 2) > 0
            %PV
            FS = F_cstr{i}'*S{i};
            F_cstr{model.parent(i)} = Xup{i}' * (F_cstr{i} - U{i}/d{i}*(FS'));
            E_cstr{model.parent(i)} = E_cstr{i} - FS/d{i}*(FS');
            E_conf{model.parent(i)} = E_conf{i} + F_cstr{i}' * (c{i} + S{i}/d{i}*(u{i} - S{i}'*(IA{i}*c{i})));
        end
    else
        
        if size(A_con, 2) > 0
            %PV
            FS = F_cstr{i}'*S{i};
            F_cstr0 = Xup{i}' * (F_cstr{i} - U{i}/d{i}*(FS'));
            E_cstr0 = E_cstr{i} - FS/d{i}*(FS');
            E_conf0 = E_conf{i} + F_cstr{i}' * (c{i} + S{i}/d{i}*(tau(i) - S{i}'*(pA{i} + IA{i}*c{i})));
        end
    end
    
end

if size(A_con, 2) > 0

    E_cstr0 = E_cstr0 + eye(size(A_con,2), size(A_con,2))*1e-6;
    b = F_cstr0'*(-a_grav) + E_conf0 - A_con'*Xa{7}*[m_con];%Xa{7}(1:3, 1:3)*[m_con];%m_con; %m_con needs to be in body frame. Otherwise use Xa{7}*[m_con];
    Lchol = cholesky(-E_cstr0);
    %inv_dual_hess = inv(E_cstr0);
    nu = back_sub(Lchol', forward_sub(Lchol,b)); 
    %nu = -inv_dual_hess * b;
    % nu = -solve(E_cstr0, b);
    %nu = fmax(-500, fmin(nu, 500));
else
    nu = 0;
end

con_torque = SX(model.NB, 1);

for i = 1:model.NB
    if model.parent(i) == 0
        a{i} = Xup{i} * -a_grav + c{i};
    else
        a{i} = Xup{i} * a{model.parent(i)} + c{i};
    end
    if size(A_con, 2) > 0
        con_torque(i) = -S{i}' * F_cstr{i} * nu;
    else
        con_torque(i) = 0;
    end
    qdd(i,1) = (u{i} - U{i}'*a{i} + con_torque(i))/d{i};
    a{i} = a{i} + S{i}*qdd(i);
end

%In world frame final segment has acceleration

%final_seg = Xa{7} *a{7};

%tau_apply = SX.zeros(model.NB, 1);
% for i = 1:model.NB
%     tau_apply(i) = con_torque(i) + tau(i);
% end



tau_apply = SX.zeros(model.NB, 1);
p_con_T{model.NB} = zeros(6,1);
tau_apply(model.NB) = con_torque{model.NB} + tau(model.NB);
for i = model.NB-1:-1:1
    p_con_T{i} = Xup{i+1}' * (p_con_T{i+1} + U{i+1} * (-S{i+1}'*p_con_T{i+1} + con_torque{i+1}) /d{i+1});
    con_torque(i) = con_torque(i) + S{i}' * p_con_T{i};
    tau_apply(i) = con_torque(i) + tau(i);
end


ee_acc_f_ext_fun = 0;
ee_acc_X_fun = 0;
q_acc_X_fun = 0;
tau_X_fun = 0;

if nargin >= 5 && size(f_ext,1) >0
    ee_acc_f_ext_fun = jacobian(a{7}, f_ext{7});

end

if size(A_con, 2) > 0
    q_acc_X_fun = jacobian([qdd;con_torque; nu], [m_con; tau]);
    tau_X_fun = q_acc_X_fun(8:end,:);
    q_acc_X_fun = q_acc_X_fun(1:7,:);
end



