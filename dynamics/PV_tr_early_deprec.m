%PV implementation for kinematic trees

function  [qdd, nu, tau_apply, a, Xa, v, q_acc_X_fun, tau_X_fun] = PV_tr_early( model, q, qd, tau, f_ext, A_con, m_con, EEs )

import casadi.*;

% FDab  Forward Dynamics via Articulated-Body Algorithm
% FDab(model,q,qd,tau,f_ext,grav_accn)  calculates the forward dynamics of
% a kinematic tree via the articulated-body algorithm.  q, qd and tau are
% vectors of joint position, velocity and force variables; and the return
% value is a vector of joint acceleration variables.  f_ext is an optional
% argument specifying the external forces acting on the bodies.  It can be
% omitted if there are no external forces.  The format of f_ext is
% explained in the source code of apply_external_forces.

%A_con and m_con define the motion constraints on the end-effector.

n = size(q, 1);
a_grav = get_gravity(model);

qdd = SX.zeros(model.NB,1);
%qdd = zeros(model.NB,1);
F_cstr0 = [];
nu = SX.zeros(6*length(A_con),1);
nu_assigned = zeros(6*length(A_con),1);
nu_val_indices = {};
counter_rank = {};
for i = 1:model.NB
    counter_rank{i} = 0;
end
nu_val_counter = 0;

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
if length(EEs) > 0
    
    for i = 1:length(EEs)
        a_grav_ee = Xa{EEs{i}} * (-a_grav);
        F_cstr{EEs{i}} = A_con{i};
        nu_val_indices{EEs{i}} = nu_val_counter + 1 : nu_val_counter + size(F_cstr{EEs{i}}, 2);
        nu_val_counter = nu_val_counter + size(F_cstr{EEs{i}}, 2);
        E_cstr{EEs{i}} = zeros (size(A_con{i},2), size(A_con{i},2)); %SX(size(A_con,2), size(A_con,2));
        E_conf{EEs{i}} = F_cstr{EEs{i}}' * (-a_grav_ee) - A_con{i}'*Xa{EEs{i}}*[m_con{i}];
        counter_rank{EEs{i}} = 6;
    end
end

for i = model.NB:-1:1
    U{i} = IA{i} * S{i};
    d{i} = S{i}' * U{i};
    u{i} = tau(i) - S{i}'*pA{i};
    
    Ia = IA{i} - U{i}/d{i}*U{i}';
    Ia = casadi_symmetric(Ia);
    pa = pA{i} + Ia*c{i} + U{i} * u{i}/d{i};
    if model.parent(i) ~= 0
        IA{model.parent(i)} = IA{model.parent(i)} + Xup{i}' * Ia * Xup{i};
        IA{model.parent(i)} = casadi_symmetric(IA{model.parent(i)});
%         for sim_i = 2:6
%             for sim_j = 1:sim_i-1
%                 IA{model.parent(i)}(sim_i, sim_j) = IA{model.parent(i)}(sim_j, sim_i);
%             end
%         end
        pA{model.parent(i)} = pA{model.parent(i)} + Xup{i}' * pa;
        
        
        if size( F_cstr{i}, 2) > 0
            %PV
            FS = F_cstr{i}'*S{i};
            if size(F_cstr{model.parent(i)}, 1) == 0 && counter_rank{i} > 0
                nu_val_indices{model.parent(i)} = nu_val_indices{i};
                F_cstr{model.parent(i)} = Xup{i}' * (F_cstr{i} - U{i}/d{i}*(FS'));
                E_cstr{model.parent(i)} = E_cstr{i} - FS/d{i}*(FS');
                E_cstr{model.parent(i)} = casadi_symmetric(E_cstr{model.parent(i)});
                E_conf{model.parent(i)} = E_conf{i} + F_cstr{i}' * (c{i} + S{i}/d{i}*(tau(i) - S{i}'*(pA{i} + IA{i}*c{i})));
                counter_rank{model.parent(i)} = counter_rank{i} - 1;
                if counter_rank{i} == 1
                    inv_E{model.parent(i)} = -inv(E_cstr{model.parent(i)} + SX.eye(size(F_cstr{model.parent(i)}, 2))*1e-6);
                    inv_E{model.parent(i)} = casadi_symmetric(inv_E{model.parent(i)});
                    IA{model.parent(i)} = IA{model.parent(i)} + F_cstr{model.parent(i)}*inv_E{model.parent(i)}*F_cstr{model.parent(i)}';
                    IA{model.parent(i)} = casadi_symmetric(IA{model.parent(i)});
                    E_conf{model.parent(i)}  = E_conf{model.parent(i)};
                    pA{model.parent(i)} = pA{model.parent(i)} + F_cstr{model.parent(i)}*(inv_E{model.parent(i)}*(E_conf{model.parent(i)}));
                end
%             else
%                 nu_val_indices{model.parent(i)} = [nu_val_indices{model.parent(i)}, nu_val_indices{i}];
%                 F_cstr{model.parent(i)} = [F_cstr{model.parent(i)}, Xup{i}' * (F_cstr{i} - U{i}/d{i}*(FS'))];
%                 zero_off_diag_LM = zeros(size(E_cstr{model.parent(i)}, 1), size(E_cstr{i}, 2));
%                 E_cstr{model.parent(i)} = [E_cstr{model.parent(i)}, zero_off_diag_LM; zero_off_diag_LM', E_cstr{i} - FS/d{i}*(FS')];
%                 E_cstr{model.parent(i)} = casadi_symmetric(E_cstr{model.parent(i)});
%                 E_conf{model.parent(i)} = [E_conf{model.parent(i)}; E_conf{i} + F_cstr{i}' * (c{i} + S{i}/d{i}*(tau(i) - S{i}'*(pA{i} + IA{i}*c{i})))]; 
            end
        end
    else
        
        if size(F_cstr{i}, 2) > 0
            %PV
            
            if size(F_cstr0, 2) == 0
                FS = F_cstr{i}'*S{i};
                F_cstr0 = Xup{i}' * (F_cstr{i} - U{i}/d{i}*(FS'));
                E_cstr0 = E_cstr{i} - FS/d{i}*(FS');
                E_cstr0 = casadi_symmetric(E_cstr0);
                E_conf0 = E_conf{i} + F_cstr{i}' * (c{i} + S{i}/d{i}*(tau(i) - S{i}'*(pA{i} + IA{i}*c{i})));
                
            else
                FS = F_cstr{i}'*S{i};
                F_cstr0 = [F_cstr0, Xup{i}' * (F_cstr{i} - U{i}/d{i}*(FS'))];
                zero_off_diag_LM = cs.MX.zeros(size(E_cstr0,1), size(E_cstr{i},2));
                E_cstr0 = [E_cstr0, zero_off_diag_LM; zero_off_diag_LM', E_cstr{i} - FS/d{i}*(FS')];
                E_cstr0 = casadi_symmetric(E_cstr0);
                E_conf0 = [E_conf0; E_conf{i} + F_cstr{i}' * (c{i} + S{i}/d{i}*(tau(i) - S{i}'*(pA{i} + IA{i}*c{i})))];
            end        
        end
    end
    
end

% if size(A_con, 2) > 0
% 
%     E_cstr0 = E_cstr0;% + eye(size(A_con,2), size(A_con,2))*1e-6;
%     b = F_cstr0'*(-a_grav) + E_conf0; % - A_con'*Xa{n}*[m_con];%Xa{n}(1:3, 1:3)*[m_con];%m_con; %m_con needs to be in body frame. Otherwise use Xa{n}*[m_con];
%     Lchol = cholesky(-E_cstr0);
%     %inv_dual_hess = inv(E_cstr0);
%     nu = back_sub(Lchol', forward_sub(Lchol,b)); 
%     %nu = -inv_dual_hess * b;
%     % nu = -solve(E_cstr0, b);
%     %nu = fmax(-500, fmin(nu, 500));
% else
%     nu = 0;
% end

con_torque = SX(model.NB, 1);

for i = 1:model.NB
    if model.parent(i) == 0
        a{i} = Xup{i} * -a_grav + c{i};
    else
        a{i} = Xup{i} * a{model.parent(i)} + c{i};
    end
    if size(F_cstr{i}, 2) > 0 && sumsqr(nu_assigned(nu_val_indices{i})) > 0
        con_torque(i) = -S{i}' * F_cstr{i} * nu(nu_val_indices{i});
    else
        con_torque(i) = 0;
    end
    qdd(i,1) = (u{i} - U{i}'*a{i} + con_torque(i))/d{i};
    a{i} = a{i} + S{i}*qdd(i);
    if size(F_cstr{i},2) > 0 && sumsqr(nu_assigned(nu_val_indices{i})) == 0
            nu(nu_val_indices{i}) = inv_E{i}*(E_conf{i} + F_cstr{i}'*a{i});
            nu_assigned(nu_val_indices{i}) = 1;
    end
end

%In world frame final segment has acceleration

%final_seg = Xa{n} *a{n};

%tau_apply = SX.zeros(model.NB, 1);
% for i = 1:model.NB
%     tau_apply(i) = con_torque(i) + tau(i);
% end



tau_apply = SX.zeros(model.NB, 1);
% p_con_T{model.NB} = zeros(6,1);
% tau_apply(model.NB) = con_torque{model.NB} + tau(model.NB);
% for i = model.NB-1:-1:1
%     p_con_T{i} = Xup{i+1}' * (p_con_T{i+1} + U{i+1} * (-S{i+1}'*p_con_T{i+1} + con_torque{i+1}) /d{i+1});
%     con_torque(i) = con_torque(i) + S{i}' * p_con_T{i};
%     tau_apply(i) = con_torque(i) + tau(i);
% end


ee_acc_f_ext_fun = 0;
ee_acc_X_fun = 0;
q_acc_X_fun = 0;
tau_X_fun = 0;

% if nargin >= 5 && size(f_ext,1) >0
%     ee_acc_f_ext_fun = jacobian(a{n}, f_ext{n});
% 
% end

if size(A_con, 2) > 0
    q_acc_X_fun = jacobian([qdd;con_torque; nu], [tau; q]);
    tau_X_fun = q_acc_X_fun(n+1:end,:);
    q_acc_X_fun = q_acc_X_fun(1:n,:);
end



