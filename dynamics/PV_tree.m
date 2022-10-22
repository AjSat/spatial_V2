%PV implementation for kinematic trees

function  [qdd, nu, a_ee] = PV_tree( model, q, qd, tau, f_ext, K_con, k_con, Soft)

import casadi.*;

bfgs = false;

if strcmp(class(q), 'casadi.MX')
    cs = MX;
    csX = @MX;
else
    cs = SX;
    csX = @SX;
end

n = size(q, 1);
a_grav = get_gravity(model);

qdd = cs.zeros(model.NB,1);


for i = 1:model.NB
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
    vJ = S{i}*qd(i);
    Xup{i} = XJ * model.Xtree{i};
    new_rot = XJ(1:3, 1:3)*model.Xtree{i}(1:3, 1:3);
    off_diag = XJ(4:6, 1:3)*model.Xtree{i}(1:3, 1:3) + XJ(4:6, 4:6)*model.Xtree{i}(4:6, 1:3);
    Xup{i} = [new_rot, csX(3,3); off_diag, new_rot];
    
    if model.parent(i) == 0
        v{i} = vJ;
        c{i} = zeros(size(a_grav));
    else
        v{i} = Xup{i}*v{model.parent(i)} + vJ;
        c{i} = crm(v{i}) * vJ;
    end
    
    IA{i} = model.I{i};
    pA{i} = crf(v{i}) * (model.I{i} * v{i});
    
    if nargin >= 8 && ~isempty(Soft{i})
        IA{i} = IA{i} + Soft{i}.Ki'*Soft{i}.Ri*Soft{i}.Ki;
        if strcmp(class(cs), 'casadi.SX')
            IA{i} = casadi_symmetric(IA{i});
        end
        pA{i} = pA{i} + Soft{i}.Ki'*Soft{i}.Ri*Soft{i}.ki;
    end
    
end

for i = 1:model.NB
    KA{i} = K_con{i};
    m_i = size(K_con{i}, 1);
    if ~bfgs
        LA{i} = csX(m_i, m_i);
    else
        HA{i} = cs.eye(m_i)*1e6;
    end
    lA{i} = -k_con{i};
    
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

KA0 = [];
if ~bfgs
    LA0 = [];
else
    HA0 = [];
end
lA0 = [];

% if size(KA{model.NB}, 1) > 0
% % KA{model.NB} = KA{model.NB}*Xa{n}';
% a_grav_ee = Xa{n} * (-a_grav);
% lA{model.NB} = lA{model.NB} + KA{model.NB} * (-a_grav_ee);
% end

for i = model.NB:-1:1
    U{i} = IA{i} * S{i};
    d{i} = S{i}' * U{i};
    u{i} = tau(i) - S{i}'*pA{i};
    
    Ia = IA{i} - U{i}/d{i}*U{i}';
    if strcmp(class(cs), 'casadi.SX')
        Ia = casadi_symmetric(Ia);
    end
    pa = pA{i} + Ia*c{i} + U{i} * (u{i}/d{i});
    
    if model.parent(i) ~= 0
        
        IA{model.parent(i)} = IA{model.parent(i)} + Xup{i}' * Ia * Xup{i};
        if strcmp(class(cs), 'casadi.SX')
            IA{model.parent(i)} = casadi_symmetric(IA{model.parent(i)});
        end
        pA{model.parent(i)} = pA{model.parent(i)} + Xup{i}' * pa;
        
        
        if size( KA{i}, 1) > 0
            %PV
            FS = KA{i}*S{i};
            KA{model.parent(i)} = [KA{model.parent(i)};  (KA{i} - FS*(U{i}'/d{i}))*Xup{i}];
            
            if ~bfgs
                zero_off_diag_LM = csX(size(LA{model.parent(i)}, 1), size(KA{i}, 1));
                LA{model.parent(i)} = [LA{model.parent(i)}, zero_off_diag_LM; zero_off_diag_LM', LA{i} + FS/d{i}*(FS')];
                if strcmp(class(cs), 'casadi.SX')
                    LA{model.parent(i)} = casadi_symmetric(LA{model.parent(i)});
                end
            else
                zero_off_diag_LM = csX(size(HA{model.parent(i)}, 1), size(KA{i}, 1));
                inv_vec = HA{i}*FS;
                H_prop = HA{i} - (inv_vec/(d{i} + FS'*inv_vec))*inv_vec';
                HA{model.parent(i)} = [HA{model.parent(i)}, zero_off_diag_LM; zero_off_diag_LM', H_prop];
                if strcmp(class(cs), 'casadi.SX')
                    HA{model.parent(i)} = casadi_symmetric(HA{model.parent(i)});
                end
            end
            lA{model.parent(i)} = [lA{model.parent(i)}; lA{i} + KA{i} * (c{i} + S{i}/d{i}*(tau(i) - S{i}'*(pA{i} + IA{i}*c{i})))];
        end
        
    else
        
        % ugly repetition of code because matlab does not support zero
        % indexing
        
        if size( KA{i}, 1) > 0
            %PV
            FS = KA{i}*S{i};
            KA0 = [KA0;  (KA{i} - FS/d{i}*U{i}')*Xup{i}];
            if ~bfgs
                zero_off_diag_LM = csX(size(LA0, 1), size(KA{i}, 1));
                LA0 = [LA0, zero_off_diag_LM; zero_off_diag_LM', LA{i} + FS/d{i}*(FS')];
                if strcmp(class(cs), 'casadi.SX')
                    LA0 = casadi_symmetric(LA0);
                end
            else
                zero_off_diag_LM = csX(size(HA0, 1), size(KA{i}, 1));
                inv_vec = HA{i}*FS;
                H_prop = HA{i} - (inv_vec/(d{i} + FS'*inv_vec))*inv_vec';
                HA0 = [HA0, zero_off_diag_LM; zero_off_diag_LM', H_prop];
                if strcmp(class(cs), 'casadi.SX')
                    HA0 = casadi_symmetric(HA0);
                end
            end
            lA0 = [lA0; lA{i} + KA{i} * (c{i} + S{i}/d{i}*(tau(i) - S{i}'*(pA{i} + IA{i}*c{i})))];
        end
    end
end

if ~isempty(lA0)
    b = lA0 - KA0*(a_grav);
    
    if ~bfgs
        if strcmp(class(cs), 'casadi.MX')
            nu = solve(LA0 + cs.eye(size(LA0,1))*1e-12, b, 'ldl');
        else
            Lchol_osim = cholesky(LA0 + cs.eye(size(LA0,1))*1e-6);
            nu = back_sub(Lchol_osim', forward_sub(Lchol_osim,b));
            
%             invosim_fun = Function('f_osim', {q}, {Lchol_osim});
%             invosim_fun.generate(strcat(strcat('osim_iiwa_6D_', num2str(size(KA0,1))), '.c'), struct('with_header', true))
            
        end
    else
        nu = HA0*b;
        
        
%         invosim_fun = Function('f_osim', {q}, {HA0});
%         invosim_fun.generate('osim_iiwa_6D_.c', struct('with_header', true))
    end
else
    nu = 0;
end

con_torque = csX(model.NB, 1);

if ~isempty(lA0)
    lambda{1} = nu;
    con_torque(1) = -S{1}' * KA{1}' * lambda{1};
end
for i = 1:model.NB
    
    if model.parent(i) == 0
        a{i} = Xup{i} * -a_grav + c{i};
    else
        a{i} = Xup{i} * a{model.parent(i)} + c{i};
    end
    if size(KA{i}, 2) > 0 && i > 1
        mi = size(KA{i}, 1);
        lambda{i} = lambda{model.parent(i)}(end - mi + 1 : end);
        lambda{model.parent(i)} = lambda{model.parent(i)}(1 : end - m_i);
        con_torque(i) = -S{i}' * KA{i}' * lambda{i};
    elseif i > 1
        con_torque(i) = 0;
    end
    qdd(i,1) = (u{i} - U{i}'*a{i} + con_torque(i))/d{i};
    a{i} = a{i} + S{i}*qdd(i);
end

a_ee = Xa{n}\a{n};

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


ee_acc_f_ext_fun = 0;
ee_acc_X_fun = 0;
q_acc_X_fun = 0;
tau_X_fun = 0;

% if nargin >= 5 && size(f_ext,1) >0
%     ee_acc_f_ext_fun = jacobian(a{n}, f_ext{n});
%
% end

% if size(KA{1}, 1) > 0
%     q_acc_X_fun = jacobian([qdd;con_torque; nu], [tau; q]);
%     tau_X_fun = q_acc_X_fun(n+1:end,:);
%     q_acc_X_fun = q_acc_X_fun(1:n,:);
% end



