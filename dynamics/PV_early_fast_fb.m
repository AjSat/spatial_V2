% Implements the PV early variant where the dual Hessian is not propagated
% and is eliminated at every joint.

function  [qdd, nu, xd_fb, a_ee] = PV_early_fast_fb( model, x_fb, q, qd, tau, f_ext, K_con, k_con, Soft)


import casadi.*;


if strcmp(class(q), 'casadi.MX')
    cs = MX;
    csX = @MX;
else
    cs = SX;
    csX = @SX;
end

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

qdd = cs.zeros(model.NB,1);

if prop_a_grav
    a_grav_links{1} = Xup{1}*a_grav;
end


for i = 2:model.NB
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
    vJ = S{i}*qd(i);
    Xup{i} = XJ * model.Xtree{i};
    
    if model.parent(i) == 0
        v{i} = vJ;
        c{i} = zeros(size(a_grav));
    else
        v{i} = Xup{i}*v{model.parent(i)} + vJ;
        c{i} = crm(v{i}) * vJ;
    end
    
    if prop_a_grav
            a_grav_links{i} = Xup{i}*a_grav_links{model.parent(i)};
    end
    
    IA{i} = model.I{i};
    pA{i} = crf(v{i}) * (model.I{i} * v{i});
    
    if nargin > 8 && ~isempty(Soft{i})
        IA{i} = IA{i} + Soft{i}.Ki'*Soft{i}.Ri*Soft{i}.Ki;
        if strcmp(class(cs), 'casadi.SX')
            IA{i} = casadi_symmetric(IA{i});
        end
        if prop_a_grav
            Soft{i}.ki = Soft{i}.ki - Soft{i}.Ki*a_grav_links{i};
        end
        pA{i} = pA{i} + Soft{i}.Ki'*Soft{i}.Ri*Soft{i}.ki;
    end
    
end

counter = {};

constrained_links = [];

for i = 1:model.NB
    KA{i} = K_con{i};
    m_i = size(K_con{i}, 1);
    LA{i} = csX(m_i, m_i);
    lA{i} = -k_con{i};
    
    if m_i > 0
        counter{i} = size(KA{i}, 1);
        constrained_links = [constrained_links; i];
    else
        counter{i} = -1;
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

KA0 = [];
LA0 = [];
lA0 = [];

% if size(KA{model.NB}, 1) > 0
% % KA{model.NB} = KA{model.NB}*Xa{n}';
% a_grav_ee = Xa{n} * (-a_grav);
% lA{model.NB} = lA{model.NB} + KA{model.NB} * (-a_grav_ee);
% end

for i = model.NB:-1:2
    
%     if size( KA{i}, 1) > 0 && counter{i} == 0            
%             early_Lchol{i} = cholesky(LA{i} + 1e-6*SX.eye(size(KA{i},1)));
%             IA{i} = IA{i} + KA{i}'*(back_sub(early_Lchol{i}', forward_sub(early_Lchol{i}, KA{i})));
%             pA{i} = pA{i} + KA{i}'*(back_sub(early_Lchol{i}', forward_sub(early_Lchol{i}, lA{i})));
%     end
    
    
    U{i} = IA{i} * S{i};
    d{i} = S{i}' * U{i};
      
    u{i} = tau(i) - S{i}'*pA{i};
    
    Ia = IA{i} - U{i}/d{i}*U{i}';
    if strcmp(class(cs), 'casadi.SX')
        Ia = casadi_symmetric(Ia);
    end
    pa = pA{i} + Ia*c{i} + U{i} * (u{i}/d{i});
    
    if model.parent(i) ~= 0
        
        if size( KA{i}, 1) > 0
            %PV
            FS = KA{i}*S{i};
            FS_norm_sqr = FS'*FS;
            usvd = FS/sqrt(FS_norm_sqr);
            Sigma{i} = FS_norm_sqr/d{i};
            e1 = csX(size(KA{i}, 1), 1);
            e1(1) = 1;
            w{i} = usvd + (usvd(1)/abs(usvd(1))*sqrt(usvd'*usvd))*e1;
            wi_norm_sqr{i} = w{i}'*w{i};
            
            KAnew = KA{i} - FS*(U{i}'/d{i});
            lAnew = lA{i} + KA{i} * (c{i} + S{i}/d{i}*(tau(i) - S{i}'*(pA{i} + IA{i}*c{i})));
            
            Kitemp{i} = KAnew - 2*(w{i}/wi_norm_sqr{i})*(w{i}'*KAnew);
            litemp{i} = lAnew - 2*(w{i}/wi_norm_sqr{i})*(w{i}'*lAnew);
           
%             LA{i} =  FS/d{i}*(FS'); means FS/sqrt(d(i)) is the rank-1
%             addition matrix.
            KA{model.parent(i)} = [KA{model.parent(i)};  (Kitemp{i}(2:end,:))*Xup{i}];
            lA{model.parent(i)} = [lA{model.parent(i)}; litemp{i}(2:end)];
            
            Ia = Ia + (Kitemp{i}(1,:)'*Kitemp{i}(1,:))/Sigma{i};
            pa = pa + (litemp{i}(1)/Sigma{i})*Kitemp{i}(1,:)';
        end
        
        IA{model.parent(i)} = IA{model.parent(i)} + Xup{i}' * Ia * Xup{i};
        if strcmp(class(cs), 'casadi.SX')
            IA{model.parent(i)} = casadi_symmetric(IA{model.parent(i)});
        end
        pA{model.parent(i)} = pA{model.parent(i)} + Xup{i}' * pa;        
%     else
%         
%         % ugly repetition of code because matlab does not support zero
%         % indexing        
% %         if size( KA{i}, 1) > 0 && counter{i} > 0
% %             %PV
% %             FS = KA{i}*S{i};
% %             KA0 = [KA0;  (KA{i} - FS/d{i}*U{i}')*Xup{i}];
% % 
% %             zero_off_diag_LM = csX(size(LA0, 1), size(KA{i}, 1));
% %             LA0 = [LA0, zero_off_diag_LM; zero_off_diag_LM', LA{i} + FS/d{i}*(FS')];
% %             if strcmp(class(cs), 'casadi.SX')
% %                 LA0 = casadi_symmetric(LA0);
% %             end
% %             lA0 = [lA0; lA{i} + KA{i} * (c{i} + S{i}/d{i}*(tau(i) - S{i}'*(pA{i} + IA{i}*c{i})))];
% %         end
    end
end

if size(lA{1}, 1) > 0
%     b = lA0 - KA0*(a_grav);   
%         if strcmp(class(cs), 'casadi.MX')
%             nu = solve(LA0 + cs.eye(size(LA0,1))*1e-12, b, 'ldl');
%         else
%             Lchol_osim = cholesky(LA0 + cs.eye(size(LA0,1))*1e-8);
%             nu = back_sub(Lchol_osim', forward_sub(Lchol_osim,b));
%         end
else
    Lchol = cholesky(IA{1});
    a{1} = -back_sub(Lchol', forward_sub(Lchol,pA{1}));
    nu = 0;
%     a_fun = Function('f_rob_dyn', {x_fb, q(2:end), qd(2:end), tau(2:end)}, {a{1}});
%     a_fun.n_instructions
end
lambda{model.NB} = [];
con_torque = csX(model.NB, 1);
for i = 2:model.NB
    

    a{i} = Xup{i} * a{model.parent(i)} + c{i};
    if size(KA{i}, 1) > 0  
        mi = size(KA{i}, 1);
        lambda_opt = (Kitemp{i}(1,:)*(a{i} - c{i}) + litemp{i}(1))/Sigma{i};
        if mi >  1
            lambda_concat = [lambda_opt; lambda{model.parent(i)}(end - mi + 2 : end)];
            lambda{model.parent(i)} = lambda{model.parent(i)}(1 : end - mi + 1);
        else
            lambda_concat = lambda_opt;
        end
        lambda{i} = lambda_concat - 2*(w{i}/wi_norm_sqr{i})*(w{i}'*lambda_concat);
        con_torque(i) = -S{i}' * KA{i}' * lambda{i};
    elseif i > 1
        con_torque(i) = 0;
    end
    qdd(i,1) = (u{i} - U{i}'*a{i} + con_torque(i))/d{i};
    a{i} = a{i} + S{i}*qdd(i);
end

nu = [];
for j = constrained_links
    nu = [nu, lambda{j}];
end
nu = nu';

a_ee = Xa{model.NB}\a{model.NB};

%In world frame final segment has acceleration
%final_seg = Xa{n} *a{n};

%tau_apply = SX.zeros(model.NB, 1);
% for i = 1:model.NB
%     tau_apply(i) = con_torque(i) + tau(i);
% end
a1_fb = Xup{1}\a{1} + a_grav;
xd_fb = [x_fb(8:end); a1_fb];


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



