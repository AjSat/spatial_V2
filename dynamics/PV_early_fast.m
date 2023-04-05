% Implements the PV early variant where the dual Hessian is not propagated
% and is eliminated at every joint.

function  [qdd, nu, a_ee, Xee] = PV_early_fast( model, q, qd, tau, f_ext, K_con, k_con, Soft)


import casadi.*;


if strcmp(class(q), 'casadi.MX')
    cs = MX;
    csX = @MX;
else
    cs = SX;
    csX = @SX;
end

n = size(q, 1);
a_grav = get_gravity(model);
prop_a_grav = true;
qdd = cs.zeros(model.NB,1);

model.parent(1) = model.NB+1;


for i = 1:model.NB
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
    vJ = S{i}*qd(i);
    Xup{i} = XJ * model.Xtree{i};
    
    if model.parent(i) == model.NB+1
        v{i} = vJ;
        c{i} = zeros(size(a_grav));
        if prop_a_grav
            a_grav_links{i} = Xup{i}*a_grav;
        end
    else
        v{i} = Xup{i}*v{model.parent(i)} + vJ;
        c{i} = crm(v{i}) * vJ;
        if prop_a_grav
            a_grav_links{i} = Xup{i}*a_grav_links{model.parent(i)};
        end
    end
    
    IA{i} = model.I{i};
    pA{i} = crf(v{i}) * (model.I{i} * v{i});
    
    if nargin > 8 && ~isempty(Soft{i})
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
    LA{i} = csX(m_i, m_i);
    lA{i} = -k_con{i};
    if m_i > 0 && prop_a_grav
        lA{i} = lA{i} + KA{i}*a_grav_links{i};
    end
    
end

parent = model.parent;
for i = 1:length(parent)
    if parent(i) == model.NB + 1
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

for i = model.NB:-1:1    
    
    U{i} = IA{i} * S{i};
    d{i} = S{i}' * U{i};
    u{i} = tau(i) - S{i}'*pA{i};
    
    Ia = IA{i} - U{i}/d{i}*U{i}';
    if strcmp(class(cs), 'casadi.SX')
        Ia = casadi_symmetric(Ia);
    end
    pa = pA{i} + Ia*c{i} + U{i} * (u{i}/d{i});
    
    
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
    
    if model.parent(i) ~= model.NB + 1
        
        
        IA{model.parent(i)} = IA{model.parent(i)} + Xup{i}' * Ia * Xup{i};
        if strcmp(class(cs), 'casadi.SX')
            IA{model.parent(i)} = casadi_symmetric(IA{model.parent(i)});
        end
        pA{model.parent(i)} = pA{model.parent(i)} + Xup{i}' * pa;
    end

end

%todo: check if there is an unresolved constraint at the base and raise
%error.

a{model.NB + 1} = -a_grav;

lambda{model.NB} = [];
con_torque = csX(model.NB, 1);
for i = 1:model.NB
    
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

nu = lambda{model.NB}; %todo: fix this to get multipliers from all branches.
a_ee = Xa{model.NB}\a{model.NB};
Xee = Xa{model.NB};