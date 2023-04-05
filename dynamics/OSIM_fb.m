% Compute the inv_OSIM for a floating-base robot

function  [inv_OSIM, IA, KA, LA] = OSIM_fb( model, x_fb, q, K_con, k_con)

import casadi.*;

if strcmp(class(q), 'casadi.MX')
    cs = MX;
    csX = @MX;
else
    cs = SX;
    csX = @SX;
end


q = [0; q];

n = size(q, 1);

qn = x_fb(1:4);				% unit quaternion fixed-->f.b.
r = x_fb(5:7);				% position of f.b. origin
Xup{1} = plux( rq(qn), r );		% xform fixed --> f.b. coords


IA{1} = model.I{1};
KA{1} = [];
LA{1} = [];

for i = 2:model.NB
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
    Xup{i} = XJ * model.Xtree{i};
    IA{i} = model.I{i};   
end


for i = 1:model.NB
    KA{i} = K_con{i};
    m_i = size(K_con{i}, 1);
    LA{i} = csX(m_i, m_i);    
end

parent = model.parent;
    for i = 1:length(parent)
        if parent(i) == 0
            Xa{i} = Xup{i};
        else
            Xa{i} = Xup{i} * Xa{parent(i)};
        end
    end

for i = model.NB:-1:2
    U{i} = IA{i} * S{i};
    d{i} = S{i}' * U{i};
    
    Ia = IA{i} - (U{i}/d{i})*U{i}';
    if strcmp(class(cs), 'casadi.SX')
        Ia = casadi_symmetric(Ia);
    end

        IA{model.parent(i)} = IA{model.parent(i)} + Xup{i}' * Ia * Xup{i};
        if strcmp(class(cs), 'casadi.SX')
            IA{model.parent(i)} = casadi_symmetric(IA{model.parent(i)});
        end     
        
        if size( KA{i}, 1) > 0
            %PV          
            FS = KA{i}*S{i};
            KA{model.parent(i)} = [KA{model.parent(i)};  (KA{i} - FS/d{i}*U{i}')*Xup{i}];
            zero_off_diag_LM = csX(size(LA{model.parent(i)}, 1), size(KA{i}, 1));
            LA{model.parent(i)} = [LA{model.parent(i)}, zero_off_diag_LM; zero_off_diag_LM', LA{i} + FS/d{i}*(FS')];
            if strcmp(class(cs), 'casadi.SX')
                LA{model.parent(i)} = casadi_symmetric(LA{model.parent(i)});
            end
        end
    
end

if isempty(LA{1})
    inv_OSIM = [];
else
    inv_OSIM = LA{1} + KA{1}*inv(IA{1})*KA{1}';
    inv_OSIM = casadi_symmetric(inv_OSIM);
    %invosim_fun = Function('f_osim', {q(2:end)}, {casadi_symmetric(inv(full_osim))});
    
    % Uncomment to create a function
%     invosim_fun = Function('f_osim', {q(2:end)}, {cholesky(inv_OSIM)});
%     invosim_fun.n_instructions
    %invosim_fun.generate(strcat(strcat('osiminv_atlas_18D_', num2str(size(LA{1},1))), '.c'), struct('with_header', true))
end

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



