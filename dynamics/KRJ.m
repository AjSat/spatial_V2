function Omega = KRJ(model,q, x_fb, K_con)
% The function takes the robot model, current joint position and Cartesian
%constraints as input and returns the OSIM

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

n = size(q, 1);
qn = x_fb(1:4);				% unit quaternion fixed-->f.b.
r = x_fb(5:7);				% position of f.b. origin
Xup{1} = plux( rq(qn), r );		% xform fixed --> f.b. coords


IA{1} = model.I{1};
KA{1} = [];
LA{1} = [];
lA{1} = [];


for i = 2:model.NB
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
    Xup{i} = XJ * model.Xtree{i};
    IA{i} = model.I{i};
    
end

IA{1} = model.I{1};

nu_vals = [];
counter = 0;
for i = 1:model.NB
    KA{i} = K_con{i};
    m_i = size(K_con{i}, 1);
    LA{i} = csX(m_i, m_i);
    
    if m_i > 0
        counter = counter + size(LA{1}, 1);
        nu_vals = [nu_vals, counter];
        
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


for i = model.NB:-1:2
    U{i} = IA{i} * S{i};
    d{i} = S{i}' * U{i};
    par_i = model.parent(i);
    
    P{i} = SX.eye(6) - (U{i}/d{i})*S{i}';
    
    
    
    IA{model.parent(i)} = IA{model.parent(i)} + Xup{i}' * P{i}*IA{i} * Xup{i};
    if strcmp(class(cs), 'casadi.SX')
        IA{model.parent(i)} = casadi_symmetric(IA{model.parent(i)});
    end   
end

Omega = {};

for i = 2:model.NB
   Omega{i,i} = casadi_symmetric(S{i}/d{i}*S{i}');
end
Omega{1,1} = casadi_symmetric(inv(IA{1}));

for i = 2:model.NB
    j = model.parent(i);
    Omega{i,i} = Omega{i, i} + P{i}'*Xup{i}*Omega{j, j}*Xup{i}'*P{i};
    Omega{i,i} = casadi_symmetric(Omega{i,i});
end

for i = 1:model.NB-1
    for j = i+1:model.NB
        
        if ancest(model, i, j)
            Omega{i,j} = Omega{i, model.parent(j)}*Xup{j}'*P{j};
        else
            Omega{i, j} = P{i}'*Xup{i}*Omega{model.parent(i), model.parent(j)}*Xup{j}'*P{j};
        end
        
        %Omega{i,j} = casadi_symmetric(Omega{i, j});
        Omega{j, i} = Omega{i, j}';



    end
end

M = Omega;

end

function [anc] = ancest(model, i, j)
    
    % return whether i is j's ancestor
    anc = false;
    while model.parent(j) > 0
       if model.parent(j) == i
           anc = true;
           break;
       else 
           j = model.parent(j);
       end
    end

end

