function [Omega] = EFP(model, q, x_fb, K_con, gca)
% Implements Patrick Wensing's EFP algorithm to compute the OSIM matrix
%   

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


for i = 2:model.NB
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
    Xup{i} = XJ * model.Xtree{i};
    IA{i} = model.I{i};
    
end

IA{1} = model.I{1};


parent = model.parent;
for i = 1:length(parent)
    if parent(i) == 0
        Xa{i} = Xup{i};
    else
        Xa{i} = Xup{i} * Xa{parent(i)};
    end
    
    KA{i} = [];
    cons{i} = [];
end

for i = 1:model.NB
    KA{i, i} = K_con{i};
    m_i = size(K_con{i}, 1);
    
    if m_i > 0
        cons{i} = [i];
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
    for j = cons{i}
        KA{par_i, j} = KA{i, j}*P{i}'*Xup{i};
        cons{par_i} = [cons{par_i}, j];
    end
end

Omega = {};

for i = 2:model.NB
   Omega{i,i} = casadi_symmetric(S{i}/d{i}*S{i}');
end
Omega{1,1} = casadi_symmetric(inv(IA{1}));

for j = cons{1}
    Omega{1, j} = Omega{1,1}*KA{1,j}';
end

for i = 2:model.NB
    par_i = model.parent(i);
    for j = cons{i}
        Omega{i,j} = Omega{i, i}*KA{i,j}' + P{i}'*Xup{i}*Omega{par_i, j};
    end
end

cons_1 = fliplr(cons{1});
for i = 1:length(cons_1)-1
    for j = i+1:length(cons_1)
        Omega{cons_1(i), cons_1(j)} = KA{gca(i, j), cons_1(i)} * Omega{gca(i, j), cons_1(j)};
        Omega{cons_1(j), cons_1(i)} = (Omega{cons_1(i), cons_1(j)})';
    end
end

% for i = 1:model.NB-1
%     for j = i+1:model.NB
%         
%         if ancest(model, i, j)
%             Omega{i,j} = Omega{i, model.parent(j)}*Xup{j}'*P{j};
%         else
%             Omega{i, j} = P{i}'*Xup{i}*Omega{model.parent(i), model.parent(j)}*Xup{j}'*P{j};
%         end
%         
%         %Omega{i,j} = casadi_symmetric(Omega{i, j});
%         Omega{j, i} = Omega{i, j}';
% 
% 
% 
%     end
% end

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
    
    Omega = [];
end

