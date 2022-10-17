function  [X, v] = FKfb( model, xfb, q, qd)

% ver.

import casadi.*

a_grav = get_gravity(model);

qn = xfb(1:4);				% unit quaternion fixed-->f.b.
r = xfb(5:7);				% position of f.b. origin
Xup{6} = plux( rq(qn), r );		% xform fixed --> f.b. coords

vfb = xfb(8:end);
v{6} = Xup{6} * vfb;			% f.b. vel in f.b. coords

IA{6} = model.I{6};
pA{6} = crf(v{6}) * model.I{6} * v{6};

X{6} = Xup{6};

for i = 7:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i-6) );
  vJ = S{i}*qd(i-6);
  Xup{i} = XJ * model.Xtree{i};
  X{i} = Xup{i}*X{model.parent(i)};
  v{i} = Xup{i}*v{model.parent(i)} + vJ;
end

% for i = 7:model.NB
%    v{i} = Xup{i} *v{i}; 
% end

