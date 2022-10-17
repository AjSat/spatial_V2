function  [tau, vb, ab, zwang] = ID_bf( model, q, qd, qdd, f_ext, qdd_free )

import casadi.*

%Modifying Featherstone code below to implement ID in base-frame rather
%than the body frame.

% ID  Inverse Dynamics via Recursive Newton-Euler Algorithm
% ID(model,q,qd,qdd,f_ext) calculates the inverse dynamics of a kinematic
% tree via the recursive Newton-Euler algorithm.  q, qd and qdd are vectors
% of joint position, velocity and acceleration variables; and the return
% value is a vector of joint force variables.  f_ext is an optional
% argument specifying the external forces acting on the bodies.  It can be
% omitted if there are no external forces.  The format of f_ext is
% explained in the source code of apply_external_forces.

a_grav = get_gravity(model);
 
tau = SX(model.NB,1);

zwang = 0;

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
  
  vJ = S{i}*qd(i);
  Xup{i} = sparsify(XJ * model.Xtree{i});
  if model.parent(i) == 0
    v{i} = vJ;
    a{i} = Xup{i}*(-a_grav) + S{i}*qdd(i);
    if nargin == 6
        ab_free{i} = Xup{i}*(-a_grav) + S{i}*qdd_free(i);
    end
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    a{i} = Xup{i}*a{model.parent(i)} + S{i}*qdd(i) + crm(v{i})*vJ;
    if nargin == 6
        ab_free{i} = Xup{i}*ab_free{model.parent(i)} + S{i}*qdd(i) + crm(v{i})*vJ;
    end
  end
  f{i} = model.I{i}*a{i} + sparsify(crf(v{i}))*model.I{i}*v{i};
  
  if nargin == 6
      zwang = zwang + (a{i} - ab_free{i})' *  model.I{i} * (a{i} - ab_free{i});
  end
end

if nargin == 5 && nargin ~= 6
  f = apply_external_forces( model.parent, Xup, f, f_ext );
end

for i = model.NB:-1:1
  tau(i,1) = S{i}' * f{i};
  if model.parent(i) ~= 0
    f{model.parent(i)} = f{model.parent(i)} + Xup{i}'*f{i};
  end
end

vb = cell_to_mat(v);
ab = cell_to_mat(a);
