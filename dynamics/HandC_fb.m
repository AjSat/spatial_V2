function  [H,C, X_up_all, v, avp, a_grav_links] = HandC_fb( model, x_fb, q, qd, f_ext )

% Adapted Featherstone's HandC function for floating-base robots. Assumes
% that the first link in the model is the floating base.



model.parent(1) = model.NB + 1;
import casadi.*

a_grav = get_gravity(model);
q = [0; q];
qd = [0; qd];

qn = x_fb(1:4);
r = x_fb(5:7);
Xup_fb = plux(rq(qn), r);
Xup{1} = SX.eye(6);

v_fb = x_fb(8:end);
v{1} = Xup{1}*v_fb;
S{1} = SX.eye(6);

avp{1} = Xup_fb * -a_grav;
fvp{1} = model.I{1}*avp{1} + crf(v{1})*(model.I{1}*v{1});
a_grav_links{1} = Xup_fb * -a_grav;

C = SX(model.NB + 5, 1);
H = SX(model.NB + 5, model.NB + 5);

for i = 2:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
  vJ = S{i}*qd(i);
  Xup{i} = XJ * model.Xtree{i};
  v{i} = Xup{i}*v{model.parent(i)} + vJ;
  avp{i} = Xup{i}*avp{model.parent(i)} + crm(v{i})*vJ;
  a_grav_links{i} = Xup{i} * a_grav_links{model.parent(i)};
  fvp{i} = model.I{i}*avp{i} + crf(v{i})*(model.I{i}*v{i});
end

if nargin == 5
  fvp = apply_external_forces( model.parent, Xup, fvp, f_ext );
end

for i = model.NB:-1:2
  C(i+5,1) = S{i}' * fvp{i};
  fvp{model.parent(i)} = fvp{model.parent(i)} + Xup{i}'*fvp{i};
end
fvp{model.NB + 1} = Xup{1}'*fvp{1};
C(1:6, 1) = S{1}'*fvp{model.NB + 1};

IC = model.I;				% composite inertia calculation
for i = model.NB:-1:2
    IC{model.parent(i)} = casadi_symmetric(IC{model.parent(i)} + Xup{i}'*IC{i}*Xup{i});
end
IC{model.NB + 1} = casadi_symmetric(Xup{1}'*IC{1}*Xup{1});

for i = 2:model.NB
  fh = IC{i} * S{i};
  H(i+5,i+5) = S{i}' * fh;
  j = i;
  while model.parent(j) ~= model.NB + 1
    fh = Xup{j}' * fh;
    j = model.parent(j);
    if j > 1
        H(i+5,j+5) = S{j}' * fh;
        H(j+5,i+5) = H(i+5,j+5);
    end
  end
  fh = Xup{1}' * fh;
  H(1:6,i+5) = S{1}' * fh;
  H(i+5,1:6) = H(1:6,i+5);
end

X_up_all = [];
for i = 1:model.NB
    X_up_all = [X_up_all, Xup{i}];
end

H(1:6, 1:6) = IC{model.NB + 1};

