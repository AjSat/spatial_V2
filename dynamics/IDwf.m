function  tau = IDwf(model, q, qd, qdd, f_ext )

import casadi.*

% ID  Inverse Dynamics via Recursive Newton-Euler Algorithm
% ID(model,q,qd,qdd,f_ext) calculates the inverse dynamics of a kinematic
% tree via the recursive Newton-Euler algorithm.  q, qd and qdd are vectors
% of joint position, velocity and acceleration variables; and the return
% value is a vector of joint force variables.  f_ext is an optional
% argument specifying the external forces acting on the bodies. Computing 
% this in the base frame rather than the body frame of the robot.

a_grav = get_gravity(model);
 
tau = SX.zeros(model.NB,1);

% T_ee = A_base;
% Ti = cell(n+1,1);
% Ti_mat = SX.zeros(4, (n+1)*4);
% for i = 1:n
%     T_ee = T_ee*A{i};
%     %tranformation due to the joint motion
%     if norm([robot.Bodies{1,i+first_joint_no-1}.Joint.JointAxis] - [0 0 1]) == 0
%     T_j = [cos(q(i)) -sin(q(i)) 0 0; sin(q(i)) cos(q(i)) 0 0; 0 0 1 0; 0 0 0 1];
%     elseif norm([robot.Bodies{1,i+first_joint_no-1}.Joint.JointAxis] - [0 1 0]) == 0
%     %for y-axis joints
%     T_j = [cos(q(i)) 0 sin(q(i)) 0; 0 1 0 0; -sin(q(i)) 0 cos(q(i)) 0; 0 0 0 1]; 
%     else
%         error('unexpected joint axis')
%     end
%     %T_j = axang2rotm([robot.Bodies{1,i+first_joint_no-1}.Joint.JointAxis, q(i)]);
%     T_ee = T_ee*T_j;
%     Ti{i} = T_ee;
%     Ti_mat(:,(i-1)*4+1:i*4) = Ti{i};
% end

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
%   Xtree2{i} = model.Xtree{i};
%   Xtree2{i}(4:6, 1:3) = -Xtree2{i}(1:3, 1:3)'*model.Xtree{i}(4:6, 1:3)*Xtree2{i}(1:3, 1:3)';
%   Xtree2{i}(1:3, 1:3) = Xtree2{i}(1:3, 1:3)';
%   Xtree2{i}(4:6, 4:6) = Xtree2{i}(4:6, 4:6)';
%   Xtree{i} = SX(6,6);
%   Xtree{i}(1:3, 1:3) = model.Tf{i}(1:3,1:3);
%   Xtree{i}(4:6, 4:6) = model.Tf{i}(1:3,1:3);
%   Xtree{i}(4:6, 1:3) = model.Tf{i}(1:3,1:3)*skew(model.Tf{i}(1:3,1:3)'*model.Tf{i}(1:3,4));
%   Xdown{i} = Xtree{i} * XJ';

   
  
  
  if model.parent(i) == 0
    Tj{i} = model.Tf{i};
    Tj{i}(1:3, 1:3) = Tj{i}(1:3, 1:3)*sparsify([cos(q(i)), -sin(q(i)), 0; sin(q(i)), cos(q(i)), 0; 0, 0, 1]);
    Xwf{i} = [Tj{i}(1:3,1:3), SX(3,3); skew(Tj{i}(1:3,4))*Tj{i}(1:3,1:3), Tj{i}(1:3,1:3)];
    S{i} = Xwf{i}*S{i};
    vJ = S{i}*qd(i);
    v{i} = vJ;
    a{i} = -a_grav + S{i}*qdd(i);
  else
    Tj{i} = Tj{model.parent(i)}*model.Tf{i};
    Tj{i}(1:3, 1:3) = Tj{i}(1:3, 1:3)*sparsify([cos(q(i)), -sin(q(i)), 0; sin(q(i)), cos(q(i)), 0; 0, 0, 1]);
    Xwf{i} = [Tj{i}(1:3,1:3), SX(3,3); skew(Tj{i}(1:3,4))*Tj{i}(1:3,1:3), Tj{i}(1:3,1:3)];
    Xwf{i}(4:6, 1:3) = skew(Tj{i}(1:3,4))*Tj{i}(1:3,1:3);
    S{i} = Xwf{i}*S{i};
    vJ = S{i}*qd(i);
    v{i} = v{model.parent(i)} + vJ;
    a{i} = a{model.parent(i)} + S{i}*qdd(i) + crm(v{i})*vJ;
  end
  
  Xwf_inv{i} = SX(6,6);
  Xwf_inv{i}(1:3, 1:3) = Xwf{i}(1:3, 1:3)';
  Xwf_inv{i}(4:6, 4:6) = Xwf{i}(1:3, 1:3)';
  Xwf_inv{i}(4:6, 1:3) = -Xwf{i}(1:3, 1:3)'*skew(Tj{i}(1:3,4));  %-Xwf{i}(1:3, 1:3)'*Xwf{i}(4:6, 1:3)*Xwf{i}(1:3, 1:3)';
  
%   Xwf_conj = SX(6,6);
%   Xwf_conj(1:3, 1:3) = Xwf{i}(1:3, 1:3);
%   Xwf_conj(4:6, 4:6) = Xwf{i}(1:3, 1:3);
%   Xwf_conj(1:3, 4:6) = Xwf{i}(4:6, 1:3);
  
  f{i} = [Xwf{i}(1:3, 1:3), Xwf{i}(4:6, 1:3); SX(3,3), Xwf{i}(1:3, 1:3)]*(model.I{i}*(Xwf_inv{i}*a{i})) + crf(v{i})*([Xwf{i}(1:3, 1:3), Xwf{i}(4:6, 1:3); SX(3,3), Xwf{i}(1:3, 1:3)]*(model.I{i}*(Xwf_inv{i}*v{i})));
end

if nargin == 5
  for i = 1:model.NB
      if size(f_ext{i}) > 0
        f{i} = f{i} + f_ext{i};
      end
  end
end

for i = model.NB:-1:1
  tau(i,1) = S{i}' * f{i};
  if model.parent(i) ~= 0
    f{model.parent(i)} = f{model.parent(i)} + f{i};
  end
end
