%Velocity dependence function (Hill-Katz model)

%u is the contraction velocity= -LCdot

function out=g1(u,Vm,A,gmax)  

% The parameters below are passed to the function, reference values given.

%A = 0.25;			% Hill constant
%gmax = 1.5;         % Maximum ecc to isometric force ratio

for idx = 1:length(u)


if u(idx) < 0
      out(idx)=(A*Vm - A*Vm*gmax + u(idx)*gmax + A*u(idx)*gmax)./(A*Vm + u(idx) + A*u(idx) - A*Vm*gmax);     %CE lengthens (Katz model)
else 
      out(idx) =	(A*Vm - A*u(idx))./(A*Vm + u(idx));       % CE shortens (Hill model)
end
end

out = out';


    