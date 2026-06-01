%Returns the length dependence factor given the
%contractile element length LC - Gaussian model

function out=f1(LC,W,Lo)
%W is the width parameter, Lo the optimal length.
out = exp(-((LC-Lo)/(W*Lo)).^2);
