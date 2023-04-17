
function xi = makeTwist(delta)
% from vector 9x1 to se2(3) 5x5
R_lie = skew(delta(1:3));
v_lie = [delta(4);delta(5);delta(6)];
p_lie = [delta(7);delta(8);delta(9)];
xi = [R_lie, v_lie, p_lie;...
    zeros(1,3),  0,     0;...
    zeros(1,3),  0,     0];
end