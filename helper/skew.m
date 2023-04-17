
function skew = skew(u)
skew = [
    0 -u(3) u(2)
    u(3) 0 -u(1)
    -u(2) u(1) 0
    ];
end


% same as wedge function