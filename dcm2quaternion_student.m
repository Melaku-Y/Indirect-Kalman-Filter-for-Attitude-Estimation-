function q = dcm2quaternion_student(C)

f = zeros(4,1);

f(1) = 0.25 * (1 + trace(C));
f(2) = 0.25 * (C(1,1) - C(2,2) - C(3,3) + 1);
f(3) = 0.25 * (C(2,2) - C(1,1) - C(3,3) + 1);
f(4) = 0.25 * (C(3,3) - C(2,2) - C(1,1) + 1);

[maxf,index ] = max(f);

q = zeros(4,1);
if ( index == 1 )
    q(1) = sqrt(f(1));
    q(2) = (C(2,3) - C(3,2)) / (4*q(1));
    q(3) = (C(1,3) - C(3,1)) / (-4*q(1));
    q(4) = (C(1,2) - C(2,1)) / (4*q(1));
elseif ( index == 2 )
    q(2) = sqrt(f(2));
    q(1) = (C(2,3) - C(3,2)) / (4*q(2));
    q(3) = (C(1,2) + C(2,1)) / (4*q(2));
    q(4) = (C(1,3) + C(3,1)) / (4*q(2));
elseif ( index == 3 )
    q(3) = sqrt(f(3));
    q(1) = (C(1,3) - C(3,1)) / (-4*q(3));
    q(2) = (C(1,2) + C(2,1)) / (4*q(3));
    q(4) = (C(2,3) + C(3,2)) / (4*q(3));
elseif ( index == 4 )
    q(4) = sqrt(f(4));
    q(1) = (C(1,2) - C(2,1)) / (4*q(4));
    q(2) = (C(1,3) + C(3,1)) / (4*q(4));
    q(3) = (C(2,3) + C(3,2)) / (4*q(4));
end