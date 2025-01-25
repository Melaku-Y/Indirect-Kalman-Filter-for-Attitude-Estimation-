function [M] = quaternion2dcm_student(q) 
% Question 1 from equation_1

% Input: q is a 4x1 quaternion [q1; q2; q3; q4]
% Output: DCM is the 3x3 direction cosine matrix
q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);

M = [-1+2*(q0^2+q1^2), 2*(q1*q2+q0*q3), 2*(q1*q3-q0*q2);
         2*(q1*q2-q0*q3), -1+2*(q0^2+q2^2), 2*(q2*q3+q0*q1);
         2*(q1*q3+ q0*q2), 2*(q2*q3-q0*q1), -1+2*(q0^2+q3^2)];
end

