function M = SAR_PointingVectors( u_phaseCenter)
% M = [ los | right_crossrange | (los x rc) downward ]
% unit normal to LOS & cross range
u_vertical = [0 0 1];
u_test = cross(u_phaseCenter,u_vertical);
u_crossRange = u_test ./ norm(u_test);
u_phaseCenterNormal = cross(u_phaseCenter,u_crossRange);
u_phaseCenterNormal = u_phaseCenterNormal ./ norm(u_phaseCenterNormal);

M = [u_phaseCenter; u_crossRange;  u_phaseCenterNormal;];


% =====================================
% M (row vectors) represents the unit vectors of the measurement frame
% in the basis vector set
% Left multiply of M rotates "from" basis "to" measurement frame
% =====================================

end