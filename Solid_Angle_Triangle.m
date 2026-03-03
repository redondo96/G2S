function omega=Solid_Angle_Triangle(p,va,vb,vc)

% Compute relative vectors and precompute norms
r_pa=p-va; norm_r_pa=norm(r_pa);
r_pb=p-vb; norm_r_pb=norm(r_pb);
r_pc=p-vc; norm_r_pc=norm(r_pc);

% Compute cross product and dot products
cross_r_pb_r_pc = cross(r_pb, r_pc);
dot_r_pa_r_pb = dot(r_pa, r_pb);
dot_r_pa_r_pc = dot(r_pa, r_pc);
dot_r_pb_r_pc = dot(r_pb, r_pc);

numerator = dot(r_pa,cross_r_pb_r_pc);
denominator = norm_r_pa * norm_r_pb * norm_r_pc + ...
              dot_r_pa_r_pb * norm_r_pc + ...
              dot_r_pa_r_pc * norm_r_pb + ...
              dot_r_pb_r_pc * norm_r_pa;

omega = 2 * atan2(numerator,denominator);

end