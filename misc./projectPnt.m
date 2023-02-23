function [Q, dist] = projectPnt(pnt, V2, V1)

v = (V2- V1)/norm(V2-V1);
Q = dot(pnt-V1, v)*v + V1;

dist = norm(pnt-Q);


end