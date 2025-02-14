function p_u = interp2grid_toy(u,Coord_1,Coord_2)
%function p_u = interp2grid_toy(u,Coord_1,Coord_2)
%INPUT:
%u:         vector. function values on Coordinate 1.
%Coord_1:   vector or structure. Coodinate 1 that u is on.   
%Coord_2:   vector or structure. Coodinate 2 that u is interpolated to.
%
%OUTPUT:
%p_u:       vector. function values u is interpolated to Coordinate 2.
%
% Feb-14,2025.

if  isstruct(Coord_1) % if Coord_1 is a mesh structure
    Coord_1 = Coord_1.Coordinates; % get the x, y coordinate
end
if  isstruct(Coord_2) % if Coord_2 is a mesh structure
    Coord_2 = Coord_2.Coordinates; % get the x, y coordinate
end

F = griddedInterpolant(Coord_1,u); % linear,natural
p_u = F(Coord_2);
end 