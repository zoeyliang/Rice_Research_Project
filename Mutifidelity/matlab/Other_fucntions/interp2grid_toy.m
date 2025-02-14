function p_u = interp2grid_toy(u,Coord_1,Coord_2)

if  isstruct(Coord_1) % if Coord_1 is a mesh structure
    Coord_1 = Coord_1.Coordinates; % get the x, y coordinate
end
if  isstruct(Coord_2) % if Coord_2 is a mesh structure
    Coord_2 = Coord_2.Coordinates; % get the x, y coordinate
end

F = griddedInterpolant(Coord_1,u); % linear,natural
p_u = F(Coord_2);
end 