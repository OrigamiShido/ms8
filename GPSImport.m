function sc=GPSImport(sc,orbit)
%计算真近地点
e=orbit(1);
M=orbit(5);
fun=@(E) E-e*sind(E)-M;
E=fzero(fun,M);
theta=atand((sqrt(1-e^2)*sind(E))/(cosd(E)-e));
satellite(sc,orbit(2),e,orbit(6),orbit(3),orbit(4),theta);

end
