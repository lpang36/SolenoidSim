%Parameters
theta = 74/180*pi; %Relative angle of two coils in xy-plane
xshift = 0.09; %Distance between starts of each coil, x direction
yshift = 0.07; %Distance between starts of each coil, y direction
L = 0.02; %Coil length, wound
n = 320; %Number of coils
R = 0.069; %Coil radius
t1 = [0:1:n*2*pi]; %Parameter of first coil
t2 = [0:1:n*2*pi]; %Parameter of second coil
mu = 4*pi*10.^-7; %Magnetic constant
I = 1; %Current in coil
const = L/(n*2*pi);

%Find first coil points
x1 = const*t1;
y1 = R*cos(t1);
z1 = R*sin(t1);
%Find first coil position differential
dx1 = const;
dy1 = -R*sin(t1);
dz1 = R*cos(t1);
%Initialize second coil points
x2 = const*t2;
y2 = R*cos(t2);
z2 = R*sin(t2);
%Rotate second coil using rotation matrix
x2 = x2*cos(theta)-y2*sin(theta);
y2 = x2*sin(theta)+y2*cos(theta);
x2 = x2+xshift; %Shift in x direction
y2 = y2+yshift; %Shift in y direction
%Find second coil position differential
dx2 = const*cos(theta)-R*sin(t2)*(-sin(theta));
dy2 = const*sin(theta)-R*sin(t2)*cos(theta);
dz2 = R*cos(t2);

%Display coils
scatter3(x1,y1,z1,'filled')
hold on
scatter3(x2,y2,z2,'filled')
%Iterate through points
Px = [-0.05:0.03:L+xshift+0.05];
Py = [-0.05:0.03:L+yshift+0.05];
Pz = [-R-0.02:0.03:R+0.02];
Bfinal = [];
Pfinal = [];
for a = 1:numel(Px)
	for b = 1:numel(Py)
		for c = 1:numel(Pz)
			%Compute B field
			P = [Px(a),Py(b),Pz(c)];
			B = [0,0,0];
			for i=1:numel(t1) %First coil B field
				r = [x1(i),y1(i),z1(i)]-P; %Position vector
				D = sqrt(r(1).^2+r(2).^2+r(3).^2); %Distance between point and coil
				dB = mu*I/4/pi/D.^3*cross([dx1,dy1(i),dz1(i)],r); %Biot-Savart law
				B = B+dB;
			end
			for i=1:numel(t2) %Second coil B field
				r = [x2(i),y2(i),z2(i)]-P; %Position vector
				D = sqrt(r(1).^2+r(2).^2+r(3).^2); %Distance between point and coil
				dB = mu*I/4/pi/D.^3*cross([dx2(i),dy2(i),dz2(i)],r); %Biot-Savart law
				B = B+dB; 
			end
      		Pfinal = cat(1,Pfinal,P);
      		Bfinal = cat(1,Bfinal,B);
		end
	end
end
%Display B field
quiver3(Pfinal(:,1),Pfinal(:,2),Pfinal(:,3),5*Bfinal(:,1),5*Bfinal(:,2),5*Bfinal(:,3),'Color','red','Autoscale','off')