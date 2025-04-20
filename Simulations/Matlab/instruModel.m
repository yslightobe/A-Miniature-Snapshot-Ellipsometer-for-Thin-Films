function Iout = instruModel(d1,d2,Iin,Ms,theta,num,configuration)
% 角度， 角度
T = @(theta)[1 0 0 0; 0 cosd(2*theta) -sind(2*theta) 0;
        0 sind(2*theta) cosd(2*theta) 0; 0 0 0 1];
Md = @(d,theta)T(-theta)*[1 0 0 0; 0 1 0 0; 
    0 0 cosd(d) sind(d); 0 0 -sind(d) cosd(d)]*T(theta);
P = @(theta)T(-theta)*0.5*[1 1 0 0; 1 1 0 0; 0 0 0 0; 0 0 0 0]*T(theta);
Sout = zeros(4,num);
for i = 1:num
    Sin = [Iin(i); 0; 0; 0];
    Sout(:,i) = P(configuration.A_azimuth)*Md(d2(1,i),configuration.r2_azimuth)*...
        T(-theta)*Ms(:,:,i)*T(theta)*Md(d1(1,i),configuration.r1_azimuth)*P(configuration.P_azimuth)*Sin;
end
Iout  = Sout(1,:);
end