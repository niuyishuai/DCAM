function draw3d(f,interval)
figure;
x=interval(1):0.1:interval(2);
y=interval(3):0.1:interval(4);
z=zeros(length(x),length(y));
for i=1:length(x)
    for j=1:length(y)
        z(j,i)=f.eval([x(i);y(j)]);
    end
end
surf(x,y,z);
title('dc function');
xlabel('x_1');
ylabel('x_2');
zlabel('f(x)');
end