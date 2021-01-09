function draw2d(f,vars,interval)
figure;
x=interval(1):0.1:interval(2);
y=zeros(size(x));
for i=1:length(x)
    assign(vars,x(i));
    y(i)=value(f);
end
plot(x,y);
title('dc function');
xlabel('x');
ylabel('objectives');
end