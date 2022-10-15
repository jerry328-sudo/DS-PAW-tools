x = linspace(0,3*pi,200);
y = cos(x) + rand(1,200);
c = linspace(1,10,length(x));
%这样RGB值就是变量了
scatter(x,y,[],c)