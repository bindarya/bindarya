%Лабораторная работа по решению ур-ия теплопроводности в матлаб
%Сандитова Биндарья 7112
%Вариант 9
%зададим константы, отрезки в которых будем рассматривать задачу
t_0 = 0;
t_e = 8;
x_0 = 0;
x_e = 1;
x_steps = 20;
h = (x_e-x_0)/(x_steps-1);

x = h:h:x_e-h;
u_0 = exp(-x.^2).*sin(x);

u = ode45(@eq9, [t_0 t_e], u_0);
solution = u.y';

siz = size(solution);
t = linspace(t_0, t_e, siz(1))';
tau = t(2)-t(1);
% a = 1/4
tau = tau/4;
solution_true = zeros(siz(1), x_steps-2);

for i = 1:siz(1)
    solution_true(i, :) = exact(t(i), x);
end

error =  sqrt(sum((solution-solution_true).^2, 2));

figure(2)
subplot(1, 2, 1)
[X,T]=meshgrid(x,t);

plot3(T, X, solution, 'c');
hold on;
plot3(T, X, solution_true, 'r');
hold on;
xlabel('t')
ylabel('x')
zlabel('u')
title('Точное решение и решение методом прямых')

subplot(1, 2, 2)
plot(t, error, '-');
xlabel('t')
ylabel('error')
title('Ошибка')