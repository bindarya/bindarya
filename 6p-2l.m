%Лабораторная работа по решению ур-ия теплопроводности в матлаб
%Сандитова Биндарья 7112
%Вариант 9
%зададим константы, отрезки в которых будем рассматривать задачу
t_0 = 0;
t_e = 8;
x_0 = 0;
x_e = 1;
%коэффициент температуропроводности
a = 0.25;
%зададим шаг
t_steps = 100;
x_steps = 20;
xhi = 1;
h = (x_e-x_0)/(x_steps-1);
tau =(t_e-t_0)/(t_steps-1);

%подсчитаем число Куранта
disp('Kurant:');
disp((tau*a)/(h^2));
%зададим нач условия
u_c = h:h:x_e-h;
u_c = exp(-u_c.^2).*sin(u_c);
%создадим и заполним матрицу A
A = zeros(x_steps-2, x_steps-2);
for i = 1:x_steps-3
    A(i, i) = h^2+2*xhi*tau*a;
    A(i+1, i) = -xhi*tau*a;
    A(i, i+1) = -xhi*tau*a;
end
A(x_steps-2, x_steps-2) = h^2+2*xhi*tau*a;

x = h:h:(x_e-h);
%создадим пустые массивы, в которые будем ложить полученные решения
solution = zeros(t_steps-1, x_steps-2);
solution_true = zeros(t_steps-1, x_steps-2);
%В цикле по времени найдем значения u
for t = tau:tau:t_e
    f = zeros(1, x_steps-2)';
    for i = 2:x_steps-3
        f(i) = (h^2-2*(1-xhi)*tau*a)*u_c(i)+(1-xhi)*tau*a*u_c(i+1)+(1-xhi)*tau*a*u_c(i-1);
    end
    u_0 = exact(t-tau, x_0);
    u_e = exact(t-tau, x_e);
    u_c0 = exact(t, x_0);
    u_ce = exact(t, x_e);
    f(1) = (h^2-2*(1-xhi)*tau*a)*u_c(1)+(1-xhi)*tau*a*u_c(2)+(1-xhi)*tau*a*u_0+xhi*tau*a*u_c0;
    f(x_steps-2) = (h^2-2*(1-xhi)*tau*a)*u_c(x_steps-2)+(1-xhi)*tau*a*u_e+(1-xhi)*tau*a*u_c(x_steps-3)+xhi*tau*a*u_ce;
    u_c = linsolve(A, f)';
    solution(int16(t/tau), :) = u_c;
    u_true = zeros(1, x_steps-2)';
    for i = 1:x_steps-2
        u_true(i) = exact(t, x(i));
    end
    solution_true(int16(t/tau), :) = u_true;
end
%подсчитаем ошибку
error =  sqrt(sum((solution-solution_true).^2, 2));

t = tau:tau:t_e;
%построим графики
figure(1)
subplot(1, 2, 1)
[X,T]=meshgrid(x,t);

plot3(T, X, solution, 'c');
hold on;
plot3(T, X, solution_true, 'r');
hold on;
xlabel('t')
ylabel('x')
zlabel('u')
title('Неявный метод')

subplot(1, 2, 2)
plot(t, error, '-');
xlabel('t')
ylabel('error')
title('Ошибка')