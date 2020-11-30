close all;
clear variables;
 
adc0 = load('calibzero.dat');
adc68 =  load('calib.dat');
adc0(:,1) = []
adc68(:,1) = []

pa0 = ones(length(adc0), 1)*0;
pa68 = ones(length(adc68), 1) * 68;

cAdc = [adc0; adc68];
cPa = [pa0; pa68];

c = polyfit(cAdc, cPa, 1);

cFigure = figure('Name', 'Калибровка', 'NumberTitle', 'off'); hold all;

plot(adc0,	pa0, '.', 'MarkerSize',	20);
plot(adc68,	pa68,	'.', 'MarkerSize',	20);

plot(cAdc, polyval(c, cAdc));

legend('0 Па', '68 Па', 'Калибровочная зависимость', 'Location', 'northwest');

grid on;
xlabel('Отсчёты АЦП');
ylabel('ΔP, Па');
title('Калибровка измерительной системы');
text(mean(cAdc) * 0.8, mean(cPa) * 0.3, ['ΔP (adс) = ', num2str(c(1)), ' * adс + ', num2str(c(2)), ' [Па]']);

saveas(cFigure, "Калибровка2.png");



mm01 = load('01mm.dat');
mm11 = load('11mm.dat');
mm21 = load('21mm.dat');
mm31 = load('31mm.dat');
mm41 = load('41mm.dat'); 
mm51 = load('51mm.dat');
mm61 = load('61mm.dat');
mm71 = load('71mm.dat');

dx = 0.25; % MM
x = mm01(:, 1)' * 0.25;

p01 = polyval(c, mm01(:, 2));
p11 = polyval(c, mm11(:, 2));
p21 = polyval(c, mm21(:, 2));
p31 = polyval(c, mm31(:, 2));
p41 = polyval(c, mm41(:, 2));
p51 = polyval(c, mm51(:, 2));
p61 = polyval(c, mm61(:, 2));
p71 = polyval(c, mm71(:, 2));

pressure = [p01, p11, p21, p31, p41, p51, p61, p71]; 
zNames = {'1 мм'; '11 мм'; '21 мм'; '31 мм';'41 мм'; '51 мм'; '61 мм'; '71 мм'};
z = [1 11 21 31 41 51 61 71];

f2 = figure();
hold on;
grid on;
title({ 'Сечения затопленной струи', 'на разном расстоянии от сопла'});
ylabel('\DeltaP, Пa');
xlabel( 'Расстояние вдоль сечения струи, мм');

for i = 1:size(pressure, 2)
    plot(x, pressure(:, i), 'DisplayName', zNames{i});
end

legend('Location', 'NorthWest');
saveas( f2,'pressure.png');



f3 = figure();

hold on;
grid on;
title({ 'Центрированые сечения затопленной струи', 'на разном расстоянии от сопла'});   
ylabel('\DeltaP, Пa');
xlabel({'Pасстояние вдоль сечения струи', 'относительно её центра, мм'});

xCentered = zeros(size(pressure));
offset = 50;

for i = 1:size(pressure, 2)
    right = x(find(pressure(:, i) > offset, 1, 'last'));
    left = x(find(pressure(:, i) > offset, 1));

    center = left + (right - left) / 2;
    xCentered(:, i) = x - center;

    plot(xCentered(:, i), pressure(:, i), 'DisplayName', zNames{i});
end

legend ('Location', 'NorthWest');
saveas (f3, 'centered.png');





f4 = figure();

sum = 0;
for i = 1:length(adc0)
    sum = sum + adc0(i);
end
sum = sum/length(adc0);
Pst = polyval(c, sum);

speed = pressure;
for i = 1:size(pressure, 2)
     speed(:, i) = sqrt(abs(pressure(:, i) - Pst)/(0.6));
end

zn = speed;
for i=1:8
    zn(:,i) = 10*(i-1)+1
end

grid on;

xCentered = zeros(size(speed));
offset = 10;



for i = 1:size(speed, 2)
    right = x(find(speed(:, i) > offset, 1, 'last'));
    left = x(find(speed(:, i) > offset, 1));

    center = left + (right - left) / 2;
    xCentered(:, i) = x - center;
end
surf(zn, xCentered, speed)
axis([ 0 80 -20 20 0 30])
title({ 'Распределение скоростей потока', 'в затопленной струе'});   
ylabel('x, мм');
xlabel({'z, мм'});
zlabel({'Скорость потока, м/с'})

colorbar;
saveas (f4, 'speed.png');



f5 = figure();

hold on;
grid on;
title({ 'Расход затопленной струи'});
ylabel('Q, м^3/c');
xlabel( 'z, мм');

Q =  zeros(8);

for i = 1:8
    y0 = speed(:, i).*abs(xCentered(:, i));
    x0 = xCentered(:, i).*10^(-3);
    Q(i) = 2*1.2*pi*trapz(x0,y0);
end

Q = Q(:,1)'/2;

plot(z, Q, 'LineWidth', 1.5);
saveas( f5,'q.png');