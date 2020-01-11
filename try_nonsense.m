t = [0:0.1:10*pi];


subplot 311
plot(t, sin(t));
subplot 312
plot(t, sin(3*t), 'r');
subplot 313
plot(t, sin(5*t), 'k');

figure
subplot 311
plot(t, sin(t)+sin(3*t));
subplot 312
plot(t, sin(t)+sin(3*t)+sin(5*t));
subplot 313
plot(t, sin(t)+sin(3*t)+sin(5*t)+sin(7*t));

answer = 0;
for loop = 1:2:1000
    answer =  answer +[ cos(loop*t)]
end
plot(answer)