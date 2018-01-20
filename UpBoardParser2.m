%% UpBoardParser
close all, clear all,clc
[file,path,~] = uigetfile('.txt','Choose a Log File');
file(end-2:end) = 'txt';
log = readtable(file,'ReadVariableNames',1,'Delimiter','\t');
%% ParseTable
index=1;
for i=1:size(log)
    if strcmp(log.PACKET{i},'ADCS')
            ADCSINST(index,2)=log.Q1(i);
            index=index+1;
            ADCSINST(index,1)=log.TIMESECONDS(i);
    end 
end
ADCS_Q1=ADCS_Q1';
%{
t = log.time/1000;
acc = [log.accelx,log.accely,log.accelz];
acc = acc * ;
gyr = [log.gyrx,log.gyry,log.gyrz];
gyr = gyr * ;
mag = [log.magx,log.magy,log.magz];
mag = mag * ;
press = log.pressure;
temp = log.temperature;
alt = log.alt;

%flight_range = log.time;
clear log;

flIndex  = size(t);
flight_range = 1:flIndex(1);
acc = acc(flight_range,:);
gyr = gyr(flight_range,:);
mag = mag(flight_range,:);
%alt = alt(flight_range);
%press = press(flight_range);
%temp = temp(flight_range);
t = t(flight_range);

f1 = figure;
    plot(t,acc.*0.000732);

    legend('x','y','z')
    
f2 = figure;
    plot(t,gyr.*0.00875);

    legend('x','y','z')
%{
f3 = figure;
    plot(t,alt);
    
    ylabel('m asl')

f4 = figure;
    plot(t,temp);

    ylabel('deg C')
%}    
f5= figure;
    plot(t,mag.*0.00029);

    legend('x','y','z')  
    
return
window = 10;
coeff = ones(window,1)./window;
a = 1;
alt_filt = filter(coeff,a,alt);
speed = (alt_filt(4001:end)-alt_filt(4000:(end-1)))./(t(4001:end)-t(4000:(end-1)));
plot(t(4001:end),speed)
%}