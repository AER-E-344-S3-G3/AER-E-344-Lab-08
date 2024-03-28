close,clc,clear;
%%Kyle Ostendorf Lab 8
%% counting varibles
X = [0:10,15:5:65];
goodTAP = [3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50,55,60,65];
y = 4:4:4*38;
rho = 1.2;

IN0 = readmatrix("0in.csv");

%% Reading Cell array 
for i = 1:length(X)

MAREIX = readmatrix(sprintf("%din.csv",X(i)));

TotalPressure(i,:) = mean([MAREIX(:,2:17),MAREIX(:,35:50),MAREIX(:,68:73)],"omitmissing");

Staticpressure(i,:) = mean(MAREIX(:,74));
% figure(i)
% bar(TotalPressure(i,:))

end

%% subtracting the total - static

BadTAP = [1,2,21,22,35];  

for i =1:22
goodTAP = 1:38;
goodTAP(BadTAP) = [];
goodData = TotalPressure(i,goodTAP);
interpdata(i,:) = interp1(goodTAP,goodData,BadTAP,"linear","extrap");
TotalPressure1(i,:) = [interpdata(i,1:2),goodData([1:18]),interpdata(i,[3:4]),goodData([19:30]),interpdata(i,[5]),goodData([31:33])];
end
% for i = 1:length(X)
% figure(i)
% bar(TotalPressure1(i,:))
% end

Q = TotalPressure1 - Staticpressure;

V = sqrt(2*Q/rho);

for i = 1:length(X)
figure(i)
plot(V(:,i))
end


