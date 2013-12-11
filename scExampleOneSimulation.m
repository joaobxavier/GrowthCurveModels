% initial parameters
p.muMax = 0.3555;
p.ni0   = 0.27743;
p.Yc    = 0.4906;
p.Yn    = 4.9042;
p.f     = 223.3470;
p.ii0   = 3527.45087; % initial value of internal iron
p.Yi    = 0.00666;
p.qG    = p.muMax / p.ii0;
p.kd    = 0.0046;

p.qR = 10288.342088;
p.qRb = 25164.535010;
p.ki1 = 1777.979240;
p.ki2 = 700.244794;
p.kg = 53590.000000;
p.h1 = 5.354021;
p.h2 = 1.025880;
p.qRb2 = 916.652509;

cmap = jet(4);


load('backbone');

x0  = 28e-5; % OD
gi0 = 7.32e4;
simulation = 1;

C0Value = 49;
N0Value = 0.5;
I0Value = 8.30988; % initial value of external iron


y0 = [x0, C0Value, N0Value, I0Value, p.ni0, p.ii0, gi0];

myModel = ModelWithGfp(p);
load('backbone');


myModel = myModel.setData(timehoursBackBone, odBackBone, gfpBackBone);
myModel = myModel.setInitialValues(y0);
myModel = myModel.solveModel;

%%
figure(1);
myModel.plotDataAndModel([1 0 0], 'dlogxdt');







