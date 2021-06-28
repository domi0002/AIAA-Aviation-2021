clear all
clc

% Kinematics from Roh and Gharib[2019]
% Digitized and simulated here
% D. CHANDAR- QUB 2021.

% Time   yLE   Angle(deg)

mydata2=[0.000000000000	0.000000000000		7.134831460674
0.000169902913	0.001500000000		7.347623723332
0.000412621359	0.002666666667		7.651612668554
0.000606796117	0.003833333333		7.894803826235
0.000825242718	0.005000000000		8.176445757343
0.000995145631	0.006166666667		8.401959296258
0.001237864078	0.007333333333		8.724121494518
0.001432038835	0.008666666667		8.998605584181
0.001626213592	0.009666666667		9.276160708483
0.001820388350	0.011000000000		9.553152613413
0.002014563107	0.011666666667		9.781144322247
0.002160194175	0.012833333333		9.952138104166
0.002378640777	0.014166666667		10.208628777044
0.002621359223	0.015500000000		10.578706227970
0.002791262136	0.016666666667		10.857996072903
0.003009708738	0.017666666667		11.186199290071
0.003131067961	0.018333333333		11.344526865733
0.003276699029	0.019166666667		11.534519957048
0.003446601942	0.020666666667		11.756178564017
0.003616504854	0.022166666667		12.002759899666
0.003762135922	0.023333333333		12.262670448183
0.003932038835	0.024500000000		12.565899422048
0.004029126214	0.025666666667		12.739173121654
0.004174757282	0.027000000000		12.996866230130
0.004320388350	0.028500000000		13.230039569553
0.004514563107	0.029833333333		13.540937354917
0.004660194175	0.031333333333		13.774110694340
0.004805825243	0.032500000000		14.028733501449
0.004927184466	0.033666666667		14.245325624618
0.005097087379	0.035000000000		14.548554598482
0.005242718447	0.036000000000		14.808465146999
0.005388349515	0.037333333333		15.049783314050
0.005533980583	0.038666666667		15.282956653474
0.005655339806	0.040166666667		15.477267769126
0.005873786408	0.041166666667		15.836696848020
0.005995145631	0.042833333333		16.050439074832
0.006165048544	0.044500000000		16.349678193777
0.006359223301	0.045500000000		16.691665757028
0.006504854369	0.046833333333		16.971110141394
0.006626213592	0.048333333333		17.208601505345
0.006844660194	0.049000000000		17.355447944180
0.007038834951	0.049500000000		17.414924042136
0.007305825243	0.050333333333		17.526219825838
0.007548543689	0.050666666667		17.689071046509
0.007839805825	0.051166666667		17.884492511850
0.008033980583	0.051666666667		17.956382755076
0.008252427184	0.052166666667		18.011345041952
0.008495145631	0.052333333333		18.072414249955
0.008762135922	0.053000000000		18.213592232952
0.009004854369	0.053500000000		18.365586705978
0.009271844660	0.054166666667		18.532780625868
0.009538834951	0.054666666667		18.679893458071
0.009733009709	0.055166666667		18.786289589556
0.009975728155	0.055500000000		18.919284753090
0.010291262136	0.055666666667		18.839314933893
0.010606796117	0.055666666667		18.710046906789
0.010873786408	0.055666666667		18.530910564050
0.011189320388	0.055666666667		18.319203977360
0.011383495146	0.055666666667		18.230011915521
0.011674757282	0.055666666667		18.098478236911
0.011966019417	0.055500000000		17.978160794319
0.012208737864	0.055500000000		17.886964110504
0.012475728155	0.055500000000		17.772243689169
0.012742718447	0.055666666667		17.607249688297
0.013033980583	0.055500000000		17.421654363795
0.013300970874	0.055666666667		17.247494030577
0.013567961165	0.055333333333		16.934247845530
0.013859223301	0.054500000000		16.585150150175
0.014174757282	0.054000000000		16.197021406679
0.014417475728	0.053333333333		15.898460835447
0.014660194175	0.052500000000		15.599900262838
0.014927184466	0.052166666667		15.321565846752
0.015291262136	0.051500000000		14.897614112617
0.015533980583	0.050833333333		14.544769799534
0.015898058252	0.050166666667		14.122992118803
0.016286407767	0.049000000000		13.625668157912
0.016699029126	0.048333333333		13.093051161832
0.017014563107	0.047333333333		12.867996187582
0.017451456311	0.046000000000		12.534938123239
0.017864077670	0.045166666667		12.174567389246
0.018155339806	0.044333333333		11.895910855748
0.018398058252	0.043666666667		11.663697078471
0.018713592233	0.042833333333		11.449596638203
0.019126213592	0.040666666667		11.128068070702
0.019563106796	0.040000000000		10.664190723884
0.020024271845	0.039166666667		10.268770130574
0.020291262136	0.037500000000		10.037778530726
0.020582524272	0.035166666667		9.600264535480
0.020849514563	0.033833333333		9.241325108792
0.021116504854	0.032000000000		8.956335472617
0.021456310680	0.030166666667		8.554343477763
0.021747572816	0.027833333333		8.200883603389
0.021990291262	0.025833333333		7.915893967640
0.022305825243	0.024166666667		7.522362822667
0.022597087379	0.022333333333		7.138485872167
0.022864077670	0.020333333333		6.762299552415
0.023179611650	0.018166666667		6.433879677493
0.023373786408	0.016833333333		6.319883822489
0.023834951456	0.015833333333		5.985198624530
0.024223300971	0.014333333333		5.852991709421
0.024660194175	0.012333333333		5.674206888354
0.025048543689	0.011000000000		5.455566321907
0.025485436893	0.009833333333		5.195497023408
0.025873786408	0.008666666667		5.084955588397
0.026189320388	0.006833333333		4.992290490367
0.026601941748	0.005500000000		4.805951112035
0.027111650485	0.004666666667		4.636285953382
0.027475728155	0.004166666667		4.550561797753
0.027936893204	0.003500000000		4.546993745647
0.028495145631	0.002500000000		4.438014166301
0.029126213592	0.001833333333		4.386912576256
0.029733009709	0.000833333333		4.382022471910
0.030242718447	0.000333333333		4.687847714738
0.030970873786	0.000000000000		5.091473946143
0.031432038835	0.000000000000		5.410153266874
0.031868932039	0.000000000000		5.769240209655
0.032305825243	0.000000000000		5.963927514933
0.032694174757	0.000000000000		6.339443271584
0.033252427184	0.000000000000		6.778117158814
];


mydata=[mydata2;mydata2];

% Initial points on the plate (Used for visualization only)
x(1)=0;
y(1)=0;

x(2)=0.635;
y(2)=0;

x(3)=0.635;
y(3)=0.01;

x(4)=0;
y(4)=0.01;

x(5)=x(1);
y(5)=y(1);

plot(x,y,'r');
axis ([-1  1 -0.5 0.5])
axis equal

tt= mydata(:,1);
yy= mydata(:,2);
aa= mydata(:,3);

% Initial rotation:
rmat=[cosd(aa(1)) -sind(aa(1)); sind(aa(1)) cosd(aa(1))];
for i = 1:5
    vec=rmat*[x(i);y(i)];
    x(i)=vec(1);
    y(i)=vec(2);
end


plot(x,y,'r');
axis ([0.6  0.635 0.08 0.096])
axis equal


T=1;
dt=0.001;

% Stroke angle
sAngle=35;

% Frequency (Hz)
ff=30;
 
for t=1:size(mydata,1)
   xnew=x-yy(t)/tand(sAngle);
   ynew=y+yy(t);
   
   % Leading edge.
   xrot=xnew(2);
   yrot=ynew(2);
   
   % Add rotation about xrot  and yrot(so that the LE moves along the
   % stroke plane
   
   theta=aa(t);
   rmat=[cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    for i = 1:5
        vec=rmat*[xnew(i)-xrot;ynew(i)-yrot];
        xnew(i)=xrot+vec(1);
        ynew(i)=yrot+vec(2);
    end
   
   
    
   plot(xnew,ynew,'r');
   
   axis equal
   axis ([-0.2  0.85 -0.15 0.2])
   pause(0.001);
   
end
    
