x = 0:10;
y = sin(x);
xx = 0:.25:10;
yy = spline(x,y,xx);
figure
plot(x,y,'o',xx,yy)

yy1 = [0.0
0.24605683812186901
0.4778380394757801
0.6810679672937752
0.8414709848078965
0.9469588960370484
0.9941932681475852
0.9820231090927241
0.9092974268256817
0.7772292977896724
0.5964880723878996
0.38010716951356427
0.1411200080598672
-0.10707405010535628
-0.3496118712157342
-0.5712643785302602
-0.7568024953079282
-0.8929611253145011
-0.9743310943428202
-0.9974672086924958
-0.9589242746631385
-0.8577626816395358
-0.703065151347184
-0.5064199885967562
-0.27941549819892586
-0.03431878254146248
0.21388786567947882
0.44954335619064656
0.6569865987187891
0.8220868617888941
0.9368348491189058
0.9947516232250073
0.9893582466233818
0.9174103756080012
0.7886020415839928
0.6158618697342726
0.4121184852417566
0.18889885023636985
-0.04787672536392673
-0.2936895946941633
-0.5440211108893698]';

%p = polyfit(x,y,7);
%y1 = polyval(p,xx);
hold on
plot(xx,yy1,'-sb')