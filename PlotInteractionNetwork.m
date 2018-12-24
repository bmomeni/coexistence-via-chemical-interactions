function chk = PlotInteractionNetworkExMT4_StrengthBasalFitness(rb,rint,At,Bt,rth,cc)

disp('local copy')
A = At';
B = Bt';
[Ns, Nc] = size(A);

xT = 15/Ns*12*max(Ns,Nc);
yT = 15/Ns*4*max(Ns,Nc);
xs = linspace(-xT,xT,Ns);
xc = linspace(-xT,xT,Nc);
r = 120/max(Ns,Nc);
q = 180/max(Ns,Nc);
ys = yT*ones(size(xs));
yc = -yT*ones(size(xc));
s = 17;
pp = 12;
pt = 6;

normB = 3/max(max(B));
rc = rint.*A;

figure(cc)
for kk = 1:Ns
    plot(xs(kk),ys(kk),'o','MarkerSize',r,'color',(1-0.9*rb(kk)/max(rb))*[1 1 1])
    hold on
end
axis equal
plot(xc,yc,'ms','MarkerSize',q)
axis([-1.3*xT 1.3*xT -1.8*yT 1.8*yT])

for ct = 1:Ns,
    text(xs(ct)-7.5,ys(ct)+25,num2str(ct))
end
for ct = 1:Nc,
    text(xc(ct)-12,yc(ct)-25,strcat('C',num2str(ct)))
end
for i1 = 1:Ns,
    for i2 = 1:Nc,

        dx = (xc(i2)-xs(i1))/sqrt((xs(i1)-xc(i2))^2+(ys(i1)-yc(i2))^2);
        dy = (yc(i2)-ys(i1))/sqrt((xs(i1)-xc(i2))^2+(ys(i1)-yc(i2))^2);

        x1 = xs(i1)+s*dx; y1 = ys(i1)+s*dy;
        x2 = xc(i2)-s*dx; y2 = yc(i2)-s*dy;
        cl = [0 0 0];
        if rc(i1,i2)>rth,
            plot([x1 x2],[y1 y2],'color',cl,'linewidth',0.5)
            plot([x1 x2-0.5*(x2-x1)],[y1 y2-0.5*(y2-y1)],'color',cl,'linewidth',abs(rc(i1,i2))*20)
            plot([x1-pt*dy+pp*dx x1 x1+pt*dy+pp*dx],[y1+pt*dx+pp*dy y1 y1-pt*dx+pp*dy],'color',cl,'linewidth',abs(rc(i1,i2))*20)
        end
        if rc(i1,i2)<-rth,
            plot([x1 x2],[y1 y2],'color',cl,'linewidth',0.5)
            plot([x1 x2-0.5*(x2-x1)],[y1 y2-0.5*(y2-y1)],'color',cl,'linewidth',abs(rc(i1,i2))*20)
            plot([x1-pt*dy x1+pt*dy],[y1+pt*dx y1-pt*dx],'color',cl,'linewidth',abs(rc(i1,i2))*20)
        end
    end
end

for i1 = 1:Nc,
    for i2 = 1:Ns,

        dx = (xs(i2)-xc(i1))/sqrt((xc(i1)-xs(i2))^2+(yc(i1)-ys(i2))^2);
        dy = (ys(i2)-yc(i1))/sqrt((xc(i1)-xs(i2))^2+(yc(i1)-ys(i2))^2);

        x1 = xc(i1)+s*dx; y1 = yc(i1)+s*dy;
        x2 = xs(i2)-s*dx; y2 = ys(i2)-s*dy;
        cl = [0 0 0];
        if B(i2,i1)>0,
            plot([x1 x2],[y1 y2],'color',cl,'linewidth',0.5)
            plot([x1 x2-0.5*(x2-x1)],[y1 y2-0.5*(y2-y1)],'color',cl,'linewidth',abs(B(i2,i1))*normB)
            plot([x1-pt*dy+pp*dx x1 x1+pt*dy+pp*dx],[y1+pt*dx+pp*dy y1 y1-pt*dx+pp*dy],'color',cl,'linewidth',abs(B(i2,i1))*normB)
        end
    end
end
chk = 0;

return

