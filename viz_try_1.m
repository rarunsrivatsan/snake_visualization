
%%
hold on 
state=[0,0,0,0,0,0,0,0];
for ii=1:20
state=[state,rand*pi/180,rand*pi/180];
drawColor=[0.2 0.1 0.3];
LINK_LENGTH=7.9;
LINK_RADIUS=3;
drawType=1;
Tregister=eye(4);
linkStartDraw=0;
axis([0 150 -50 50 -50 50]);
[h, snakePoints] = drawState(state,drawColor,LINK_LENGTH,LINK_RADIUS,drawType,Tregister,linkStartDraw);
pause(1);
end