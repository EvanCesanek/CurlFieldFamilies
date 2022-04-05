function wl_draw_rectangle_gl(pos,width,height,rotate,col,pivotUV,varargin)
%WL_DRAW_RECTANGLE draws a sphere on the active psychtoolbox screen
%   WL_DRAW_RECTANGLE(POS,WIDTH,HEIGHT,COL) draws a sphere at POS ([X Y Z]) of
%   scalar WIDTH AND HEIGHT and color COL ([R G B])
%
%   WL_DRAW_SPHERE(POS,RADIUS,COL,ALPHA) Draws a rectangle as
%   per above with aditional customisation parameters such as the
%   transparency level ALPHA (goes from 0 to 1)
global GL

pSet = inputParser;
addParameter(pSet,'Alpha',1);
parse(pSet,varargin{:});
v2struct(pSet.Results);

glPushMatrix();
glTranslatef(pos(1), pos(2), pos(3));%-10);
glRotatef(rotate,0,0,1);
%glMaterialfv(GL.FRONT_AND_BACK,GL.AMBIENT, [col Alpha]);
%glMaterialfv(GL.FRONT_AND_BACK,GL.DIFFUSE, [col Alpha]);
glColor4fv([col Alpha]);
glRectf(-width*pivotUV(1),-height*pivotUV(2),width*(1-pivotUV(1)),height*(1-pivotUV(2)));
glPopMatrix();