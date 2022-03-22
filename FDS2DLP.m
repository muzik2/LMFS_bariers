function [h,hx,hy] = FDS2DLP(x,y,vx,vy)
    rx = x-vx;
    ry = y-vy;
    r = sqrt(rx^2+ry^2);
    rx = rx/r; % dr/dx
    ry = ry/r; % dr/dy
    h = log(1.0/r)/(2*pi); % u*
    hx = -(1.0/2/pi/r)*rx; % du*/dx
    hy = -(1.0/2/pi/r)*ry; % du*/dy
end