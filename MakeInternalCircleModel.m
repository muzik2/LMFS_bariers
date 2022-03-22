function [xc,yc,xcf,ycf] = MakeInternalCircleModel(xr,yr,r,f,nr)
    arguments
        xr (1,1) {mustBeNumeric, mustBeFinite} = 0;
        yr (1,1) {mustBeNumeric, mustBeFinite} = 0;
        r  (1,1) {mustBeNumeric, mustBeFinite} = 1;
        f  (1,1) {mustBeNumeric, mustBeFinite} = 2/3;
        nr (1,1) {mustBeNumeric, mustBeFinite} = 11;
    end
    
    t = linspace(0,2*pi,nr+1);
    t = t(1:end-1);
    
    xc = xr+r*cos(t);
    yc = yr+r*sin(t);
    
    xcf = xr+(f*r)*cos(t);
    ycf = yr+(f*r)*sin(t);

    %{
    plot(xc,yc,'or',xcf,ycf,'k+');
    daspect([1,1,1]);
    axis([xr-r-0.1,xr+r+0.1,yr-r-0.1,yr+r+0.1]);
    %}
end