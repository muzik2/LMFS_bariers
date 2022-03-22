function [b] = IsInCircle(xr,yr,r,x,y)
    arguments
      xr (1,1) {mustBeNumeric,mustBeFinite} = 0;
      yr (1,1) {mustBeNumeric,mustBeFinite} = 0;
      r  (1,1) {mustBeNumeric,mustBeFinite} = 1;
      x  (1,:) {mustBeNumeric,mustBeFinite} = [0.5,0,1];
      y  (1,:) {mustBeNumeric,mustBeFinite} = [0.5,0.5,1];
    end
    b = ((xr-x).^2 + (yr-y).^2)<=r*r;
end