function [str] = timeToString(time)

s = seconds(25.6)
s.Format = 'hh:mm:ss'


  time = round(time);
  min = floor(time / 60);
  sec = time - 60 * min;
  min = num2str(min);
  sec = num2str(sec);
  if length(sec) == 1
    sec = strcat('0', sec);
  end
  str = strcat(min, ':', sec);
end

