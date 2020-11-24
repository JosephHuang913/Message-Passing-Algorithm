close all; 
clear all;

load s_random.log;

for i=0:1:408
    tone(i+1,:) = floor(s_random(i*20+1:20*i+20, 2)/8);
end

l=0;
for i=1:1:409
    for j=1:1:20
        for k=1:1:20
            if (tone(i,j) == tone(i,k) && j ~= k)
                l=l+1;
                error(l,1) = i;
                error(l,2) = j;
                error(l,3) = k;
                error(l,4) = tone(i,j);
            end
        end
    end
end
